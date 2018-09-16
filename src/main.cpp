#include <stdio.h>
#include <cstdlib>
#include <time.h>
#include <math.h>
#include <iostream>
#include <tiledb/tiledb>
#include <arrayfire.h>
#include <af/util.h>

using namespace af;
using namespace tiledb;

// Name of array.
std::string skeleton("skeleton");
std::string dilated("dilated");

void create_array() {
  // Create a TileDB context.
  Context ctx;

  // The array will be 4x4 with dimensions "rows" and "cols", with domain [1,4].
  Domain domain(ctx);
  domain.add_dimension(Dimension::create<int>(ctx, "rows", {{0, 3}}, 4))
      .add_dimension(Dimension::create<int>(ctx, "cols", {{0, 3}}, 4))
      .add_dimension(Dimension::create<int>(ctx, "pages", {{0, 3}}, 4));

  // The array will be sparse.
  ArraySchema schema(ctx, TILEDB_SPARSE);
  schema.set_domain(domain).set_order({{TILEDB_ROW_MAJOR, TILEDB_ROW_MAJOR}});

  // Add a single attribute "a" so each (i,j) cell can store an integer.
  schema.add_attribute(Attribute::create<int>(ctx, "voxel"));

  // Create the (empty) array on disk.
  Array::create(skeleton, schema);
  Array::create(dilated, schema);
}

void write_array() {
  Context ctx;

  // Write some simple data to cells (1, 1), (2, 4) and (2, 3).
  std::vector<int> coords = {0, 0, 0, 1, 3, 0, 1, 2, 0};
  std::vector<int> sk_data = {1, 4, 3};
  std::vector<int> di_data = {10, 40, 30};

  // Open the array for writing and create the query.
  Array sk_array(ctx, skeleton, TILEDB_WRITE);
  Query sk_query(ctx, sk_array, TILEDB_WRITE);
  sk_query.set_layout(TILEDB_UNORDERED)
      .set_buffer("voxel", sk_data)
      .set_coordinates(coords);

  // Perform the write and close the array.
  sk_query.submit();
  sk_array.close();

  // Open the array for writing and create the query.
  Array di_array(ctx, dilated, TILEDB_WRITE);
  Query di_query(ctx, di_array, TILEDB_WRITE);
  di_query.set_layout(TILEDB_UNORDERED)
      .set_buffer("voxel", di_data)
      .set_coordinates(coords);

  // Perform the write and close the array.
  di_query.submit();
  di_array.close();
}

void read_array(std::string array_name) {
  Context ctx;

  // Prepare the array for reading
  Array array(ctx, array_name, TILEDB_READ);

  // Slice only rows 1, 2 and cols 2, 3, 4
  const std::vector<int> subarray = {0, 1, 1, 3, 0, 1};

  // Prepare the vector that will hold the result.
  // We take an upper bound on the result size, as we do not
  // know a priori how big it is (since the array is sparse)
  auto max_el = array.max_buffer_elements(subarray);
  std::vector<int> data(max_el["voxel"].second);
  std::vector<int> coords(max_el[TILEDB_COORDS].second);

  // Prepare the query
  Query query(ctx, array, TILEDB_READ);
  query.set_subarray(subarray)
      .set_layout(TILEDB_ROW_MAJOR)
      .set_buffer("voxel", data)
      .set_coordinates(coords);

  // Submit the query and close the array.
  query.submit();
  array.close();

  // Print out the results.
  auto result_num = (int)query.result_buffer_elements()["voxel"].second;
  for (int r = 0; r < result_num; r++) {
    int i = coords[3 * r], j = coords[3 * r + 1], k = coords[3 * r + 2];
    int a = data[r];
    std::cout << "Cell (" << i << ", " << j << ", " << k << ") has data " << a << "\n";
  }
}

// creates a 3D spherical mask of given radius (number of voxels)
array create_mask(int radius)
{
    int sz = 2*radius + 1;
    int c = radius; //center coordinate; x = y = z = c;
    int sqrad = (radius-1)*(radius-1), sqdist;
    
    array mask = constant(0, sz, sz, sz, u8);
    for (int x = 0; x < sz; x++)
        for (int y = 0; y < sz; y++)
            for (int z = 0; z < sz; z++)
            {
                sqdist = (x-c)*(x-c) + (y-c)*(y-c) + (z-c)*(z-c);
                if (sqdist <= sqrad) mask(x, y, z) = 1;
            }
    return mask;
}
/*
array create_mask(int radius)
{
    int sz = 2*radius + 1;
    int c = radius; //center coordinate; x = y = z = c;
    array mask = constant(0, sz, sz, sz, u8);
    int R2 = radius*radius;
    for (int z = 0; z <= c; z++)
    {
        int r2 = R2 - (radius-z)*(radius-z);
        for (int y = 0; y <= c; y++)
        {
            int l2 = r2 - (radius-y)*(radius-y);
            if (l2 < 0) continue;
            int l = (int)round(std::sqrt(l2));
            //std::cout << y << " " << z << " " << l << std::endl;
            if (l >= 0 && l <= c) mask(seq(c-l, c+l), y, z) = 1;
            if (l >= 0 && l <= c) mask(seq(c-l, c+l), sz-y-1, z) = 1;
            if (l >= 0 && l <= c) mask(seq(c-l, c+l), y, sz-z-1) = 1;
            if (l >= 0 && l <= c) mask(seq(c-l, c+l), sz-y-1, sz-z-1) = 1;
        }
    }
    return mask;
}
*/

// dilates 3D array block using the 3D mask
array dilate_with_mask(array block, array mask, bool useFFT = true, bool usePad = true)
{
    array out;
    if (useFFT)
    {
        if (usePad)
            out = fftConvolve3(block, mask, AF_CONV_EXPAND); // with padding, expand
        else
            out = fftConvolve3(block, mask); // without padding, same
        out = out >= 1;
    }
    else
        out = dilate3(block, mask); // can't use padding
        
    return out;
}


static void dilate_test()
{
    int block_size = 256;
    int mask_radius = 15;
    
    array block = constant(0, block_size, block_size, block_size, u8);
    
    array mask = create_mask(mask_radius);
    //af_print(mask);
    
    // create a diagonal line for testing
    for (int i = 3; i < block_size-3; i++)
        block(i, i, i) = 1;
    
    clock_t s1, e1, s2, e2;
    int N = 1;
    
    array di;
    s1 = clock();
    for (int i = 0; i < N; i++)
    {
        di = dilate_with_mask(block, mask, false);
    }
    e1 = clock();
    double elapsed1 = (double)e1 - (double)s1;
    
    array fc;
    s2 = clock();
    for (int i = 0; i < N; i++)
    {
        fc = dilate_with_mask(block, mask);
    }
    e2 = clock();
    double elapsed2 = (double)e2 - (double)s2;
    
    //for (int i = 0; i < block_size; i += 64)
    //    saveImage("testf.png", fc.slice(i));
        
    //saveImage("testf050.png", fc.slice(50));
    
    std::cout << "Dilate3 time = " << elapsed1/N << " ms." << std::endl;
    std::cout << "FFT convolve3 time = " << elapsed2/N << " ms." << std::endl;
    
    //std::cout << di.dims() << std::endl;
    //std::cout << fc.dims() << std::endl;
}

int main(int argc, char** argv)
{
    int device = argc > 1 ? atoi(argv[1]) : 0;

    try {
        af::info();
        af::setDevice(device);
        printf("** Dilate using ArrayFire **\n\n");
        dilate_test();

    } catch (af::exception& e) {
        fprintf(stderr, "%s\n", e.what());
        throw;
    }
    

    Context ctx;

    if (Object::object(ctx, skeleton).type() != Object::Type::Array) {
        create_array();
    }
    
    write_array();

    read_array(skeleton);
    read_array(dilated);

    return 0;
}
