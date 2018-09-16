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

std::string attri("voxel"); // attribute name

// creates TileDB array, arguments specify numbers per dimension
void create_tiledb_array(const std::string array_name, const int voxel_dim, const int tile_dim)
{
    // Create a TileDB context.
    Context ctx;

    // The array will be 4x4 with dimensions "rows" and "cols", with domain [1,4].
    Domain domain(ctx);
    domain.add_dimension(Dimension::create<int>(ctx, "rows", {{0, voxel_dim-1}}, tile_dim))
        .add_dimension(Dimension::create<int>(ctx, "cols", {{0, voxel_dim-1}}, tile_dim))
        .add_dimension(Dimension::create<int>(ctx, "pages", {{0, voxel_dim-1}}, tile_dim));

    // The array will be sparse.
    ArraySchema schema(ctx, TILEDB_SPARSE);
    schema.set_domain(domain).set_order({{TILEDB_ROW_MAJOR, TILEDB_ROW_MAJOR}});

    // Add a single attribute "a" so each (i,j) cell can store an integer.
    schema.add_attribute(Attribute::create<char>(ctx, attri));

    // Create the (empty) array on disk.
    Array::create(array_name, schema);
}

void write_tiledb_array(const std::string array_name, std::vector<char>& data,
                                std::vector<int>& coords)
{
    // Create a TileDB context.
    Context ctx;

    // Open the array for writing and create the query.
    Array array(ctx, array_name, TILEDB_WRITE);
    Query query(ctx, array, TILEDB_WRITE);
    query.set_layout(TILEDB_UNORDERED)
        .set_buffer(attri, data)
        .set_coordinates(coords);

    // Perform the write and close the array.
    query.submit();
    array.close();
}

void read_tiledb_array(const std::string array_name, const std::vector<int> slice_coords,
                        std::vector<char>& data, std::vector<int>& coords)
{
    // Create a TileDB context.
    Context ctx;

    // Prepare the array for reading
    Array array(ctx, array_name, TILEDB_READ);

    // Prepare the vector that will hold the result.
    // We take an upper bound on the result size, as we do not
    // know a priori how big it is (since the array is sparse)
    auto max_el = array.max_buffer_elements(slice_coords);
    data.resize(max_el[attri].second);
    coords.resize(max_el[TILEDB_COORDS].second);

    // Prepare the query
    Query query(ctx, array, TILEDB_READ);
    query.set_subarray(slice_coords)
        .set_layout(TILEDB_ROW_MAJOR)
        .set_buffer(attri, data)
        .set_coordinates(coords);

    // Submit the query and close the array.
    query.submit();
    array.close();

    // Print out the results.
    std::cout << array_name << " data:" << std::endl;
    auto result_num = (int)query.result_buffer_elements()[attri].second;
    std::cout << "Number of entries = " << result_num << std::endl;
    data.resize(result_num);
    coords.resize(result_num*3);
    //for (int r = 0; r < result_num; r++) {
    //    std::cout << "(" << coords[3*r] << ", " << coords[3*r+1] << ", "
    //         << coords[3*r+2] << ") : " << data[r] << std::endl;
    //}
}

// creates a 3D spherical mask of given radius (number of voxels)
array create_mask(const int radius)
{
    int sz = 2*radius + 1;
    int c = radius; //center coordinate; x = y = z = c;
    int sqrad = (radius-1)*(radius-1), sqdist;
    
    array mask = constant(0, sz, sz, sz, b8);
    for (int x = 0; x < sz; x++)
        for (int y = 0; y < sz; y++)
            for (int z = 0; z < sz; z++)
            {
                sqdist = (x-c)*(x-c) + (y-c)*(y-c) + (z-c)*(z-c);
                if (sqdist <= sqrad) mask(x, y, z) = true;
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

array create_af_array(const int dim, const std::vector<int>& coords)
{
    array block = constant(0, dim, dim, dim, b8);
    for (int i = 0; i < coords.size()/3; i++)
        block(coords[3*i], coords[3*i+1], coords[3*i+2], 0) = true;
    return block;
}

void get_true_coords(const array& block, std::vector<int>& coords, int offset = 0)
{
    dim4 dims = block.dims();
    array nz = where(block);
    int idx, x, y, z;
    for (int i = 0; i < nz.dims()[0]; i++)
    {
        //TODO: check if linear index to array index is correct
        idx = sum<int>(nz(i));
        x =  idx % dims[0];
        idx /= dims[0];
        y = idx % dims[1];
        idx /= dims[1];
        z = idx;
        coords.insert(coords.end(), {x+offset, y+offset, z+offset});
    }
}
/*
void get_true_coords(const array& block, std::vector<int>& coords, int offset = 0)
{
    dim4 dims = block.dims();
    for (int i = 0; i < dims[0]; i++)
        for (int j = 0; j < dims[1]; j++)
            for (int k = 0; k < dims[2]; k++)
                if (block(i, j, k).scalar<bool>())
                //if (sum<int>(block(i, j, k)) > 0)
                    coords.insert(coords.end(), {i+offset, j+offset, k+offset});
}*/

// dilates 3D array block using the 3D mask
array dilate_with_mask(const array& block, const array& mask, bool usePad = true)
{
    array out;
    if (usePad)
        out = fftConvolve3(block, mask, AF_CONV_EXPAND); // with padding, expand
    else
        out = fftConvolve3(block, mask); // without padding, same
    out = out >= 1;
    return out;
}


static void dilate_test()
{
    int block_size = 51;
    int mask_radius = 3;
    
    array block = constant(0, block_size, block_size, block_size, b8);
    
    array mask = create_mask(mask_radius);
    //af_print(mask);
    
    // create a diagonal line for testing
    for (int i = 3; i < block_size-3; i++)
        block(i, i, i) = true;
    
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
        
    //saveImage("testf050.png", fc.slice(25));
    
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
        printf("** Dilate using ArrayFire and TileDB **\n");
        //dilate_test();

    } catch (af::exception& e) {
        fprintf(stderr, "%s\n", e.what());
        throw;
    }
    
    int voxel_dim = 256;
    int tile_dim = 64;
    
    // Name of array.
    std::string skeleton("Skeleton");
    std::string dilated("Dilated");

    Context ctx;

    if (Object::object(ctx, skeleton).type() != Object::Type::Array) {
        create_tiledb_array(skeleton, voxel_dim, tile_dim);
        create_tiledb_array(dilated, voxel_dim, tile_dim);
    }
    
    std::vector<int> coords;
    for (int i = 20; i < 236; i++)
        coords.insert(coords.end(), {i, i, i});
    std::vector<char> data(coords.size()/3, '1');
    write_tiledb_array(skeleton, data, coords);

    std::vector<char> sk_data, di_data;
    std::vector<int> sk_coords, di_coords;
    std::vector<int> slice_coords = {0, tile_dim-1, 0, tile_dim-1, 0, tile_dim-1};
    read_tiledb_array(skeleton, slice_coords, sk_data, sk_coords);
    
    array block = create_af_array(tile_dim, sk_coords);
    array mask = create_mask(10);
    array out = dilate_with_mask(block, mask, false);
    
    get_true_coords(out, di_coords, 0);
    std::vector<char> datad(di_coords.size()/3, '1');
    write_tiledb_array(dilated, datad, di_coords);
        
    read_tiledb_array(dilated, slice_coords, di_data, di_coords);

    return 0;
}
