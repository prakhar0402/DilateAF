#include <foam.h>

Foam::Foam()
{
    attri = "voxel";
    skeleton = "Skeleton";
    dilated = "Dilated";
    
    voxel_dim = 256;
    tile_dim = 64;
    offset = 0;
    
    float tmp[] = {0.0, 0.0, 0.0, 10.0, 10.0, 10.0};
    bb.assign(tmp, tmp+6);

    Context ctx;

    if (Object::object(ctx, skeleton).type() != Object::Type::Array)
    {
        create_tiledb_array(skeleton, voxel_dim, tile_dim);
        create_tiledb_array(dilated, voxel_dim, tile_dim);
    }
}

// dilates 3D array block using the 3D mask
array Foam::dilate_block(const array& block, const array& mask, bool usePad = true)
{
    array out;
    if (usePad)
        out = fftConvolve3(block, mask, AF_CONV_EXPAND); // with padding, expand
    else
        out = fftConvolve3(block, mask); // without padding, same
    out = out >= 1;
    return out;
}

void Foam::dilate(const array& mask)
{
    std::vector<char> sk_data, di_data;
    std::vector<int> sk_coords, di_coords;
    array block, di_block;
    std::vector<int> tile_coords(6, 0);
    for (int i = 0; i < voxel_dim; i += tile_dim)
        for (int j = 0; j < voxel_dim; j += tile_dim)
            for (int k = 0; k < voxel_dim; k += tile_dim)
            {
                std::cout << "(" << i << ", " << j << ", " << k << ")" << std::endl;
                tile_coords = {i, i+tile_dim-1, j, j+tile_dim-1, k, k+tile_dim-1};
                read_tiledb_array_relative(skeleton, tile_coords, sk_data, sk_coords);
                
                if (sk_coords.size() == 0) continue;
                
                std::cout << sk_coords.size()/3 << std::endl;
                block = create_af_array(tile_dim, sk_coords);
                di_block = dilate_block(block, mask, true);
                
                get_true_coords(di_block, di_coords, tile_coords, mask.dims());
                di_data.resize(di_coords.size()/3, '1');
                write_tiledb_array(dilated, di_data, di_coords);
            }
}


array Foam::create_af_array(const int dim, const std::vector<int>& coords)
{
    array block = constant(0, dim, dim, dim, b8);
    for (int i = 0; i < coords.size()/3; i++)
        block(coords[3*i], coords[3*i+1], coords[3*i+2], 0) = true;
    return block;
}

void Foam::get_true_coords(const array& block, std::vector<int>& coords, std::vector<int>& tile_coords, dim4 mask_dims)
{
    dim4 dims = block.dims();
    array nz = where(block);
    int idx, x, y, z;
    for (int i = 0; i < nz.dims()[0]; i++)
    {
        //TODO: check if linear index to array index is correct
        idx = sum<int>(nz(i));
        z =  idx % dims[0] + tile_coords[4] - (mask_dims[2]-1)/2;
        idx /= dims[0];
        y = idx % dims[1] + tile_coords[2] - (mask_dims[1]-1)/2;
        idx /= dims[1];
        x = idx + tile_coords[0] - (mask_dims[0]-1)/2;
        if (x < 0 || x >=voxel_dim || y < 0 || y >= voxel_dim
                || z < 0 || z >= voxel_dim) continue;
        //std::cout << "(" << x << "," << y << "," << z << ") ";
        coords.insert(coords.end(), {x, y, z});
    }
    std::cout << std::endl;
}
/*
void Foam::get_true_coords(const array& block, std::vector<int>& coords)
{
    dim4 dims = block.dims();
    for (int i = 0; i < dims[0]; i++)
        for (int j = 0; j < dims[1]; j++)
            for (int k = 0; k < dims[2]; k++)
                if (block(i, j, k).scalar<bool>())
                //if (sum<int>(block(i, j, k)) > 0)
                    coords.insert(coords.end(), {i, j, k});
}*/

void Foam::spatial_to_voxel(const std::vector<float>& s_coords, std::vector<int>& v_coords)
{
    for (int i = 0; i < 3; i++)
        v_coords[i] = offset + (int) round((voxel_dim-2*offset)*
                                (s_coords[i] - bb[i])/(bb[i+3] - bb[i]));
}

void Foam::voxel_to_spatial(const std::vector<int>& v_coords, std::vector<float>& s_coords)
{
    for (int i = 0; i < 3; i++)
        s_coords[i] = bb[i] + (bb[i+3] - bb[i])*(v_coords[i]-offset)/(voxel_dim-2*offset);
}

// creates TileDB array, arguments specify numbers per dimension
void Foam::create_tiledb_array(const std::string array_name, const int voxel_dim, const int tile_dim)
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


void Foam::write_tiledb_array(const std::string array_name, std::vector<char>& data,
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

void Foam::read_tiledb_array_relative(const std::string array_name, const std::vector<int> slice_coords,
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
    for (int r = 0; r < result_num; r++) {
        coords[3*r] -= slice_coords[0];
        coords[3*r+1] -= slice_coords[2];
        coords[3*r+2] -= slice_coords[4];
    }
}

void Foam::read_tiledb_array(const std::string array_name, const std::vector<int> slice_coords,
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
