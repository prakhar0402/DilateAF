#include <foam.h>

Foam::Foam() {
    attri = "Voxel";
    skeleton = "Skeleton";
    dilation = "Dilation";
    
    set_dims(256, 64);
    set_mask_radius(5);
    
    //float tmp[] = {0.0, 0.0, 0.0, 10.0, 10.0, 10.0};
    //bb.assign(tmp, tmp+6);
}

Foam::Foam(int v_dim, int t_dim, int rad) {
    attri = "Voxel";
    skeleton = "Skeleton";
    dilation = "Dilation";
    
    set_dims(v_dim, t_dim);
    set_mask_radius(rad);
    
    //float tmp[] = {0.0, 0.0, 0.0, 10.0, 10.0, 10.0};
    //bb.assign(tmp, tmp+6);
}

void Foam::reset() {
    reset_skeleton();
    reset_dilation();
}

void Foam::reset_skeleton() {
    Context ctx;
    if (Object::object(ctx, skeleton).type() == Object::Type::Array) {
        system(("rm -r " + skeleton).c_str());
    }
    create_skeleton(voxel_dim, tile_dim);
}

void Foam::reset_dilation() {
    Context ctx;
    if (Object::object(ctx, dilation).type() == Object::Type::Array) {
        system(("rm -r " + dilation).c_str());
    }
    create_dilation(voxel_dim, tile_dim);
}

void Foam::set_dims(int v_dim, int t_dim) {
    voxel_dim = v_dim;
    tile_dim = t_dim;
    reset();
}

void Foam::set_mask_radius(int rad) {
    radius = rad;
    create_mask();
}

void Foam::create_mask() {
    int sz = 4*radius + 1;
    int c = 2*radius; //center coordinate; x = y = z = c;
    int sqrad = (radius-1)*(radius-1), sqdist;
    
    mask = constant(0, sz, sz, sz, b8);
    for (int x = radius; x < sz-radius; x++)
        for (int y = radius; y < sz-radius; y++)
            for (int z = radius; z < sz-radius; z++) {
                sqdist = (x-c)*(x-c) + (y-c)*(y-c) + (z-c)*(z-c);
                if (sqdist <= sqrad) mask(x, y, z) = true; // set TRUE if inside sphere
            }
}
/*
array Foam::create_mask() {
    int sz = 2*radius + 1;
    int c = radius; //center coordinate; x = y = z = c;
    mask = constant(0, sz, sz, sz, u8);
    int R2 = radius*radius;
    for (int z = 0; z <= c; z++) {
        int r2 = R2 - (radius-z)*(radius-z);
        for (int y = 0; y <= c; y++) {
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
}
*/

array Foam::create_af_array(const std::vector<int>& coords,
            const std::vector<int>& tile_coords) {
    array block = constant(0, tile_coords[1] - tile_coords[0] + 1,
                                tile_coords[3] - tile_coords[2] + 1,
                                tile_coords[5] - tile_coords[4] + 1, b8);
    for (int i = 0; i < coords.size()/3; i++)
        block(coords[3*i]-tile_coords[0], coords[3*i+1]-tile_coords[2],
                    coords[3*i+2]-tile_coords[4], 0) = true;
    return block;
}

array Foam::dilate_af_array(const array& block, bool usePad) {
    array out;
    if (usePad)
        out = fftConvolve3(block, mask, AF_CONV_EXPAND); // with padding, expand
    else
        out = fftConvolve3(block, mask); // without padding, same
    out = out >= 1;
    return out;
}

void Foam::get_nz_coords(const array& block, std::vector<int>& coords,
                const std::vector<int>& tile_coords) {
    coords.clear();
    dim4 dims = block.dims();
    array nz = where(block); // finds local linear index of all true voxels
    int idx, x, y, z;
    for (int i = 0; i < nz.dims()[0]; i++) {
        // converts local linear index to global TileDB array index
        idx = sum<int>(nz(i));
        x =  idx % dims[0] + tile_coords[0] - 2*radius;
        idx /= dims[0];
        y = idx % dims[1] + tile_coords[2] - 2*radius;
        idx /= dims[1];
        z = idx + tile_coords[4] - 2*radius;
        
        // skip if index out of bounds
        if (x < 0 || x >=voxel_dim || y < 0 || y >= voxel_dim
                || z < 0 || z >= voxel_dim) continue;
        
        coords.insert(coords.end(), {x, y, z});
    }
}

void Foam::dilate() {
    std::vector<char> sk_data(1, '1'), di_data(1, '1'); // dummy initialization
    // TileDB bug: sometimes throws an error if uninitialized
    std::vector<int> sk_coords(3, 0), di_coords(3, 0);
    array block, di_block;
    std::vector<int> tile_coords(6, 0);
    for (int i = 0; i < voxel_dim; i += tile_dim)
        for (int j = 0; j < voxel_dim; j += tile_dim)
            for (int k = 0; k < voxel_dim; k += tile_dim) {
                // reads one tile of TileDB array 'skeleton'
                tile_coords = {i, i+tile_dim-1, j, j+tile_dim-1, k, k+tile_dim-1};
                read_skeleton(tile_coords, sk_data, sk_coords);
                
                if (sk_coords.size() == 0) continue; // skips if tile is empty
                
                // creates a temporary AF array and dilates it
                block = create_af_array(sk_coords, tile_coords);
                di_block = dilate_af_array(block, true);
                
                // saves the dilated AF array in TileDB array 'dilation'
                get_nz_coords(di_block, di_coords, tile_coords);
                di_data.resize(di_coords.size()/3, '1');
                write_dilation(di_data, di_coords);
            }
}

void Foam::save_slice_image(const char* filename, const int height, const int axis) {
    std::vector<int> slice_coords;
    if (axis == 0)
        slice_coords = {height, height, 0, voxel_dim-1, 0, voxel_dim-1};
    else if (axis == 1)
        slice_coords = {0, voxel_dim-1, height, height, 0, voxel_dim-1};
    else
        slice_coords = {0, voxel_dim-1, 0, voxel_dim-1, height, height};
    
    std::vector<char> data(1, '1'); // dummy initialization
    // TileDB bug: sometimes throws an error if uninitialized
    std::vector<int> coords(3, 0);
    read_dilation(slice_coords, data, coords);
    array block = create_af_array(coords, slice_coords);
    
    if (axis == 0)
        saveImage(filename, reorder(block, 1, 2, 0));
    else if (axis == 1)
        saveImage(filename, reorder(block, 2, 0, 1));
    else
        saveImage(filename, block);
}

void Foam::create_skeleton(const int voxel_dim, const int tile_dim) {
    create_tiledb_array(skeleton, voxel_dim, tile_dim);
}

void Foam::create_dilation(const int voxel_dim, const int tile_dim) {
    create_tiledb_array(dilation, voxel_dim, tile_dim);
}

void Foam::create_tiledb_array(const std::string& array_name, const int
            voxel_dim, const int tile_dim) {
    // Create a TileDB context.
    Context ctx;

    // The array will be 3D with "rows", "cols", and "pages".
    Domain domain(ctx);
    domain.add_dimension(Dimension::create<int>(ctx, "rows", {{0, voxel_dim-1}}, tile_dim))
        .add_dimension(Dimension::create<int>(ctx, "cols", {{0, voxel_dim-1}}, tile_dim))
        .add_dimension(Dimension::create<int>(ctx, "pages", {{0, voxel_dim-1}}, tile_dim));

    // The array will be sparse.
    ArraySchema schema(ctx, TILEDB_SPARSE);
    schema.set_domain(domain).set_order({{TILEDB_ROW_MAJOR, TILEDB_ROW_MAJOR}});

    // Add a single attribute 'attri' so each (i,j,k) cell can store a character.
    schema.add_attribute(Attribute::create<char>(ctx, attri));

    // Create the (empty) array on disk.
    Array::create(array_name, schema);
}

void Foam::write_skeleton(std::vector<char>& data, std::vector<int>& coords) {
    write_tiledb_array(skeleton, data, coords);
}

void Foam::write_dilation(std::vector<char>& data, std::vector<int>& coords) {
    write_tiledb_array(dilation, data, coords);
}

void Foam::write_tiledb_array(const std::string& array_name, std::vector<char>& data,
            std::vector<int>& coords) {
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

void Foam::read_skeleton(const std::vector<int>& tile_coords, std::vector<char>& data,
            std::vector<int>& coords) {
    read_tiledb_array(skeleton, tile_coords, data, coords);
}

void Foam::read_dilation(const std::vector<int>& tile_coords, std::vector<char>& data,
            std::vector<int>& coords) {
    read_tiledb_array(dilation, tile_coords, data, coords);
}

void Foam::read_tiledb_array(const std::string& array_name, const std::vector<int>& tile_coords,
            std::vector<char>& data, std::vector<int>& coords) {
    // Create a TileDB context.
    Context ctx;

    // Prepare the array for reading
    Array array(ctx, array_name, TILEDB_READ);

    // Prepare the vector that will hold the result.
    // We take an upper bound on the result size, as we do not
    // know a priori how big it is (since the array is sparse)
    auto max_el = array.max_buffer_elements(tile_coords);
    data.resize(max_el[attri].second);
    coords.resize(max_el[TILEDB_COORDS].second);

    // Prepare the query
    Query query(ctx, array, TILEDB_READ);
    query.set_subarray(tile_coords)
        .set_layout(TILEDB_ROW_MAJOR)
        .set_buffer(attri, data)
        .set_coordinates(coords);

    // Submit the query and close the array.
    query.submit();
    array.close();

    // Print out the results.
    /*std::cout << array_name << " data: ";
    std::cout << "(" << tile_coords[0] << ", " << tile_coords[2] << ", "
            << tile_coords[4] << ") - ";
    std::cout << "- (" << tile_coords[1] << ", " << tile_coords[3] << ", "
            << tile_coords[5] << ")" << std::endl;*/
    auto result_num = (int)query.result_buffer_elements()[attri].second;
    //std::cout << "Number of entries = " << result_num << std::endl;
    data.resize(result_num);
    coords.resize(result_num*3);
}

/*
void Foam::spatial_to_voxel(const std::vector<float>& s_coords, std::vector<int>& v_coords) {
    for (int i = 0; i < 3; i++)
        v_coords[i] = offset + (int) round((voxel_dim-2*offset)*
                                (s_coords[i] - bb[i])/(bb[i+3] - bb[i]));
}

void Foam::voxel_to_spatial(const std::vector<int>& v_coords, std::vector<float>& s_coords) {
    for (int i = 0; i < 3; i++)
        s_coords[i] = bb[i] + (bb[i+3] - bb[i])*(v_coords[i]-offset)/(voxel_dim-2*offset);
}
*/
