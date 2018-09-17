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

class Foam
{
private:
    // Name of TileDB attribute
    std::string attri;
    
    int voxel_dim;
    int tile_dim;
    int offset;
    
    // spatial bounding box
    std::vector<float> bb;
    
public:
    // Name of TileDB arrays
    std::string skeleton;
    std::string dilated;
    
    Foam();
    
    void dilate(const array& mask);// dilates 3D array block using the 3D mask
    array dilate_block(const array& block, const array& mask, bool usePad);
    array create_af_array(const int dim, const std::vector<int>& coords);
    void get_true_coords(const array& block, std::vector<int>& coords,
                            std::vector<int>& tile_coords, dim4 mask_dims);
    
    void spatial_to_voxel(const std::vector<float>& s_coords, std::vector<int>& v_coords);
    void voxel_to_spatial(const std::vector<int>& v_coords, std::vector<float>& s_coords);
    
    void create_tiledb_array(const std::string array_name, const int voxel_dim,
                                const int tile_dim);
    void write_tiledb_array(const std::string array_name, std::vector<char>& data,
                                std::vector<int>& coords);
    void read_tiledb_array(const std::string array_name, const std::vector<int> slice_coords,
                                std::vector<char>& data, std::vector<int>& coords);
    void read_tiledb_array_relative(const std::string array_name, const std::vector<int> slice_coords,
                                std::vector<char>& data, std::vector<int>& coords);
};
