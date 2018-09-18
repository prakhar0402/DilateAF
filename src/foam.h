#include <iostream>
#include <math.h>
#include <tiledb/tiledb>
#include <arrayfire.h>
#include <af/util.h>

using namespace af;
using namespace tiledb;

// This class creates two TileDB arrays - one for skeleton structure and the other for its dilation
// It uses ArrayFire to dilate skeleton in chunks efficiently and unions the output
class Foam {
private:
    // name of TileDB attribute
    std::string attri;
    
    // name of TileDB arrays
    std::string skeleton;
    std::string dilation;
    
    // number of voxels along each axis in the array
    int voxel_dim; //TODO: create cuboid voxel array with different edge lengths
    
    // number of voxels along each axis in a tile
    int tile_dim;
    
    // radius of mask used for dilation
    int radius;
    
    // mask used for dilation
    array mask;
    
    // creates a 3D spherical mask of given radius (number of voxels)
    void create_mask();
    
    // returns a 3D boolean AF array with global coordinates 'coords' set to 'true'
    array create_af_array(const std::vector<int>& coords,
            const std::vector<int>& tile_coords);
    
    // returns dilated AF array - with or without padding
    array dilate_af_array(const array& block, bool usePad = true);
    
    // gets non-zero (true) voxel coordinates from AF array
    // uses tile_coords and mask_dims to transform local to global coordinates
    void get_nz_coords(const array& block, std::vector<int>& coords,
            const std::vector<int>& tile_coords);
    
    // spatial bounding box
    //std::vector<float> bb;
    
    // creates TileDB array
    void create_tiledb_array(const std::string& array_name, const int voxel_dim,
            const int tile_dim);
    
    // write data in TileDB array
    void write_tiledb_array(const std::string& array_name, std::vector<char>& data,
            std::vector<int>& coords);
    
    // reads data from a fragmant (tile_coords) of TileDB array
    void read_tiledb_array(const std::string& array_name, const std::vector<int>& tile_coords,
            std::vector<char>& data, std::vector<int>& coords);
    
public:
    Foam();
    Foam(int v_dim, int t_dim, int rad);
    
    // deletes old TileDB arrays and creates new ones
    void reset();
    void reset_skeleton();
    void reset_dilation();
    
    // setters
    void set_dims(int v_dim, int t_dim);
    void set_mask_radius(int rad);
    
    // dilates the skeleton
    void dilate();
    
    // saves image of a sliced plane at the specified height along an axis in dilation
    void save_slice_image(const char* filename, const int height, const int axis = 2);
    
    // create, write, and read TileDB arrays 'skeleton' and 'dilation'
    void create_skeleton(const int voxel_dim, const int tile_dim);
    void create_dilation(const int voxel_dim, const int tile_dim);
    void write_skeleton(std::vector<char>& data, std::vector<int>& coords);
    void write_dilation(std::vector<char>& data, std::vector<int>& coords);
    void read_skeleton(const std::vector<int>& tile_coords, std::vector<char>& data,
            std::vector<int>& coords);
    void read_dilation(const std::vector<int>& tile_coords, std::vector<char>& data,
            std::vector<int>& coords);

    //void spatial_to_voxel(const std::vector<float>& s_coords, std::vector<int>& v_coords);
    //void voxel_to_spatial(const std::vector<int>& v_coords, std::vector<float>& s_coords);
};
