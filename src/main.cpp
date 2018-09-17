#include <iostream>
#include <math.h>
#include <tiledb/tiledb>
#include <arrayfire.h>
#include <af/util.h>

#include <foam.h>

using namespace af;
using namespace tiledb;

int main(int argc, char** argv) {
    int device = argc > 1 ? atoi(argv[1]) : 0;

    try {
        af::info();
        af::setDevice(device);

    } catch (af::exception& e) {
        fprintf(stderr, "%s\n", e.what());
        throw;
    }
    
    std::cout << "** Dilate using ArrayFire and TileDB **" << std::endl;
    
    // creates Foam object with voxel_dim = 256x256x256, tile_dim = 64x64x64,
    // and mask radius = 5(voxels) for dilation
    Foam foam(256, 64, 5);
    
    //foam.set_mask_radius(10);
    
    //TODO: implement line drawing
    // creates an arbitrary skeleton for testing
    std::vector<int> coords;
    for (int i = 0; i < 256; i++)
        coords.insert(coords.end(), {62, i, i});
    std::vector<char> data(coords.size()/3, '1');
    foam.write_skeleton(data, coords);
    
    // dilate skeleton
    foam.dilate();
    std::cout << "** Dilation complete! **" << std::endl;
    
    // testing
    std::vector<char> di_data;
    std::vector<int> di_coords;
    std::vector<int> tile_coords = {64, 127, 64, 127, 64, 127};
    
    foam.save_tile_image_Z(tile_coords, 0, "z064.png");
    foam.save_tile_image_Z(tile_coords, 31, "z095.png");
    foam.save_tile_image_Z(tile_coords, 63, "z127.png");
    
    return 0;
}
