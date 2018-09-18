#include <iostream>
#include <math.h>
#include <time.h>
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
    for (int i = 20; i < 236; i++) {
        coords.insert(coords.end(), {i, i, 32});
    }
    std::vector<char> data(coords.size()/3, '1');
    foam.write_skeleton(data, coords);
    
    // dilate skeleton
    clock_t start = clock();
    foam.dilate();
    float elapsed = (float)(clock() - start) / (float)CLOCKS_PER_SEC;
    std::cout << "** Dilation complete in " << elapsed << " secs! **" << std::endl;
    
    // testing
    std::vector<int> tile_coords = {0, 255, 0, 255, 29, 35};
    foam.save_tile_image_Z(tile_coords, 0, "z029.png");
    foam.save_tile_image_Z(tile_coords, 3, "z032.png");
    foam.save_tile_image_Z(tile_coords, 6, "z035.png");
    
    return 0;
}
