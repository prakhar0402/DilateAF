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
    // save images at different heights along different axes
    foam.save_slice_image("slice_Z_029.png", 29, 2);
    foam.save_slice_image("slice_Z_032.png", 32, 2);
    foam.save_slice_image("slice_Z_035.png", 35, 2);
    foam.save_slice_image("slice_X_033.png", 33, 0);
    foam.save_slice_image("slice_Y_033.png", 33, 1);
    
    return 0;
}
