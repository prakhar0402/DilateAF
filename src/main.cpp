#include <stdio.h>
#include <arrayfire.h>
#include <af/util.h>
#include <cstdlib>
#include <time.h>
#include <math.h>
#include <iostream>

using namespace af;

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

    return 0;
}
