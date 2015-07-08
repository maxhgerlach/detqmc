// use cnpy from Github to save an Armadillo cube (3d array) to Numpy format
#include <iostream>
#include <armadillo>
#include "cnpy.h"

arma::cube transpose_3d(const arma::cube& orig) {
    arma::cube result(orig.n_slices,
                      orig.n_cols,
                      orig.n_rows); // flipped dimensions
    for (int s = 0; s < orig.n_slices; ++s) {
        for (int r = 0; r < orig.n_rows; ++r) {
            for (int c = 0; c < orig.n_cols; ++c) {
                result(s, c, r) = orig(r, c, s);
            }
        }
    }
    return result;
}





int main() {
    int nz = 3;
    int ny = 2;
    int nx = 4;
    arma::cube test(nz, ny, nx);
    
    int i = 1;
    for (int x = 0; x < nx; ++x) {
        for (int z = 0; z < nz; ++z) {
            for (int y = 0; y < ny; ++y) {
                test(z, y, x) = i;
                ++i;
            }
        }
    }
    test.print("test");

    arma::cube transp = transpose_3d(test);
    transp.print("transp");

    // save to file
    const unsigned int shape[] = {nz,ny,nx};
    cnpy::npy_save("test.npy",transp.memptr(),shape,3,"w");

    return 0;
}
