/*
 * detmodel.cpp
 *
 *  Created on: Feb 21, 2013
 *      Author: gerlach
 */

#include "detmodel.h"
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wconversion"
#pragma GCC diagnostic ignored "-Wshadow"
#include "boost/assign/std/vector.hpp"    // 'operator+=()' for vectors
#pragma GCC diagnostic pop
#include "exceptions.h"


MatNum computePropagator(num scalar, const MatNum& matrix) {
    using namespace arma;

    VecNum eigval;
    MatNum eigvec;
    eig_sym(eigval, eigvec, matrix);

    return eigvec * diagmat(exp(-scalar * eigval)) * trans(eigvec);
}





