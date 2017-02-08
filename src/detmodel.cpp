/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0.  See the enclosed file LICENSE for a copy or if
 * that was not distributed with this file, You can obtain one at
 * http://mozilla.org/MPL/2.0/.
 *
 * Copyright 2017 Max H. Gerlach
 * 
 * */

/*
 * detmodel.cpp
 *
 *  Created on: Feb 21, 2013
 *      Author: gerlach
 */

#include "detmodel.h"



MatNum computePropagator(num scalar, const MatNum& matrix) {
    using namespace arma;

    VecNum eigval;
    MatNum eigvec;
    eig_sym(eigval, eigvec, matrix);

    return eigvec * diagmat(exp(-scalar * eigval)) * trans(eigvec);
}

MatCpx computePropagator(num scalar, const MatCpx& matrix) {
    using namespace arma;

    VecNum eigval;              // hermitian matrix has real eigenvalues
    MatCpx eigvec;
    eig_sym(eigval, eigvec, matrix); // for hermitian matrix

    return eigvec * diagmat(exp(-scalar * eigval)) * trans(eigvec);
}





