/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0.  See the enclosed file LICENSE for a copy or if
 * that was not distributed with this file, You can obtain one at
 * http://mozilla.org/MPL/2.0/.
 *
 * Copyright 2017 Max H. Gerlach
 * 
 * */

/*
 * pytools.h
 *
 *  Created on: Mar 27, 2014
 *      Author: max
 */

#ifndef PYTOOLS_H_
#define PYTOOLS_H_

#include <armadillo>

void python_matshow(const arma::mat& mat, const std::string& title = "");
void python_matshow2(const arma::mat& mat1, const std::string& title1,
		const arma::mat& mat2, const std::string& title2);



#endif /* PYTOOLS_H_ */
