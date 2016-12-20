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
