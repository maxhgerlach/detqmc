/*
 * udv.h
 *
 *  Created on: Feb 12, 2013
 *      Author: gerlach
 */

#ifndef UDV_H_
#define UDV_H_

#include <armadillo>

//matrices used in the computation of B-matrices decomposed into
//(U,d,V) = (orthogonal matrix, diagonal matrix elements, row-normalized triangular matrix.
//Initialize at beginning of simulation by member function setupUdVStorage()
template <typename ArmaMat, typename ArmaVec>
struct UdV {
	ArmaMat U;
	ArmaVec d;
	ArmaMat V;
	//default constructor: leaves everything empty
	UdV() : U(), d(), V() {}
	//specify matrix size: initialize to identity
	UdV(unsigned size) :
		U(arma::eye(size,size)), d(arma::ones(size)), V(arma::eye(size,size))
	{ }
};

template <typename ArmaMat, typename ArmaVec>
UdV<ArmaMat, ArmaVec> udvDecompose(const ArmaMat& mat) {
	typedef UdV<ArmaMat, ArmaVec> UdV;
	UdV result;

//	ArmaMat V_transpose;
//	arma::svd(result.U, result.d, V_transpose, mat, "standard");
//	result.V = V_transpose.t();			//potentially it may be advisable to not do this generally

	arma::qr(result.U, result.V, mat);
	//normalize rows of V to obtain scales in d:
	result.d = ArmaVec(mat.n_rows);
	for (unsigned rown = 0; rown < mat.n_rows; ++rown) {
		const num norm = arma::norm(result.V.row(rown), 2);
		result.d[rown] = norm;
		result.V.row(rown) /= norm;
	}

	return result;
}




#endif /* UDV_H_ */
