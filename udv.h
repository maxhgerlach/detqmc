/*
 * udv.h
 *
 *  Created on: Feb 12, 2013
 *      Author: gerlach
 */

#ifndef UDV_H_
#define UDV_H_

#include <complex>
#include <armadillo>
#include "timing.h"
#include "exceptions.h"
#include "tools.h"

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wconversion"
#pragma GCC diagnostic ignored "-Wshadow"
#include "boost_serialize_armadillo.h"
#pragma GCC diagnostic pop

//matrices used in the computation of B-matrices decomposed into
//(U,d,V) = (unitary matrix, real diagonal matrix elements >= 0, unitary matrix)
//We store U, d, V_t. V_t is the conjugate-transpose of V. The Lapack-Routines
//to compute svd(M) for a matrix M return U, d, V_t with M = U*d*V = U*d*(V_t)^(dagger).
//Initialize at beginning of simulation by member function setupUdVStorage()

template <typename num, typename num_s = double>
struct UdV {
    arma::Mat<num> U;
    arma::Col<num_s> d;
    arma::Mat<num> V_t;			//conjugate-transpose of V
    //default constructor: leaves everything empty
    UdV() : U(), d(), V_t() {}
    //specify matrix size: initialize to identity
    UdV(uint32_t size) :
        U(arma::eye(size,size)), d(arma::ones(size)), V_t(arma::eye(size,size))
    { }
private:
    //for serialization with Boost
    friend class boost::serialization::access;
    template<class Archive>
    void serialize(Archive& ar, const uint32_t /* version */) {
        ar & U & d & V_t;
    }
};

template<> inline
UdV<std::complex<double>>::UdV(uint32_t size) :
    U(arma::eye(size,size), arma::zeros(size,size)),
    d(arma::ones(size)),
    V_t(arma::eye(size,size), arma::zeros(size,size))
{ }

template <typename Val>
UdV<Val> udvDecompose(const arma::Mat<Val>& mat) {
    timing.start("udvDecompose");

    typedef UdV<Val> UdV;
    UdV result;

//  ArmaMat V_transpose;
//  arma::svd(result.U, result.d, V_transpose, mat, "standard");
//  result.V = V_transpose.t();         //potentially it may be advisable to not do this generally

    //svd-call should use divide-and-conquer algorithm.
    //mat == U * diag(d) * trans(V_t)

/*
    
    bool ok = arma::svd(result.U, result.d, result.V_t, mat, "dc");
    
    if (not ok) {
    	// try the standard method instead of divide-and-conquer
    	bool ok2 = arma::svd(result.U, result.d, result.V_t, mat, "std");
    	if (not (ok2)) {
            throw GeneralError("SVD failed (dc, then std)");
    	}
    }

*/
    

    //Use std algorithm instead -- more precise -- this leads to
    //much higher stability, with actually not much longer runtimes
    bool ok = arma::svd(result.U, result.d, result.V_t, mat, "std");
    if (not ok) {
        throw GeneralError("SVD failed (std)");
    }


    
//    print_matrix_diff(mat,
//    		(result.U * arma::diagmat(result.d) * result.V_t.t()).eval(),
//    		"SVD");

//    timing.start("qr");
//    arma::qr(result.U, result.V, mat);
//    timing.stop("qr");
//    //normalize rows of V to obtain scales in d:
//    result.d.set_size(mat.n_rows);
//    for (uint32_t rown = 0; rown < mat.n_rows; ++rown) {
//        const Val norm = arma::norm(result.V.row(rown), 2);
//        result.d[rown] = norm;
//        result.V.row(rown) /= norm;
//    }

    timing.stop("udvDecompose");

    return result;
}






#endif /* UDV_H_ */
