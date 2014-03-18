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

#include "boost_serialize_armadillo.h"

//matrices used in the computation of B-matrices decomposed into
//(U,d,V) = (unitary matrix, real diagonal matrix elements >= 0, unitary matrix)
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
    bool ok = arma::svd(result.U, result.d, result.V_t, mat);
    if (not ok) {
    	throw GeneralError("SVD failed");
    }

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
