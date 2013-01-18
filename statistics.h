/*
 * statistics.h
 *
 *  Created on: Dec 13, 2012
 *      Author: gerlach
 */

#ifndef STATISTICS_H_
#define STATISTICS_H_

#include <cmath>
#include <vector>
#include <tuple>
#include <armadillo>


//There is support for Armadillo vector/matrix values in this functions. But this requires
//passing an instance that represents 0 (correctly sized object)

template<typename T>
T variance(const std::vector<T>& numbers, T meanValue, const T& zeroValue = T()) {
	using std::pow;
	using arma::pow;
	unsigned N = numbers.size();
    T sum = zeroValue;
    for (unsigned i = 0; i < N; ++i) {
        sum += pow(numbers[i] - meanValue, 2);
    }
    return sum / (N-1);
}


//Take a vector of block values, estimate their error using standard jackknife
//use this if the average is already known
template<typename T>
T jackknife(
        const std::vector<T>& blockValues, T blockAverage, const T& zeroValue = T()) {
    using std::pow;
    using std::sqrt;
    using arma::pow;
    using arma::sqrt;
	unsigned bc = blockValues.size();
    T squaredDeviation = zeroValue;
    for (unsigned b = 0; b < bc; ++b) {
        squaredDeviation += pow(blockAverage - blockValues[b], 2);
    }
    return sqrt(double(bc - 1) / double(bc) * squaredDeviation);
}


//Take a vector of block values, calculate their average and estimate
//their error using standard jackknife
//return a tuple [average, error]
template<typename T>
std::tuple<T,T> jackknife(const std::vector<T>& blockValues, const T& zeroValue = T()) {
    T outBlockAverage = zeroValue;
	unsigned bc = blockValues.size();
    for (unsigned b = 0; b < bc; ++b) {
        outBlockAverage += blockValues[b];
    }
    outBlockAverage /= static_cast<T>(bc);

    T outBlockError = jackknife(blockValues, outBlockAverage);

    return std::make_tuple(outBlockAverage, outBlockError);
}





#endif /* STATISTICS_H_ */
