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

template<typename T>
T variance(const std::vector<T>& numbers, T meanValue, int N = 0) {
	using std::pow;
    if (N == 0) {
        N = numbers.size();
    }
    T sum = 0;
    for (int i = 0; i < N; ++i) {
        sum += pow(numbers[i] - meanValue, 2);
    }
    return sum / (N-1);
}

//Take a vector of block values, calculate their average and estimate
//their error using standard jackknife
//return a tuple [average, error]
template<typename T>
std::tuple<T,T> jackknife(const std::vector<T>& blockValues) {
    T outBlockAverage = 0;
    T outBlockError = 0;
	unsigned bc = blockValues.size();
    for (unsigned b = 0; b < bc; ++b) {
        outBlockAverage += blockValues[b];
    }
    outBlockAverage /= static_cast<T>(bc);
    T squaredDeviation = 0;
    for (unsigned b = 0; b < bc; ++b) {
        squaredDeviation += std::pow(outBlockAverage - blockValues[b], 2);
    }
    outBlockError = std::sqrt(double(bc - 1) / double(bc) * squaredDeviation);
    return std::make_tuple(outBlockAverage, outBlockError);
}

//Take a vector of block values, estimate their error using standard jackknife
//use this if the average is already known
template<typename T>
T jackknife(
        const std::vector<T>& blockValues, T blockAverage) {
    unsigned bc = blockValues.size();
    T squaredDeviation = 0;
    for (unsigned b = 0; b < bc; ++b) {
        squaredDeviation += std::pow(blockAverage - blockValues[b], 2);
    }
    return std::sqrt(double(bc - 1) / double(bc) * squaredDeviation);
}



#endif /* STATISTICS_H_ */
