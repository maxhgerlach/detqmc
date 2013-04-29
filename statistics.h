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


// Compute jackknife-block-wise estimates of the average of data
// TODO: extended version that applies a functor to each data-point
template<typename T>
std::vector<T> jackknifeBlockEstimates(const std::vector<T>& data, unsigned jkBlocks) {
	std::vector<T> blockEstimates(jkBlocks, T(0));
	unsigned jkBlockSize = data.size() / jkBlocks;

	for (unsigned i = 0; i < data.size(); ++i) {
		T value = data[i];
		unsigned curBlock = i / jkBlockSize;
		for (unsigned jb = 0; jb < jkBlocks; ++jb) {
			if (jb != curBlock) {
				blockEstimates[jb] += value;
			}
		}
	}

	unsigned jkTotalSamples = data.size() - jkBlockSize;
	for (T& blockEstimate : blockEstimates) {
		blockEstimate /= T(jkTotalSamples);
	}

	return blockEstimates;
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


//if end==0: compute average over whole vector
//else compute average for elements at start, start+1, ..., end-1
template<typename T>
T average(const std::vector<T>& vec, std::size_t start = 0, std::size_t end = 0) {
    if (end==0) {
        end = vec.size();
    }
//  assert(end > start);
    T avg = T(0);
    for (std::size_t i = start; i < end; ++i) {
        avg += T(vec[i]);
    }
    avg /= (end-start);
    return avg;
}



//integrated autocorrelation time,
//employ a self-consistent cut-off
//\tau_int = 1/2 + \sum_{k=1}^{k_max} A(k)
//k_max \approx 6*\tau_int
template<typename T>
T tauint(const std::vector<T>& data, T selfConsCutOff = T(6)) {
    std::size_t m = data.size();
    T mean = average(data);
    T var = 0;
    T result = 0.5;
    for (std::size_t t = 0; t < m - 1; ++t) {
        T autoCorr = 0;
        for (size_t k = 0; k < m - t; ++k){
            autoCorr += (data[k] - mean) * (data[k + t] - mean);
        }
        autoCorr /= (m - t);
        if (t == 0) {
            var = autoCorr;
        } else {
            result += autoCorr / var;
            if (t > selfConsCutOff * result)
                break;
        }
    }
    return result;
}

//integrated autocorrelation time,
//stop accumulating once autoCorr <= 0
template<typename T>
T tauint_stopAtZeroCrossing(const std::vector<T>& data) {
    std::size_t m = data.size();
    T mean = average(data);
    T var = 0;
    T result = 0.5;
    for (int t = 0; t < m - 1; ++t) {
        T autoCorr = 0;
        for (int k = 0; k < m - t; ++k){
            autoCorr += (data[k] - mean) * (data[k + t] - mean);
        }
        autoCorr /= (m - t);
        if (t == 0) {
            var = autoCorr;
        } else {
            if (autoCorr <= 0) {
                break;
            } else {
                result += autoCorr / var;
            }
        }
    }
    return result;
}

//faster estimation of tauint using an adaptive integration scheme (compare [Chodera2007] pg. 38)
//(also stops at the zero crossing of autoCorr)
template<typename T>
T tauint_adaptive(const std::vector<T>& data) {
    std::size_t m = data.size();
    T mean = average(data);
    T result = 0.5;

    //compute variance
    T var = 0;
    for (size_t k = 0; k < m; ++k){
        var += (data[k] - mean) * (data[k] - mean);
    }
    var /= m;

    //adaptive integration of autocorrelation function
    //high time resolution for small lag times, lower resolution for higher times in the
    //slowly decaying tail of the autocorrelation function
    size_t i = 1;
    size_t t_i = 1;
    while (t_i < m - 1) {
        //lag time for this step of the iteration
        t_i = 1 + i * (i - 1) / 2;

        //compute autocorrelation function
        T autoCorr = 0.0;
        for (size_t k = 0; k < m - t_i; ++k){
            autoCorr += (data[k] - mean) * (data[k + t_i] - mean);
        }
        autoCorr /= (m - t_i);
        autoCorr /= var;

        if (autoCorr <= 0) {
            break;
        } else {
            //weighted addition to estimate integrated autocorrelation time
            T t_next = 1 + (i+1) * (i) / 2;
            result += autoCorr * (t_next - t_i);
        }
        i = i + 1;
    }

    return result;
}







#endif /* STATISTICS_H_ */
