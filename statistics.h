/*
 * statistics.h
 *
 *  Created on: Dec 13, 2012
 *      Author: gerlach
 */

#ifndef STATISTICS_H_
#define STATISTICS_H_

#include <cmath>
#include <map>
#include <vector>
#include <tuple>
#include <functional>
#include <armadillo>




//There is support for Armadillo vector/matrix values in this functions. But this requires
//passing an instance that represents 0 (correctly sized object)

template<typename T>
T variance(const std::vector<T>& numbers, T meanValue, const T& zeroValue = T()) {
    typedef typename std::vector<T>::size_type size_type;
    using std::pow;
    using arma::pow;
    size_type N = numbers.size();
    T sum = zeroValue;
    for (size_type i = 0; i < N; ++i) {
        sum += pow(numbers[i] - meanValue, 2);
    }
    return sum / static_cast<T>(N-1);
}


// Compute jackknife-block-wise estimates of the values of a function applied element-wise
// to the data
template<typename T>
std::vector<T> jackknifeBlockEstimates(const std::function<T(T)>& func,
        const std::vector<T>& data, uint32_t jkBlocks) {
    std::vector<T> blockEstimates(jkBlocks, T(0));
    uint32_t jkBlockSize = static_cast<uint32_t>(data.size()) / jkBlocks;
    //if jkBlocks is not a divisor of data.size() --> some data at the end will be discarded
    uint32_t totalSamples = jkBlocks * jkBlockSize;

    for (uint32_t i = 0; i < totalSamples; ++i) {
        T value = func(data[i]);
        uint32_t curBlock = i / jkBlockSize;
        for (uint32_t jb = 0; jb < jkBlocks; ++jb) {
            if (jb != curBlock) {
                blockEstimates[jb] += value;
            }
        }
    }

    uint32_t jkTotalSamples = totalSamples - jkBlockSize;
    for (T& blockEstimate : blockEstimates) {
        blockEstimate /= T(jkTotalSamples);
    }

    return blockEstimates;
}



// Compute jackknife-block-wise estimates of the average of data
template<typename T>
std::vector<T> jackknifeBlockEstimates(const std::vector<T>& data, uint32_t jkBlocks) {
    return jackknifeBlockEstimates<T>([](T v) { return v; },    //identity lambda function
            data, jkBlocks);
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
    uint32_t bc = static_cast<uint32_t>(blockValues.size());
    T squaredDeviation = zeroValue;
    for (uint32_t b = 0; b < bc; ++b) {
        squaredDeviation += pow(blockAverage - blockValues[b], 2);
    }
    return sqrt((double(bc - 1) / double(bc)) * squaredDeviation);
}


//Take a vector of block values, calculate their average and estimate
//their error using standard jackknife
//return a tuple [average, error]
//template<typename T>
//std::tuple<T,T> jackknife(const std::vector<T>& blockValues, const T& zeroValue = T()) {
//    T outBlockAverage = zeroValue;
//  uint32_t bc = blockValues.size();
//    for (uint32_t b = 0; b < bc; ++b) {
//        outBlockAverage += blockValues[b];
//    }
//    outBlockAverage /= static_cast<T>(bc);
//
//    T outBlockError = jackknife(blockValues, outBlockAverage);
//
//    return std::make_tuple(outBlockAverage, outBlockError);
//}

//if end==0: compute average over whole vector
//else compute average for elements at start, start+1, ..., end-1
// - more generic version that applies a function to each element before taking the average
template<typename T>
T average(const std::function<T(T)>& func,
        const std::vector<T>& vec, std::size_t start = 0, std::size_t end = 0) {
    if (end==0) {
        end = vec.size();
    }
//  assert(end > start);
    T avg = T(0);
    for (std::size_t i = start; i < end; ++i) {
        avg += func(T(vec[i]));
    }
    avg /= static_cast<T>(end-start);
    return avg;
}


//if end==0: compute average over whole vector
//else compute average for elements at start, start+1, ..., end-1
template<typename T>
T average(const std::vector<T>& vec, std::size_t start = 0, std::size_t end = 0) {
    return average<T>( [](T v) { return v; }, vec, start, end );
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
        autoCorr /= static_cast<T>(m - t);
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
    var /= static_cast<T>(m);

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
        autoCorr /= static_cast<T>(m - t_i);
        autoCorr /= var;

        if (autoCorr <= 0) {
            break;
        } else {
            //weighted addition to estimate integrated autocorrelation time
            T t_next = 1 + static_cast<T>(i+1) * static_cast<T>(i) / 2;
            result += autoCorr * (static_cast<T>(t_next) - static_cast<T>(t_i));
        }
        i = i + 1;
    }

    return result;
}



// some extras for mrpt:

//if end==0: compute average over whole vector
//else compute average for elements at start, start+1, ..., end-1
double average(const std::vector<double>* vec, int start=0, int end=0);
double average(const std::vector<int>* vec, int start=0, int end=0);


//if end==0: compute square average over whole vector
//else compute sq. average for elements at start, start+1, ..., end-1
double sqAverage(const std::vector<double>* vec, int start=0, int end=0);

double variance(const std::vector<double>* numbers, double meanValue, int N=0);
double variance(const std::vector<int>* numbers, double meanValue, int N=0);


typedef std::map<int, double> AutoCorrMap;

//integrated autocorrelation time,
//employ a self-consistent cut-off
//\tau_int = 1/2 + \sum_{k=1}^{k_max} A(k)
//k_max \approx 6*\tau_int
//
//if points != 0 store the evaluated points of the auto-correlation function
//into that map
double tauint(const std::vector<double>* data, double selfConsCutOff = 6, AutoCorrMap* points = 0);

//integrated autocorrelation time,
//stop accumulating once autoCorr <= 0
//
//if points != 0 store the evaluated points of the auto-correlation function
//into that map
double tauint_stopAtZeroCrossing(const std::vector<double>* data, AutoCorrMap* points = 0);

//faster estimation of tauint using an adaptive integration scheme (compare [Chodera2007] pg. 38)
//(also stops at the zero crossing of autoCorr)
//
//if points != 0 store the evaluated points of the auto-correlation function
//into that map
double tauint_adaptive(const std::vector<double>* data, AutoCorrMap* points = 0);


//subsample inTimeSeries into outTimeSeries (appending to its end), taking only samples separated by sampleSize
template <typename T>
void subsample(const std::vector<T>& inTimeSeries, int sampleSize, std::vector<T>& outTimeSeries) {
    for (unsigned index = 0; index < inTimeSeries.size(); index += sampleSize) {
        outTimeSeries.push_back(inTimeSeries[index]);
    }
}

template <typename T1, typename T2>
void subsampleTypeCasting(const std::vector<T1>& inTimeSeries, int sampleSize, std::vector<T2>& outTimeSeries) {
    for (unsigned index = 0; index < inTimeSeries.size(); index += sampleSize) {
        outTimeSeries.push_back(static_cast<T2>(inTimeSeries[index]));
    }
}


template<typename T>
void findMinMaxMean(const std::vector<T>& timeSeries, T& min, T& max, T& mean) {
    T runningMin = timeSeries[0];
    T runningMax = timeSeries[0];
    T cummulativeSum = timeSeries[0];
    for (unsigned k = 1; k < timeSeries.size(); ++k) {
        if (timeSeries[k] < runningMin)
            runningMin = timeSeries[k];
        else if (timeSeries[k] > runningMax)
            runningMax = timeSeries[k];
        cummulativeSum += timeSeries[k];
    }
    min = runningMin;
    max = runningMax;
    mean = cummulativeSum / static_cast<T>(timeSeries.size());
}

template <typename T>
void findMinMaxMean(const std::vector<std::vector<T>*>& timeSeriesVector, T& min, T& max, T& mean) {
    T runningMin, runningMax, cummulatedMeans;
    int firstTimeSeries = timeSeriesVector.size();
    int numTimeSeries = 0;
    for (int k = 0; k < timeSeriesVector.size(); ++k) {
        if (timeSeriesVector[k] and timeSeriesVector[k]->size() > 0) {
            findMinMaxMean(*timeSeriesVector[k],
                           runningMin, runningMax, cummulatedMeans);
            firstTimeSeries = k;
            ++numTimeSeries;
            break;
        }
    }
    for (unsigned k = firstTimeSeries + 1; k < timeSeriesVector.size(); ++k) {
        T curMin, curMax, curMean;
        if (timeSeriesVector[k] and timeSeriesVector[k]->size() > 0) {
            ++numTimeSeries;
            findMinMaxMean(*timeSeriesVector[k], curMin, curMax, curMean);
            if (curMin < runningMin)
                runningMin = curMin;
            else if (curMax > runningMax)
                runningMax = curMax;
            cummulatedMeans += curMean;
        }
    }
    min = runningMin;
    max = runningMax;
    mean = cummulatedMeans / static_cast<T>(numTimeSeries);
}

template <typename T>
void findMinMax(const std::vector<T>& timeSeries, T& min, T& max) {
    T runningMin = timeSeries[0];
    T runningMax = timeSeries[0];
    for (unsigned k = 1; k < timeSeries.size(); ++k) {
        if (timeSeries[k] < runningMin)
            runningMin = timeSeries[k];
        else if (timeSeries[k] > runningMax)
            runningMax = timeSeries[k];
    }
    min = runningMin;
    max = runningMax;
}

template <typename T>
void findMinMax(const std::vector<std::vector<T>*>& timeSeriesVector, T& min, T& max) {
    T runningMin, runningMax;
    int firstTimeSeries = timeSeriesVector.size();
    for (unsigned k = 0; k < timeSeriesVector.size(); ++k) {
        if (timeSeriesVector[k] and timeSeriesVector[k]->size() > 0) {
            findMinMax(*timeSeriesVector[k], runningMin, runningMax);
            firstTimeSeries = k;
            break;
        }
    }
    for (unsigned k = firstTimeSeries + 1; k < timeSeriesVector.size(); ++k) {
        T curMin, curMax;
        if (timeSeriesVector[k] and timeSeriesVector[k]->size() > 0) {
            findMinMax(*timeSeriesVector[k], curMin, curMax);
            if (curMin < runningMin)
                runningMin = curMin;
            if (curMax > runningMax)
                runningMax = curMax;
        }
    }
    min = runningMin;
    max = runningMax;
}


#endif /* STATISTICS_H_ */
