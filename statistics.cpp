//Code by Max Henner Gerlach, 2010--2012,
//used for the diploma thesis "Directional Ordering in the Classical Compass Model in Two and Three Dimensions"
//contact: maxgerlach@gmail.com

/*
 * statistics.cpp
 *
 *  Created on: May 11, 2011
 *      Author: gerlach
 */

// taken some parts for SDW DQMC mrpt (2015-02-07 - )

#include <cmath>
#include "statistics.h"

using namespace std;

//if end==0: compute average over whole vector
//else compute average for elements at start, start+1, ..., end-1
double average(const std::vector<double>* vec, std::size_t start, std::size_t end) {
    if (end==0) {
        end = vec->size();
    }
//  assert(end > start);
    double avg = 0;
    for (std::size_t i = start; i < end; ++i) {
        avg += double((*vec)[i]) / (double(end)-double(start));
    }
    return avg;
}

double average(const std::vector<int>* vec, std::size_t start, std::size_t end) {
    if (end==0) {
        end = vec->size();
    }
//  assert(end > start);
    double avg = 0;
    for (std::size_t i = start; i < end; ++i) {
        avg += double((*vec)[i]) / (double(end)-double(start));
    }
    return avg;
}

//if end==0: compute square average over whole vector
//else compute sq. average for elements at start, start+1, ..., end-1
double sqAverage(const std::vector<double>* vec, std::size_t start, std::size_t end) {
    if (end==0) {
        end = vec->size();
    }
//  assert(end > start);
    double sqAvg = 0;
    for (std::size_t i = start; i < end; ++i) {
        sqAvg += std::pow(vec->at(i),2);
    }
    sqAvg /= (double(end)-double(start));
    return sqAvg;
}

double variance(const std::vector<double>* numbers, double meanValue, std::size_t N) {
    if (N == 0) {
        N = numbers->size();
    }
//  assert(N > 0);
    double sum = 0;
    for (std::size_t i = 0; i < N; ++i) {
        sum += std::pow( numbers->at(i) - meanValue, 2);
    }
    return sum / double(N-1);
}

double variance(const std::vector<int>* numbers, double meanValue, std::size_t N) {
    if (N == 0) {
        N = numbers->size();
    }
//  assert(N > 0);
    double sum = 0;
    for (std::size_t i = 0; i < N; ++i) {
        sum += std::pow( double(numbers->at(i)) - meanValue, 2);
    }
    return sum / double(N-1);
}



//integrated autocorrelation time,
//employ a self-consistent cut-off
//\tau_int = 1/2 + \sum_{k=1}^{k_max} A(k)
//k_max \approx 6*\tau_int
double tauint(const std::vector<double>* data, double selfConsCutOff, AutoCorrMap* points) {
    std::size_t m = data->size();
    double mean = average(data);
    double var = 0;
    double result = 0.5;
    for (int t = 0; t < int(m) - 1; ++t) {
        double autoCorr = 0;
        for (int k = 0; k < int(m) - t; ++k){
            autoCorr += ((*data)[k] - mean) * ((*data)[k + t] - mean);
        }
        autoCorr /= double(m - t);
        if (t == 0) {
            var = autoCorr;
        } else {
            result += autoCorr / var;
            if (t > selfConsCutOff * result)
                break;
        }
        if (points) {
            points->insert(make_pair(t, autoCorr));
        }
    }
    return result;
}

//integrated autocorrelation time,
//stop accumulating once autoCorr <= 0
double tauint_stopAtZeroCrossing(const std::vector<double>* data, AutoCorrMap* points) {
    std::size_t m = data->size();
    double mean = average(data);
    double var = 0;
    double result = 0.5;
    for (int t = 0; t < int(m) - 1; ++t) {
        double autoCorr = 0;
        for (int k = 0; k < int(m) - t; ++k){
            autoCorr += ((*data)[k] - mean) * ((*data)[k + t] - mean);
        }
        autoCorr /= (int(m) - t);
        if (t == 0) {
            var = autoCorr;
        } else {
            if (autoCorr <= 0) {
                break;
            } else {
                result += autoCorr / var;
            }
        }
        if (points) {
            points->insert(make_pair(t, autoCorr));
        }
    }
    return result;
}

//faster estimation of tauint using an adaptive integration scheme (compare [Chodera2007] pg. 38)
//(also stops at the zero crossing of autoCorr)
double tauint_adaptive(const std::vector<double>* data, AutoCorrMap* points) {
    std::size_t m = data->size();
    double mean = average(data);
    double result = 0.5;

    //compute variance
    double var = 0;
    for (std::size_t k = 0; k < m; ++k){
        var += ((*data)[k] - mean) * ((*data)[k] - mean);
    }
    var /= double(m);

    //adaptive integration of autocorrelation function
    //high time resolution for small lag times, lower resolution for higher times in the
    //slowly decaying tail of the autocorrelation function
    int i = 1;
    int t_i = 1;
    while (t_i < int(m) - 1) {
        //lag time for this step of the iteration
        t_i = 1 + i * (i - 1) / 2;

        //compute autocorrelation function
        double autoCorr = 0.0;
        for (int k = 0; k < int(m) - t_i; ++k){
            autoCorr += ((*data)[k] - mean) * ((*data)[k + t_i] - mean);
        }
        autoCorr /= double(m - t_i);
        autoCorr /= var;

        if (autoCorr <= 0) {
            break;
        } else {
            //weighted addition to estimate integrated autocorrelation time
            double t_next = 1 + (i+1) * (i) / 2;
            result += autoCorr * (t_next - t_i);
        }
        if (points) {
            points->insert(make_pair(t_i, autoCorr));
        }

        i = i + 1;
    }

    return result;
}



