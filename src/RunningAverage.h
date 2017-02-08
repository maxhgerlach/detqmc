/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0.  See the enclosed file LICENSE for a copy or if
 * that was not distributed with this file, You can obtain one at
 * http://mozilla.org/MPL/2.0/.
 *
 * Copyright 2017 Max H. Gerlach
 * 
 * */

#include <list>
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wconversion"
#pragma GCC diagnostic ignored "-Wshadow"
#include <boost/serialization/list.hpp>
#pragma GCC diagnostic pop

namespace RA {

template<typename Val>
class RunningAverage {
    int sampleSize;
    int samplesAdded;
    std::list<Val> values; //deque was problematic when serialized in a debug build
    Val runningAverage;
public:
    RunningAverage(int sampleSize_);
    void addValue(Val v);
    Val get();

    int getSamplesAdded();

private:
    friend class boost::serialization::access;

    RunningAverage() :
        sampleSize(), samplesAdded(), values(),
        runningAverage()
    {
        //private default constructor, just for serialization
    }

    template<class Archive>
    void serialize(Archive& ar, const uint32_t version) {
        (void)version;
        ar & sampleSize & samplesAdded
           & values & runningAverage;
    }
};

template<typename Val>
RunningAverage<Val>::RunningAverage(int sampleSize_) :
    sampleSize(sampleSize_), samplesAdded(0), values(),
    runningAverage(0)
{
}

template<typename Val>
void RunningAverage<Val>::addValue(Val v) {
    if (samplesAdded < sampleSize) {
        values.push_back(v);
        runningAverage += v / sampleSize;
    } else {
        runningAverage -= values.front() / sampleSize;
        values.pop_front();
        values.push_back(v);
        runningAverage += v / sampleSize;
    }
    ++samplesAdded;
}

template<typename Val>
Val RunningAverage<Val>::get() {
    return runningAverage;
}

template<typename Val>
int RunningAverage<Val>::getSamplesAdded() {
    return samplesAdded;
}

}

typedef RA::RunningAverage<double> RunningAverage;
