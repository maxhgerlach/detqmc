//Code by Max Henner Gerlach, 2010--2012,
//used for the diploma thesis "Directional Ordering in the Classical Compass Model in Two and Three Dimensions"
//contact: maxgerlach@gmail.com

#include <deque>
#include <boost/serialization/deque.hpp>

namespace RA {

template<typename Val>
class RunningAverage {
    int sampleSize;
    int samplesAdded;
    std::deque<Val> values;
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
