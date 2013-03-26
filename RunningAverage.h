//Code by Max Henner Gerlach, 2010--2012,
//used for the diploma thesis "Directional Ordering in the Classical Compass Model in Two and Three Dimensions"
//contact: maxgerlach@gmail.com

#include <queue>

namespace RA {

template<typename Val>
class RunningAverage {
    int sampleSize;
    int samplesAdded;
    std::queue<Val> values;
    Val runningAverage;
public:
    RunningAverage(int sampleSize_);
    void addValue(Val v);
    Val get();

    int getSamplesAdded();
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
        values.push(v);     //pushes to end of queue
        runningAverage += v / sampleSize;
    } else {
        runningAverage -= values.front() / sampleSize;
        values.pop();
        values.push(v);
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
