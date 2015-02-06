//Code by Max Henner Gerlach, 2010--2012,
//used for the diploma thesis "Directional Ordering in the Classical Compass Model in Two and Three Dimensions"
//contact: maxgerlach@gmail.com

// generalized / adapted for SDW DQMC (2015-02-06 - )

/*
 * histogram.h
 *
 *  Created on: Apr 19, 2011
 *      Author: gerlach
 */

#ifndef HISTOGRAM_H_
#define HISTOGRAM_H_

#include <algorithm>
#include <cmath>
#include <map>
#include <fstream>
#include <string>
#include <vector>
#include <dlib/string.h>
#include "myexception.h"
#include "logval.h"
#include "metadata.h"
#include "tools.h"

template<typename KeyType, typename ValType>
struct HistogramT {
    typedef std::map<KeyType, ValType> HistogramMap;

    HistogramMap histo;
    HistogramMap errors;
    MetadataMap meta;
    std::string headerLines;        //each line starts with ##

    //parameters of the simulation the histogram was measured at
    //-- optional, may not even be relevant (e.g. replica-specific
    //histograms..)
    int N;              //system volume
    double beta;        //inverse temperature


    ValType total;      //all histogram values summed up
    KeyType spacing;        //distance between two bins
    KeyType minBin;
    KeyType maxBin;
    unsigned binCount;

    //default copy constructor / assignment operator work well for this

    void load(const char* inFname) throw(ReadError) {
        using namespace std;
        ifstream input(inFname);
        if (not input) {
            throw ReadError(string(inFname));
        }
//      cout << inFname << endl;
        total = ValType(0);
        binCount = 0;
        headerLines = "";
        string line;
        string metadataLines;
        while (not input.eof()) {
            getline(input, line);
            line = dlib::ltrim(line);
            if (line[0] == '#') {
                if (line[1] == '#')  {
                    //comment line -- header text
                    headerLines += line;
                } else {
                    //lines starting with # contain no histogram data
                    //interpret these as configuration data (specifing inverse temperature etc)
                    line[0] = ' ';
                    line += '\n';
                    metadataLines += line;
                }
            }
            else if (line != "") {
                //skip empty lines all together
                //from the remaining parse [bin], [count] pairs
                istringstream ss(line);
                KeyType key;
                ValType value;
                ss >> key;
                ss >> value;
//              std::cout << "///" << key << "///" << value << std::endl;
                histo[key] = value;
    //          cout << line << "- " << value << endl;
                total += value;
                ++binCount;
            }
        }
        meta = parseMetadataBlock(metadataLines);
        if (meta.size()) {
            beta = dlib::sa = meta["beta"];
            spacing = dlib::sa = meta["spacing"];
            N = dlib::sa = meta["N"];
        }
        if (meta.count("min")) {
            minBin = dlib::sa = meta["min"];
        } else {
            minBin = histo.begin()->first;      //first key
        }
        if (meta.count("max")) {
            maxBin = dlib::sa = meta["max"];
        } else {
            maxBin = histo.rbegin()->first;     //last key
        }
    }

    void assignVector(const std::vector<ValType>& vec, KeyType minVal, KeyType maxVal, double beta_, int N_) {
        histo.clear();
        beta = beta_;
        N = N_;
        binCount = vec.size();
        minBin = minVal;
        maxBin = maxVal;
        const KeyType SMALL = static_cast<KeyType>(1e-10);
        KeyType binSize = (maxBin - minBin + SMALL) / binCount;
        spacing = binSize;
        total = 0;
        KeyType bin = minBin;
        for (unsigned m = 0; m < binCount; ++m) {
            histo[bin] = vec[m];
            total += histo[bin];
            bin += binSize;
        }
    }

    void assignErrorBarVector(const std::vector<ValType>& errorBars) {
        KeyType bin = minBin;
        for (unsigned m = 0; m < binCount; ++m) {
            errors[bin] = errorBars[m];
            bin += spacing;
        }
    }

    void updateMeta() {
        meta["beta"] = numToString(beta);
        meta["spacing"] = numToString(spacing);
        meta["min"] = numToString(minBin);
        meta["max"] = numToString(maxBin);
        meta["N"] = numToString(N);
    }

    void save(const char* outFname) {
        using namespace std;
        ofstream output(outFname);
        output << headerLines;
        output << metadataToString(meta, "#");
        if (errors.empty()) {
            typename HistogramMap::const_iterator p;
            for (p = histo.begin(); p != histo.end(); ++p) {
                output << p->first << '\t' << p->second << '\n';
            }
        } else {
            typename HistogramMap::const_iterator p, q;
            for (p = histo.begin(), q = errors.begin(); p != histo.end() and q != errors.end(); ++p, ++q) {
                output << p->first << '\t' << p->second;
                output << '\t' << q->second;
                output << '\n';
            }
        }
    }

    void save(const std::string& outFname) {
        save(outFname.c_str());
    }

    //uses as key the smallest value of the bin interval (add step/2 to get the intermediate value)
    void binData(std::vector<KeyType>* timeSeries, KeyType minBinTarget, KeyType step) {
        using namespace std;
        histo.clear();
        for (unsigned i = 0; i < timeSeries->size(); ++i) {
            KeyType timeSeriesValue = (*timeSeries)[i];
            KeyType histogramKey = ceil((timeSeriesValue - minBinTarget) / step) * step + minBinTarget;
            histo[histogramKey] += 1;
        }
        minBin = minBinTarget;
        spacing = step;
        maxBin = histo.end()->first;
        total = timeSeries->size();
    }
};

typedef HistogramT<double, LogVal> HistogramLog;        //logarithmic values to avoid high data range problems
typedef HistogramT<double, double> HistogramDouble;     //faster calculations


//return area (up to a factor of binsize) enclosed by both histograms
template<typename KeyType, typename ValType>
ValType calcOverlap(HistogramT<KeyType, ValType>& h1, HistogramT<KeyType, ValType>& h2) {
    using namespace std;
    typedef HistogramT<KeyType, ValType> HT;
    typedef typename HT::HistogramMap HM;
    typedef typename HM::const_iterator Iter;
    ValType count = 0;
    Iter p1, p2;
    p1 = h1.histo.begin();
    p2 = h2.histo.begin();
    while (p1 != h1.histo.end() and p2 != h2.histo.end()) {
        KeyType k1 = p1->first;
        KeyType k2 = p2->first;
        while (k2 < k1 and p2 != h2.histo.end()) {
            ++p2;
            k2 = p2->first;
        }
        if (pow(k1 - k2, 2) <= pow(h1.spacing / 2.0, 2)) {
            //energies considered equal
            double c1 = p1->second / h1.total;
            double c2 = p2->second / h2.total;
            count += std::min(c1, c2);
        }
        ++p1;
        ++p2;
    }
    return count;
}

//If the histogram has two peaks, return the absolute value of the difference
//of their heights. If the histogram has only one peak, return its height.
//   tolerance : factor between 0 and 1 to compensate noisy data.
//     e.g.: tolerance = 10% demand that there is a dip between the two peaks
//           that's at least 10% lower than both peaks -> else no dip is
//           recognized and histogramPeakDiff gives the max peak height
//This works only for non-negative histogram entries
inline double histogramPeakDiff(const HistogramDouble* histogram,
        double tolerance) {
    //ugly C++-isms:
    typedef HistogramDouble::HistogramMap::const_iterator Iter;
    const HistogramDouble::HistogramMap &histoMap = histogram->histo;

    //values of -1.0 in the following mean: "not found yet"
    double max1 = -1.0;      //first peak value
    double  dip = -1.0;      //a histogram value deemed lower than both peaks
    double max2 = -1.0;      //value of peak after dip

    double globalMaximum = (std::max_element(histoMap.begin(), histoMap.end(),
            mapValueCompare<double,double>))->second;

    //iterate over histogram from low to high keys:
    for (Iter p = histoMap.begin(); p != histoMap.end(); ++p) {
        double val = p->second;
        if (dip < 0.0) {
            //dip after first peak not located yet
            if (val > max1) {
                //found new candidate for first peak
                max1 = val;
            } else if ((max1 - val) > tolerance * globalMaximum) {
                //found "dip" (if there is no max2 > dip afterwards,
                //this is not an actual dip, but just some value lower than
                //the peak)
                dip = val;
            }
        } else {
            //max1 and dip have been set
            if (max2 < 0.0) {
                //no new rise after the dip found so far
                if (val < dip) {
                    //still in the dip between peaks, update its "depth":
                    dip = val;
                } else if ((val - dip)  > tolerance * globalMaximum) {
                    //found new ascent after dip
                    max2 = val;
                }
            } else {
                //found a new hike after the dip --> find its maximum
                if (val > max2) {
                    max2 = val;
                }
            }
        }
    }

    if (max2 >= 0.0) {
        //two peaks, return abs(their difference)
//        std::cout << "Debug: max1=" << max1 << " max2=" << max2 << std::endl;

        return std::fabs(max2 - max1);
    } else {
        //only one peak, return its value
        return max1;
    }
}

//If this is a two-peak histogram, return max(max1, max2) / (min).
//If there is only one peak, return 1.0
//tolerance: require a dip of at least (0<tolerance<1)*max(max1,max2) below max1
//(just like for the above function)
inline double histogramRelativeDip(const HistogramDouble* histogram,
        double tolerance) {
    //lots of code copied and slightly extended from the above function
    typedef HistogramDouble::HistogramMap::const_iterator Iter;
    const HistogramDouble::HistogramMap &histoMap = histogram->histo;
    double max1 = -1.0;
    double  dip = -1.0;
    double max2 = -1.0;
    //locations of the maxima
    Iter max1iter = histoMap.begin();
    Iter max2iter = histoMap.end();

    double globalMaximum = (std::max_element(histoMap.begin(), histoMap.end(),
            mapValueCompare<double,double>))->second;

    for (Iter p = histoMap.begin(); p != histoMap.end(); ++p) {
        double val = p->second;
        if (dip < 0.0) {
            if (val > max1) {
                max1 = val;
                max1iter = p;
            } else if ((max1 - val) > tolerance * globalMaximum) {
                dip = val;
//                std::cout << "Debug: dip=" << dip << "@" << p->first
//                        << std::endl;
            }
        } else {
            if (max2 < 0.0) {
                if (val < dip) {
                    dip = val;
                } else if ((val - dip) > tolerance * globalMaximum) {
                    max2 = val;
                    max2iter = p;
                }
            } else {
                if (val > max2) {
                    max2 = val;
                    max2iter = p;
                }
            }
        }
    }
//    std::cout << "Debug: max1=" << max1 << "@" << max1iter->first
//            << " max2=" << max2 << "@" << max2iter->first << std::endl;
    if (max2 >= 0.0) {
        //found two maxima, search minium in between
        double min = max1;
//        double debugMinLoc = max1iter->first;
        for (Iter p = max1iter; p != max2iter; ++p) {
            double val = p->second;
            if (val < min) {
                min = val;
//                debugMinLoc = p->first;
            }
        }
//        std::cout << "Debug: min=" << min << "@" << debugMinLoc << std::endl;
        return std::max(max1, max2) / min;
    } else {
        //only one peak
        return 1.0;
    }
}


//If this is a two-peak histogram, return the absolute value of the
//difference of the weights between of both peaks.
//Peaks are assumed to be separated at key cutOff.
//If this is a single-peak histogram, return the whole weight (~ 1.0, if
//normalized)
//Calculate weight as sum over histogram values --> approximation
//of integral if bin width is unity
inline double histogramWeightDiff(const HistogramDouble* histogram,
        double cutOff) {
    //lots of code copied and slightly extended from the above function
    typedef HistogramDouble::HistogramMap::const_iterator Iter;
    const HistogramDouble::HistogramMap &histoMap = histogram->histo;

    Iter p = histoMap.begin();
    double key = p->first;
    double val = p->second;
    double weight1 = 0.0;
    while (p != histoMap.end() and key < cutOff) {
        weight1 += val;
        ++p;
        key = p->first;
        val = p->second;
    }
    double weight2 = 0.0;
    while (p != histoMap.end()) {
        weight2 += p->second;
        ++p;
//        key = p->first;
//        val = p->second;
    }

    return std::fabs(weight2 - weight1);
}

//If this is a two-peak histogram return the key of the minimum value between
//the two peaks (tolerance is used like in the above functions)
//[ugly code doubling]
inline double histogramMinimumLocation(const HistogramDouble* histogram,
        double tolerance) {
    //lots of code copied and slightly extended from the above function
    typedef HistogramDouble::HistogramMap::const_iterator Iter;
    const HistogramDouble::HistogramMap &histoMap = histogram->histo;
    double max1 = -1.0;
    double  dip = -1.0;
    double max2 = -1.0;
    //locations of the maxima
    Iter max1iter = histoMap.begin();
    Iter max2iter = histoMap.end();

    double globalMaximum = (std::max_element(histoMap.begin(), histoMap.end(),
            mapValueCompare<double,double>))->second;

    for (Iter p = histoMap.begin(); p != histoMap.end(); ++p) {
        double val = p->second;
        if (dip < 0.0) {
            if (val > max1) {
                max1 = val;
                max1iter = p;
            } else if ((max1 - val) > tolerance * globalMaximum) {
                dip = val;
            }
        } else {
            if (max2 < 0.0) {
                if (val < dip) {
                    dip = val;
                } else if ((val - dip) > tolerance * globalMaximum) {
                    max2 = val;
                    max2iter = p;
                }
            } else {
                if (val > max2) {
                    max2 = val;
                    max2iter = p;
                }
            }
        }
    }
    double minLoc = histoMap.begin()->first; //already the result for 1-peak
    if (max2 >= 0.0) {
        //found two maxima, search minium in between
        double min = max1;
        minLoc = max1iter->first;
        for (Iter p = max1iter; p != max2iter; ++p) {
            double val = p->second;
            if (val < min) {
                min = val;
                minLoc = p->first;
            }
        }
    }
    return minLoc;
}


//to sort by beta
template<typename KeyType, typename ValType>
bool histogramCompare(const HistogramT<KeyType, ValType>& h1, const HistogramT<KeyType, ValType>& h2) {
    return (h1.beta < h2.beta);
}

template<typename KeyType, typename ValType>
bool pHistogramCompare(const HistogramT<KeyType, ValType>* h1, const HistogramT<KeyType, ValType>* h2) {
    return (h1->beta < h2->beta);
}


template bool histogramCompare<double,double>(const HistogramDouble&, const HistogramDouble&);



#endif /* HISTOGRAM_H_ */
