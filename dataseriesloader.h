//Code by Max Henner Gerlach, 2010--2012,
//used for the diploma thesis "Directional Ordering in the Classical Compass Model in Two and Three Dimensions"
//contact: maxgerlach@gmail.com

/*
 * dataseriesloader.h
 *
 *  Created on: Jan 28, 2011
 *      Author: gerlach
 */

#ifndef DATASERIESLOADER_H_
#define DATASERIESLOADER_H_

#include <string>
#include <vector>
#include <fstream>
#include <sstream>
#include <cstdlib>
#include <memory>
#include "exceptions.h"
#include "tools.h"
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wconversion"
#pragma GCC diagnostic ignored "-Wshadow"
#include "boost/algorithm/string.hpp"   //trimming
#pragma GCC diagnostic pop
#include "metadata.h"


//reads in ASCII formated data
//supports multiple columns in the format
//data0[i]  data1[i]    data2[i]    {...}   dataM[i]
//finds number of columns from the first line of data.
//also supports a block of meta data at the beginning in the format
// # key = val
//if subsample > 1: take only one from @subsample lines
//if discardData > 0: leave out the first @discardData lines
//if sizeHint > 0: preallocate space for (sizeHind - discardData) / subsample
//   samples
template <typename ValueType>
class DataSeriesLoader {
public:
    DataSeriesLoader();
    void readFromFile(const std::string& filename, uint32_t subsample = 1,
                      uint32_t discardData = 0, // data points to discard at beginning
                      uint32_t readMaxData = 0, // If > 0: maximum number of data points to read,
                                                // starting after the discarded part.
                                                // If = 0: don't stop reading ever
                      uint32_t sizeHint = 0);

    int getColumns();           // if this is 0: time series is empty

    typedef std::shared_ptr<std::vector<ValueType>> TimeSeriesPtr;
    TimeSeriesPtr getData(int column = 0);
    void deleteData();

    template<typename MetaValueType>
    void getMeta(const std::string& key, MetaValueType& value);

    MetadataMap getMetadataMap();
private:
    int columns;
    std::vector<TimeSeriesPtr> data;
    MetadataMap meta;
};

template <typename ValueType>
DataSeriesLoader<ValueType>::DataSeriesLoader() : columns(0) {
}

template <typename ValueType>
void DataSeriesLoader<ValueType>::deleteData() {
    // should not require to do anything
    data.clear();
    // destroyAll(*data);
    // destroy(data);
}

template <typename ValueType>
int DataSeriesLoader<ValueType>::getColumns() {
    return columns;
}


template <typename ValueType>
void DataSeriesLoader<ValueType>::readFromFile(
    const std::string& filename, uint32_t subsample, uint32_t discardData,
    uint32_t readMaxData, uint32_t sizeHint)
{
    using namespace std;
    ifstream input(filename.c_str());
    if (not input) {
        throw_ReadError(filename);
    }
    data.clear();
    string line;
    string configLines;
    while (getline(input, line)) {
        boost::algorithm::trim_left(line);
        if (line[0] == '#') {
            //interpret lines starting with # as meta data (only at the beginning of the file)
            line[0] = ' ';
            line += '\n';
            configLines += line;
        } else if (line != "") {
            //skip empty lines all together
            //from the remaining parse data

            //scan this first line of data to determine number of columns
            //also reads in the first data point
            columns = 0;
            stringstream ss(line);
            while (not ss.eof()) {
                ValueType val;
                ss >> val;
                data.push_back(TimeSeriesPtr(new vector<ValueType>));
                if (sizeHint > 0) {
                    data.at(columns)->
                        reserve((sizeHint - discardData) / subsample);
                }
                data.at(columns)->push_back(val);
                ++columns;
            }
            //handle the rest of the lines below
            break;
        }
    }
    meta = parseMetadataBlock(configLines);

    if (discardData > 0) {
        //the line that was already read in earlier has to be discarded:
        for (int c = 0; c < columns; ++c) data.at(c)->pop_back();

        for (uint32_t linesRead = 1; input and linesRead < discardData; ++linesRead) {
            getline(input, line);
        }
    }

    if (not input) {            // discarded all there is
        return;
    }

    //read in the remaining data
    uint32_t linesRead = 0;
    if (subsample > 1) {
        uint32_t samples = 1;
        while (getline(input, line) and (readMaxData == 0 or ((linesRead++) < readMaxData))) {
            if (samples == subsample) {
                samples = 0;
            }
            stringstream ss(line);
            for (int c = 0; c < columns; ++c) {
                ValueType val;
                ss >> val;
                if (samples == 0) {
                    data[c]->push_back(val);
                }
            }
            ++samples;
        }
    } else if (subsample == 1) {
        while (getline(input, line) and (readMaxData == 0 or ((linesRead++) < readMaxData))) {
            stringstream ss(line);
            for (int c = 0; c < columns; ++c) {
                ValueType val;
                ss >> val;
                data[c]->push_back(val);
            }
        }
    }
}

template<typename ValueType>
typename DataSeriesLoader<ValueType>::TimeSeriesPtr DataSeriesLoader<ValueType>::getData(int column) {
    return data[column];
}

template<typename ValueType>
template<typename MetaValueType>
void DataSeriesLoader<ValueType>::getMeta(const std::string& key, MetaValueType& value) {
    if (meta.count(key) == 0) {
        throw_KeyUndefined(key);
    } else {
        std::stringstream parseStream(meta[key]);
        parseStream >> value;
    }
}

template<typename ValueType>
MetadataMap DataSeriesLoader<ValueType>::getMetadataMap() {
    return meta;
}

typedef DataSeriesLoader<double> DoubleSeriesLoader;

//for doubles use a faster implementation
template<> inline
void DataSeriesLoader<double>::readFromFile(
       const std::string& filename, uint32_t subsample, uint32_t discardData,
       uint32_t readMaxData, uint32_t sizeHint) {
    using namespace std;
    ifstream input(filename.c_str());
    if (not input) {
        throw_ReadError(filename);
    }
    data.clear();
    //// this was nonsense
    // if (sizeHint > 0) {
    //     data->reserve((sizeHint - discardData) / subsample);
    // }
    string line;
    string configLines;
    while (getline(input, line)) {
        boost::algorithm::trim_left(line);
        if (line[0] == '#') {
            //interpret lines starting with # as meta data
            //(only at the beginning of the file)
            line[0] = ' ';
            line += '\n';
            configLines += line;
        } else if (line != "") {
            //skip empty lines all together
            //from the remaining parse data

            //scan this first line of data to determine number of columns
            columns = 0;
            stringstream ss(line);
            while (not ss.eof()) {
                double val;
                ss >> val;
                data.push_back(std::shared_ptr<vector<double>>(new vector<double>));
                if (sizeHint > 0) {
                    data.at(columns)->
                        reserve((sizeHint - discardData) / subsample);
                }
                data.at(columns)->push_back(val);
                ++columns;
            }
            //handle the rest of the lines below
            break;
        }
    }
    meta = parseMetadataBlock(configLines);

    if (discardData > 0) {
        //the line that was already read in earlier has to be discarded:
        for (int c = 0; c < columns; ++c) data.at(c)->pop_back();

        for (uint32_t linesRead = 1; input and linesRead < discardData; ++linesRead) {
            getline(input, line);
        }
    }

    if (not input) {            // discarded all there is
        return;
    }

    //read in the remaining data
    //NEW: faster method using strtod instead of more general, robust, slow stringstream
    //NEW2: read the whole remaining file into memory first, then parse that
    streampos startPos = input.tellg();
    input.seekg(-1, ios::end);                  //points onto the last character
    //seek last character that is not white space:
    while (isspace(input.peek())) {
        input.seekg(-1, ios::cur);
    }
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wconversion"
    //max length of data to be read in (in char's !)
    std::size_t length = input.tellg() - startPos + 1;
#pragma GCC diagnostic pop
    input.seekg(startPos);
    char* buffer = new char[length + 1];
    input.read(buffer, length);
    buffer[length] = '\0';              //add null termination
    char* tokenPointer = buffer;
    char* newTokenPointer = tokenPointer;

//  if (meta.count("samples")) {
//      //time series file tells us how many samples are contained -> preallocate vector
//      uint32_t samples = fromString<uint32_t>(meta["samples"]);
//      for (int c = 0; c < columns; ++c) {
//          (*data)[c]->resize(samples);
//      }
//      uint32_t index = 1;
//      while (tokenPointer < buffer + length) {
//          for (int c = 0; c < columns; ++c) {
//              double val = strtod(tokenPointer, &tokenPointer);   //scans the double at tokenPointer, then sets tokenPointer to point behind the double
//              (*((*data)[c]))[index] = val;
//          }
//      }
//  } else {

    // the following is zero if we have discarded initial entries, one otherwise
    uint32_t valuesRead = (uint32_t)data.at(0)->size(); 
    //number of samples not known beforehand, resize on the go:
    if (subsample > 1) {
        uint32_t samples = 1;
        while (tokenPointer < buffer + length and (readMaxData == 0 or (valuesRead < readMaxData))) {
            if (samples == subsample) {
                samples = 0;
            }
            for (int c = 0; c < columns; ++c) {
                //scans the double at tokenPointer, then sets newTokenPointer to point behind the double
                double val = strtod(tokenPointer, &newTokenPointer);
                if (tokenPointer == newTokenPointer) {
                    // this means no conversion took place
                    throw_GeneralError( "Could not convert token after value number" + numToString(valuesRead) );
                }
                tokenPointer = newTokenPointer; // next token will be after the current one
                if (samples == 0) {
                    data[c]->push_back(val);
                }
            }
            ++valuesRead;
            ++samples;
        }
    } else if (subsample == 1) {
        while (tokenPointer < buffer + length and (readMaxData == 0 or (valuesRead < readMaxData))) {
            for (int c = 0; c < columns; ++c) {
                //scans the double at tokenPointer, then sets newTokenPointer to point behind the double
                double val = strtod(tokenPointer, &newTokenPointer);
                if (tokenPointer == newTokenPointer) {
                    // this means no conversion took place
                    throw_GeneralError( "Could not convert token after value number" + numToString(valuesRead) );
                }
                tokenPointer = newTokenPointer; // next token will be after the current one
                data[c]->push_back(val);
            }
            ++valuesRead;
        }
    }
//  }
    delete[] buffer;
}


#endif /* DATASERIESLOADER_H_ */
