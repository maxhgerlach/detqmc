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
#include "exceptions.h"
#include "tools.h"
#include "boost/algorithm/string.hpp"   //trimming
#include "metadata.h"

// TODO: switch to smart (shared) pointers


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
    int columns;
    std::vector<std::vector<ValueType>*>* data;
    MetadataMap meta;
public:
    DataSeriesLoader();
    void readFromFile(const std::string& filename, unsigned subsample = 1,
            unsigned discardData = 0, unsigned sizeHint = 0) throw (ReadError);

    int getColumns();
    std::vector<ValueType>* getData(int column = 0);
    void deleteData();

    template<typename MetaValueType>
    void getMeta(const std::string& key, MetaValueType& value) throw(KeyUndefined);

    MetadataMap getMetadataMap();
};

template <typename ValueType>
DataSeriesLoader<ValueType>::DataSeriesLoader() {
    columns = 0;
    data = 0;
}

template <typename ValueType>
void DataSeriesLoader<ValueType>::deleteData() {
    destroyAll(*data);
    destroy(data);
}

template <typename ValueType>
int DataSeriesLoader<ValueType>::getColumns() {
    return columns;
}


template <typename ValueType>
void DataSeriesLoader<ValueType>::readFromFile(const std::string& filename,
        unsigned subsample, unsigned discardData, unsigned sizeHint)
        throw (ReadError) {
    using namespace std;
    ifstream input(filename.c_str());
    if (not input) {
        throw ReadError(filename);
    }
    data = new vector<vector<ValueType>*>();
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
                data->push_back(new vector<ValueType>);
                if (sizeHint > 0) {
                    data->at(columns)->
                            reserve((sizeHint - discardData) / subsample);
                }
                data->at(columns)->push_back(val);
                ++columns;
            }
            //handle the rest of the lines below
            break;
        }
    }
    meta = parseMetadataBlock(configLines);

    if (discardData > 0) {
        //the line that was already read in earlier has to be discarded:
        for (int c = 0; c < columns; ++c) data->at(c)->pop_back();

        for (unsigned linesRead = 1; linesRead < discardData; ++linesRead) {
            getline(input, line);
        }
    }

    //read in the remaining data
    if (subsample > 1) {
        unsigned samples = 1;
        while (getline(input, line)) {
            if (samples == subsample) {
                samples = 0;
            }
            stringstream ss(line);
            for (int c = 0; c < columns; ++c) {
                ValueType val;
                ss >> val;
                if (samples == 0) {
                    (*data)[c]->push_back(val);
                }
            }
            ++samples;
        }
    } else if (subsample == 1) {
        while (getline(input, line)) {
            stringstream ss(line);
            for (int c = 0; c < columns; ++c) {
                ValueType val;
                ss >> val;
                (*data)[c]->push_back(val);
            }
        }
    }
}

template<typename ValueType>
std::vector<ValueType>* DataSeriesLoader<ValueType>::getData(int column) {
    return (*data)[column];
}

template<typename ValueType>
template<typename MetaValueType>
void DataSeriesLoader<ValueType>::getMeta(const std::string& key, MetaValueType& value) throw(KeyUndefined) {
    if (meta.count(key) == 0) {
        throw KeyUndefined(key);
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
void DataSeriesLoader<double>::readFromFile(const std::string& filename,
        unsigned subsample, unsigned discardData, unsigned sizeHint)
        throw (ReadError) {
    using namespace std;
    ifstream input(filename.c_str());
    if (not input) {
        throw ReadError(filename);
    }
    data = new vector<vector<double>*>();
    if (sizeHint > 0) {
        data->reserve((sizeHint - discardData) / subsample);
    }
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
                data->push_back(new vector<double>);
                if (sizeHint > 0) {
                    data->at(columns)->
                            reserve((sizeHint - discardData) / subsample);
                }
                data->at(columns)->push_back(val);
                ++columns;
            }
            //handle the rest of the lines below
            break;
        }
    }
    meta = parseMetadataBlock(configLines);

    if (discardData > 0) {
        //the line that was already read in earlier has to be discarded:
        for (int c = 0; c < columns; ++c) data->at(c)->pop_back();

        for (unsigned linesRead = 1; linesRead < discardData; ++linesRead) {
            getline(input, line);
        }
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
    unsigned long length = input.tellg() - startPos + 1;     //length of data to be read in
    input.seekg(startPos);
    char* buffer = new char[length + 1];
    input.read(buffer, length);
    buffer[length] = '\0';              //add null termination
    char* tokenPointer = buffer;

//  if (meta.count("samples")) {
//      //time series file tells us how many samples are contained -> preallocate vector
//      unsigned samples = fromString<unsigned>(meta["samples"]);
//      for (int c = 0; c < columns; ++c) {
//          (*data)[c]->resize(samples);
//      }
//      unsigned index = 1;
//      while (tokenPointer < buffer + length) {
//          for (int c = 0; c < columns; ++c) {
//              double val = strtod(tokenPointer, &tokenPointer);   //scans the double at tokenPointer, then sets tokenPointer to point behind the double
//              (*((*data)[c]))[index] = val;
//          }
//      }
//  } else {

    //number of samples not known beforehand, resize on the go:
    if (subsample > 1) {
        unsigned samples = 1;
        while (tokenPointer < buffer + length) {
            if (samples == subsample) {
                samples = 0;
            }
            for (int c = 0; c < columns; ++c) {
                double val = strtod(tokenPointer, &tokenPointer);   //scans the double at tokenPointer, then sets tokenPointer to point behind the double
                if (samples == 0) {
                    (*data)[c]->push_back(val);
                }
            }
            ++samples;
        }
    } else if (subsample == 1) {
        while (tokenPointer < buffer + length) {
            for (int c = 0; c < columns; ++c) {
                double val = strtod(tokenPointer, &tokenPointer);   //scans the double at tokenPointer, then sets tokenPointer to point behind the double
                (*data)[c]->push_back(val);
            }
        }
    }
//  }
    delete[] buffer;
}


#endif /* DATASERIESLOADER_H_ */
