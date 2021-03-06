/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0.  See the enclosed file LICENSE for a copy or if
 * that was not distributed with this file, You can obtain one at
 * http://mozilla.org/MPL/2.0/.
 *
 * Copyright 2017 Max H. Gerlach
 * 
 * */

/*
 * dataserieswritersucc.h
 *
 *  Created on: Sep 29, 2011
 *      Author: gerlach
 */

#ifndef DATASERIESWRITERSUCC_H_
#define DATASERIESWRITERSUCC_H_

/*
 * write a container (providing an iterator) containing data to an ascii file
 * supports a header of the form:
 *
 * ## blah blah
 * ## blah
 * # foo1 = bar
 * # foo2 = 124
 * # foo4 = 28.983
 * ## palaver palaver
 * value0
 * value1
 * ...
 * etc.
 *
 * This version allows to successively add data to the file.
 *
 * Slight modification, less pointer handling [2012.12.14]
 * 
 * add writing of single data entries [2014.11.18]
 *
 */

#include <string>
#include <sstream>
#include <fstream>
#include <cassert>

#include "metadata.h"

template <class Container>
class DataSeriesWriterSuccessive {
public:
    //prepare the file, the following operations will append to it
    DataSeriesWriterSuccessive(const std::string& filename, bool appendToFile = false);

    //prepare the file at @filename,
    //copy the header from @headerSource
    DataSeriesWriterSuccessive(
            const DataSeriesWriterSuccessive& headerSource,
            const std::string& filename,
            bool appendToFile = false);

    //use the following functions to prepare the header
    template <class ValueType>
    void addMeta(const std::string& key, ValueType val);
    void addMetadataMap(const MetadataMap& meta);
    void addHeaderText(const std::string& headerText);
    //write meta data and header text to the file
    void writeHeader();

    //write the whole data from the container to the file
    void writeData(const Container& dataSeries);
    void writeData(const Container& dataSeries, uint32_t floatPrecision);
    //write a single data point to the file
    void writeData(typename Container::value_type value);
    void writeData(typename Container::value_type value, uint32_t floatPrecision);
    // or just a preformatted line:
    void writeData(const std::string& line);    
private:
    std::ofstream output;
    std::string header;
};

template <class Container>
DataSeriesWriterSuccessive<Container>::
DataSeriesWriterSuccessive(const std::string& filename, bool appendToFile)
    : output(filename.c_str(), appendToFile ? std::ios::app : std::ios::out), header("")
{
	if (not output) {
            std::cerr << "Could not open file " << filename << " for writing.\n";
            std::cerr << "Error code: " << strerror(errno) << "\n";
	}
}

template <class Container>
DataSeriesWriterSuccessive<Container>::
DataSeriesWriterSuccessive(
        const DataSeriesWriterSuccessive<Container>& headerSource,
        const std::string& filename,
        bool appendToFile)
    : output(filename.c_str(), appendToFile ? std::ios::app : std::ios::out),
      header(headerSource.header)
{}

template <class Container>
template <class ValueType>
void DataSeriesWriterSuccessive<Container>::
addMeta(const std::string& key, ValueType val) {
    std::stringstream ss;
    ss << "# " << key << " = " << val << '\n';
    header += ss.str();
}

template <class Container>
void DataSeriesWriterSuccessive<Container>::
addMetadataMap(const MetadataMap& meta) {
    header += metadataToString(meta, "# ");
}


template <class Container>
void DataSeriesWriterSuccessive<Container>::
addHeaderText(const std::string& headerText) {
    std::stringstream ss(headerText);
    std::string line;
    while (std::getline(ss, line)) {
        header += "## " + line + '\n';
    }
}

template <class Container>
void DataSeriesWriterSuccessive<Container>::
writeHeader() {
    output << header;
    output.flush();
}

template <class Container>
void DataSeriesWriterSuccessive<Container>
::writeData(const Container& data) {
    for (auto iter = data.cbegin(); iter != data.cend(); ++iter) {
        output << *iter << '\n';
    }
    output.flush();
}

template <class Container>
void DataSeriesWriterSuccessive<Container>
::writeData(const Container& data, uint32_t floatPrecision) {
    output.precision(floatPrecision);
    output.setf(std::ios::scientific, std::ios::floatfield);
    for (auto iter = data.cbegin(); iter != data.cend(); ++iter) {
        output << *iter << '\n';
    }
    output.flush();
}

template <class Container>
void DataSeriesWriterSuccessive<Container>
::writeData(typename Container::value_type value) {
    output << value << '\n';
    output.flush();
}

template <class Container>
void DataSeriesWriterSuccessive<Container>
::writeData(const std::string& line) {
    output << line << '\n';
    output.flush();
}


template <class Container>
void DataSeriesWriterSuccessive<Container>
::writeData(typename Container::value_type value, uint32_t floatPrecision) {
    output.precision(floatPrecision);
    output.setf(std::ios::scientific, std::ios::floatfield);
    output << value << '\n';
    output.flush();
}


typedef DataSeriesWriterSuccessive<std::vector<double>> DoubleVectorWriterSuccessive;
typedef DataSeriesWriterSuccessive<std::vector<int>> IntVectorWriterSuccessive;




#endif /* DATASERIESWRITERSUCC_H_ */
