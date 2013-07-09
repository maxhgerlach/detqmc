//Code by Max Henner Gerlach, 2010--2012,
//used for the diploma thesis "Directional Ordering in the Classical Compass Model in Two and Three Dimensions"
//contact: maxgerlach@gmail.com

/*
 * datamapwriter.h
 *
 *  Created on: Jan 28, 2011
 *      Author: gerlach
 */

#ifndef DATAMAPWRITER_H_
#define DATAMAPWRITER_H_

#include <memory>
#include <string>
#include <sstream>
#include <fstream>
#include <map>
#include "metadata.h"


/*
 * write a map containing key-value data to an ascii file
 * supports a header of the form:
 *
 * ## blah blah
 * ## blah
 * # foo1 = bar
 * # foo2 = 124
 * # foo4 = 28.983
 * ## palaver palaver
 * key0     data0
 * key1     data1
 * ...
 * etc.
 *
 * if errors are set the data is written out as follows instead:
 * key0     data0   error0
 * key1     data1   error1
 * ...
 *
 *
 *
 * 20121214 -- move to smart pointers (shared_ptr) to improve resource managment
 *
 */
template <typename Key, typename Value>
class DataMapWriter {
public:
    DataMapWriter();
    typedef std::shared_ptr<const std::map<Key,Value>> MapPtr;
    void setData(MapPtr dataMap);
    void setErrors(MapPtr errorMap);

    template <typename ValueType>
    void addMeta(const std::string& key, ValueType val);
    void addMetadataMap(const MetadataMap& meta);

    void addHeaderText(const std::string& headerText);

    void writeToFile(const std::string& filename);

    //direct string manipulation (no special handling of line comments or other)
    std::string getHeader();
    void addHeaderDirectly(const std::string& headerNew);
private:
    MapPtr data;
    MapPtr errors;
    std::string header;
};

template<typename Key, typename Value>
DataMapWriter<Key, Value>::DataMapWriter() :
    data(), errors(), header("") {
}

template <typename Key, typename Value>
void DataMapWriter<Key,Value>::setData(MapPtr dataMap) {
    data = dataMap;
}

template <typename Key, typename Value>
void DataMapWriter<Key,Value>::setErrors(MapPtr errorMap) {
    errors = errorMap;
}

template <typename Key, typename Value>
template <typename ValueType>
void DataMapWriter<Key,Value>::addMeta(const std::string& key, ValueType val) {
    std::stringstream ss;
    ss << "# " << key << " = " << val << '\n';
    header += ss.str();
}


template <typename Key, typename Value>
void DataMapWriter<Key,Value>::addMetadataMap(const MetadataMap& meta) {
    header += metadataToString(meta, "# ");
}


template <typename Key, typename Value>
void DataMapWriter<Key,Value>::addHeaderText(const std::string& headerText) {
    std::stringstream ss(headerText);
    std::string line;
    while (std::getline(ss, line)) {
        header += "## " + line + '\n';
    }
}

template <typename Key, typename Value>
void DataMapWriter<Key,Value>::writeToFile(const std::string& filename) {
    std::ofstream output(filename.c_str());
    output << header;
    output.precision(15);
    output.setf(std::ios::scientific, std::ios::floatfield);
    if (not errors) {
        for (auto iter = data->cbegin(); iter != data->cend(); ++iter) {
            output << iter->first << '\t' << iter->second << '\n';
        }
    } else {
        for (auto iter = data->cbegin(); iter != data->cend(); ++iter) {
            Key key = iter->first;
            output << key << '\t' << iter->second << '\t' << errors->at(key) << '\n';
        }
    }
}

template<typename Key, typename Value>
inline std::string DataMapWriter<Key, Value>::getHeader() {
    return header;
}

template<typename Key, typename Value>
inline void DataMapWriter<Key, Value>::addHeaderDirectly(
        const std::string& headerNew) {
    header += headerNew;
}

typedef DataMapWriter<double, double> DoubleMapWriter;
typedef DataMapWriter<int, double> IntDoubleMapWriter;
typedef DataMapWriter<std::string, double> StringDoubleMapWriter;




#endif /* DATAMAPWRITER_H_ */
