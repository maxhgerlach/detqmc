//Code by Max Henner Gerlach, 2010--2012,
//used for the diploma thesis "Directional Ordering in the Classical Compass Model in Two and Three Dimensions"
//contact: maxgerlach@gmail.com

/*
 * tools.h
 *
 *  Created on: Feb 10, 2011
 *      Author: gerlach
 */

#ifndef TOOLS_H_
#define TOOLS_H_

#include <vector>
#include <sstream>
#include <iostream>
#include <fstream>
#include <algorithm>
#include <cstdlib>

template<typename T> inline
std::string numToString(T i) {
    std::stringstream s;
    s << i;
    return s.str();
}

template<typename T> inline
std::string numToString(T i, unsigned floatPrecision) {
    std::stringstream s;
    s.precision(floatPrecision);
    s.setf(std::ios::scientific, std::ios::floatfield);
    s << i;
    return s.str();
}


template<typename T> inline
T fromString(const std::string& str) {
    std::stringstream s(str);
    T i;
    s >> i;
    return i;
}

template<typename T> inline
T fromString(const char* str) {
    return fromString<T>(std::string(str));
}




template<typename T> inline
void destroy(T*& pointer) {
    delete pointer;
    pointer = 0;
}

template<typename T> inline
void destroyAll(std::vector<T*>& vec) {
    for (typename std::vector<T*>::iterator iter = vec.begin(); iter != vec.end(); ++iter) {
        delete *iter;
        *iter = 0;
    }
}

template<typename T> inline
void printVector(const std::vector<T>& vec) {
    for (unsigned n = 0; n < vec.size(); ++n) {
        std::cout << vec[n] << "\t";
    }
    std::cout << std::endl;
}



//compare two entries of a std::map by their mapped values using this function
template<typename K, typename V> inline
bool mapValueCompare(const std::pair<K, V>& a, const std::pair<K, V>& b) {
    return (a.second < b.second);
}



template <typename BoostArray, typename Value>
void initArray(BoostArray& array, Value val) {
    std::fill(array.data(), array.data() + array.num_elements(), val);
}
//
//template <typename BoostArray, typename STLVector>
//void copyVectorToArray(BoostArray& array, const STLVector& vector) {
//  std::copy(vector.begin(), vector.end(), array.data());
//}

//rename a file to "$filename~"
//it is not a good idea to do this if we write to the original file name
//directly after the call (moving seems to take some time to have effect...)
inline void moveFileBackUp(const std::string& filename) {
    std::string command = "mv " + filename + " " + filename + "~" + " &>/dev/null";
    std::system(command.c_str());
}

//for small files:
inline void copyFileBackUp(const std::string& filename) {
    std::ifstream input(filename.c_str(), std::ios::binary);
    std::ofstream output((filename + "~").c_str(), std::ios::trunc | std::ios::binary);
    output << input.rdbuf();
}



#endif /* TOOLS_H_ */
