/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0.  See the enclosed file LICENSE for a copy or if
 * that was not distributed with this file, You can obtain one at
 * http://mozilla.org/MPL/2.0/.
 *
 * Copyright 2017 Max H. Gerlach
 * 
 * */

/*
 * symmat.h
 *
 *  Created on: Oct 26, 2013
 *      Author: max
 */

// simple struct holding data that can be indexed via 2 dimensions --> GenMat

// struct holding data that can be indexed via 2 dimensions, exchange of the two
// indices should yield the same data --> SymMat

// variation without entries for idx1==idx2 --> SymMatOffdiag

#include "checkarray.h"
#include <utility>
#include <stdexcept>


template<typename T, uint32_t size>
struct GenMat {
    static const uint32_t dataSize = size * size;

    checkarray<T, dataSize> data;

    T& operator()(uint32_t idx1, uint32_t idx2) {
#ifdef _GLIBCXX_DEBUG
        if (idx2 >= size or idx1 >= size) {
            throw std::out_of_range("GenMat: access outside matrix");
        }
#endif //_GLIBCXX_DEBUG
        return data[idx1 * size + idx2];
    }
};


template<typename T, uint32_t size>
struct SymMat {
    static const uint32_t dataSize = size * (size + 1) / 2;

    checkarray<T, dataSize> data;

    T& operator()(uint32_t idx1, uint32_t idx2) {
#ifdef _GLIBCXX_DEBUG
        if (idx2 >= size or idx1 >= size) {
            throw std::out_of_range("SymMat: access outside matrix");
        }
#endif //_GLIBCXX_DEBUG
        if (idx2 < idx1) {
            std::swap(idx1, idx2);
        }
        // computation of dataIndex (0-based like everything):
        // dataIndex = size + (size-1) + ... + (size-idx1+1) + idx2 - idx1
        //           = idx1 * size - \sum_{i=1}^{idx1-1} i   + idx2 - idx1
        //           = idx1 * size - (idx1-1)*idx1/2 + idx2 - idx1
        return data[idx1 * size - (idx1-1)*idx1/2 + idx2 - idx1];
    }
};


template<typename T, uint32_t size>
struct SymMatOffdiag {
    static const uint32_t dataSize = size * (size + 1) / 2 - size;  // subtract diagonal length

    checkarray<T, dataSize> data;

    T& operator()(uint32_t idx1, uint32_t idx2) {
#ifdef _GLIBCXX_DEBUG
        if (idx2 == idx1) {
            throw std::out_of_range("SymMatOffdiag: access to diagonal elements");
        }
        if (idx2 >= size or idx1 >= size) {
            throw std::out_of_range("SymMatOffdiag: access outside matrix");
        }
#endif //_GLIBCXX_DEBUG
        if (idx2 < idx1) {
            std::swap(idx1, idx2);
        }
        return data[idx1 * size - (idx1-1)*idx1/2 + idx2 - idx1 - (idx1 + 1)]; // subtract number of elements that would be on diagonal
    }
};
