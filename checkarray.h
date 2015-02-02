/*
 * checkarray.h
 *
 *  Created on: May 22, 2013
 *      Author: gerlach
 */

#ifndef CHECKARRAY_H_
#define CHECKARRAY_H_

// a variation of std::array that includes bounds checking for
// operator[] if the compilation flag -D_GLIBCXX_DEBUG is given (like
// std::vector), compare to:
// http://stackoverflow.com/questions/8181887/bound-checking-of-stdarray-in-debug-version-of-gcc

#include <array>
#include <initializer_list>

#ifdef _GLIBCXX_DEBUG
template<typename T, std::size_t N>
struct checkarray : public std::array<T,N> {
    typedef std::array<T,N> Base;
    T& operator[](std::size_t n) {
        return Base::at(n);
    }

    const T& operator[](std::size_t n) const {
        return Base::at(n);
    }

    //std::array<T,N> is an aggregate type, but because we
    //derive from it, checkarray<T,N> no longer is --> need
    //to provide explicit constructors
    checkarray(std::initializer_list<T> l) {
        //this is not very pretty and lacks compile time
        //bounds checking, but I have not found an easy solution
        if (l.size() != Base::size()) {
            throw std::out_of_range("checkarray.h: Size of initializer_list does not match size of array");
        }
        std::copy(std::begin(l), std::end(l), this->begin());
    }

    checkarray() : Base()
    { }
};
#else //_GLIBCXX_DEBUG
//This is an ugly workaround for lacking support for template-aliases
#define checkarray std::array

#endif //_GLIBCXX_DEBUG




#endif /* CHECKARRAY_H_ */
