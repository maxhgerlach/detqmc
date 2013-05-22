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

template<typename T, std::size_t N>
struct checkarray : std::array<T,N> {
#ifdef _GLIBCXX_DEBUG
	T& operator[](std::size_t n)	{
		return at(n);
	}

	const T& operator[](std::size_t n) const	{
		return at(n);
	}
#endif //_GLIBCXX_DEBUG
};




#endif /* CHECKARRAY_H_ */
