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
#include <string>
#include <sstream>
#include <iostream>
#include <fstream>
#include <algorithm>
#include <cstdlib>
#include <armadillo>
#include <unistd.h>             // gethostname, getpid, sleep
#include "boost/mpl/assert.hpp" // for not_defined, as introduced below

//instantiate this to get a nice compile-time error message stating that
//something is undefined via BOOST_MPL_ASSERT(( not_defined<T> ));
template<typename T> struct not_defined : boost::mpl::false_ { }; 


//only use with small numbers!
inline uint32_t uint_pow(uint32_t base, uint32_t exponent) {
    uint32_t result = 1;
    for (uint32_t times = 0; times < exponent; ++times) {
        result *= base;
    }
    return result;
} 

template<typename T> inline
std::string numToString(T i) {
    std::stringstream s;
    s << i;
    return s.str();
}

template<typename T> inline
std::string numToString(T i, uint32_t floatPrecision) {
    std::stringstream s;
    s.precision(int(floatPrecision));
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
    for (uint32_t n = 0; n < vec.size(); ++n) {
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


// Wrapper aroundPOSIX glob() that returns a vector<string>
std::vector<std::string> glob(const std::string& path);


// pass this as a default do-nothing function parameter in cases
// where no return value is expected
struct VoidNoOp {
    void operator()() const { };
    template<typename P1, typename... Params>
    void operator()(P1 p1, Params... parameters) {
        (void)(p1);             // we do this just to remove warnings -- we need the recursion for that
        operator()(parameters...);
    }
    // template<class A>
    // void operator()(A a) const { (void)(a); }
    // template<class A, class B>
    // void operator()(A a, B b) const { (void)(a); (void)(b); }
    // template<class A, class B, class C>
    // void operator()(A a, B b, C c) const { (void)(a); (void)(b); (void)(c); }
};



typedef double num;
typedef std::complex<num> cpx;
typedef arma::Mat<num> MatNum;
typedef arma::Mat<cpx> MatCpx;
template<class Matrix>
inline void print_matrix_diff(Matrix mat1, Matrix mat2,
		const std::string& name) {
	//MatNum diff_wrapped = arma::abs(mat1 - mat2) / arma::abs(mat1);
	MatNum diff = arma::abs(mat1 - mat2);
	using arma::max; using arma::mean;
	std::cout << name << ": mean = " << arma::mean(arma::mean(arma::abs(mat1))) << "\n";
	std::cout << name << ": max  diff = " << max(max(diff)) << "\n";
	std::cout << name << ": mean diff = " << mean(mean(diff)) << "\n";
}

template<class Matrix>
inline void print_matrix_rel_diff(Matrix mat1, Matrix mat2,
		const std::string& name) {
	MatNum diff = arma::abs(mat1 - mat2) / arma::abs(mat1);
	using arma::max; using arma::mean;
	std::cout << name << ": max  rel diff = " << max(max(diff))   << "\n";
	std::cout << name << ": mean rel diff = " << mean(mean(diff)) << "\n";
}


template<typename Matrix> inline
void debugSaveMatrix(const Matrix& matrix, const std::string& basename) {
    matrix.eval().save(basename + ".csv", arma::csv_ascii);
}

template<typename Matrix> inline
void debugSaveMatrixCpx(const Matrix& matrix, const std::string& basename) {
    MatNum r = arma::real(matrix);
    r.save(basename + "_real.csv", arma::csv_ascii);
    MatNum i = arma::imag(matrix);
    i.save(basename + "_imag.csv", arma::csv_ascii);
}

template<typename Val> inline
void debugSaveMatrixRealOrCpx(const arma::Mat<Val>& matrix, const std::string& basename) {
    (void)matrix; (void)basename;
}

template<> inline
void debugSaveMatrixRealOrCpx(const arma::Mat<num>& matrix, const std::string& basename) {
    debugSaveMatrix(matrix, basename);
}

template<> inline
void debugSaveMatrixRealOrCpx(const arma::Mat<cpx>& matrix, const std::string& basename) {
    debugSaveMatrixCpx(matrix, basename);
}


// these can be called from the debugger)
void printMatrixReal(const arma::Mat<num>& mat);
void printMatrixComplex(const arma::Mat<cpx>& mat);



// setting the real part of a vector
template<typename VectorInRealPart>
void setVectorReal(arma::Col<cpx>& out, const VectorInRealPart& in) {
    const uint32_t N = out.n_elem;
    for (uint32_t i = 0; i < N; ++i) {
        out[i].real(in[i]);
    }
}

template<typename VectorInRealPart>
void setVectorReal(arma::Col<num>& out, const VectorInRealPart& in) {
    const uint32_t N = out.n_elem;
    for (uint32_t i = 0; i < N; ++i) {
        out[i] = in[i];
    }
}

// setting the imaginary part of a vector
template<typename VectorInImagPart>
void setVectorImag(arma::Col<cpx>& out, const VectorInImagPart& in) {
    const uint32_t N = out.n_elem;
    for (uint32_t i = 0; i < N; ++i) {
        out[i].imag(in[i]);
    }
}

template<typename VectorInImagPart>
void setVectorImag(arma::Col<num>& out, const VectorInImagPart& in) {
    (void)out; (void)in;
}


// transpose a Cube in all three indices.  The (linear) representation
// in memory of the result corresponds is the C-ordered representation
// of a Fortran-ordered original.
template<typename Val>
arma::Cube<Val> transpose_3d(const arma::Cube<Val>& orig) {
    arma::Cube<Val> result(orig.n_slices,
                           orig.n_cols,
                           orig.n_rows); // flipped dimensions
    for (int s = 0; s < orig.n_slices; ++s) {
        for (int r = 0; r < orig.n_rows; ++r) {
            for (int c = 0; c < orig.n_cols; ++c) {
                result(s, c, r) = orig(r, c, s);
            }
        }
    }
    return result;
}



// given a container of numeric values, find the index of the element
// that is numerically nearest to a value to search for.  This does
// *not* require the container to be sorted.
template<typename Container, typename Value>
std::size_t findNearest(const Container& c, Value target) {
    auto it = std::min_element(std::begin(c), std::end(c),
                               [target](Value x, Value y)
                               {
                                   return std::abs(x - target) < std::abs(y - target);
                               });
    return std::distance(std::begin(c), it);
}



#endif /* TOOLS_H_ */
