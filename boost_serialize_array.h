/*
 * boost_serialize_array.h
 *
 *  Created on: Apr 17, 2013
 *      Author: max
 */

#ifndef BOOST_SERIALIZE_ARRAY_H_
#define BOOST_SERIALIZE_ARRAY_H_

//serialization support for std::array

#include <array>
#include <boost/serialization/array.hpp>

template <class Archive, class T, std::size_t N>
void serialize(Archive& ar, std::array<T,N>& a, const unsigned int /* version */) {
	ar & boost::serialization::make_array(a.data(), a.size());
}


#endif /* BOOST_SERIALIZE_ARRAY_H_ */
