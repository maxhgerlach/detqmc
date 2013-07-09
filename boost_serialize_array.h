/*
 * boost_serialize_array.h
 *
 *  Created on: Apr 17, 2013
 *      Author: max
 */

#ifndef BOOST_SERIALIZE_ARRAY_H_
#define BOOST_SERIALIZE_ARRAY_H_

//serialization support for std::array or checkarray

#include <array>
#include "checkarray.h"
#include "boost/serialization/array.hpp"

namespace boost { namespace serialization {

template <class Archive, class T, std::size_t N>
void serialize(Archive& ar, checkarray<T,N>& a, const uint32_t /* version */) {
    ar & boost::serialization::make_array(a.data(), a.size());
}

} } //namespace boost::serialization


#endif /* BOOST_SERIALIZE_ARRAY_H_ */
