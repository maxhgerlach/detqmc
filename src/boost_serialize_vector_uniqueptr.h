/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0.  See the enclosed file LICENSE for a copy or if
 * that was not distributed with this file, You can obtain one at
 * http://mozilla.org/MPL/2.0/.
 *
 * Copyright 2017 Max H. Gerlach
 * 
 * */

/*
 * boost_serialize_vector_uniqueptr.h
 *
 *  Created on: Apr 16, 2013
 *      Author: gerlach
 */

#ifndef BOOST_SERIALIZE_VECTOR_UNIQUEPTR_H_
#define BOOST_SERIALIZE_VECTOR_UNIQUEPTR_H_

//Boost serialization for std::vector<std::unique_ptr<T>>, based on
// http://stackoverflow.com/questions/13347776/boost-serialization-of-stl-collection-of-std-unique-ptrs

//NOTE: do not include boost/serialization/vector.hpp //is this note still relevant?

#include <memory>
#include <vector>
#include <boost/serialization/nvp.hpp>

namespace boost { namespace serialization {


template<class Archive, class T, class Allocator>
inline void save(Archive& ar,
                 const std::vector<std::unique_ptr<T>, Allocator>& vec,
                 const uint32_t /*version*/) {
    collection_size_type count = vec.size();
    ar << BOOST_SERIALIZATION_NVP(count);
    for (auto it = vec.begin(), end = vec.end(); it != end; ++it)
        ar << boost::serialization::make_nvp("item", (*it));
}

template<class Archive, class T, class Allocator>
inline void load(Archive& ar,
                 std::vector<std::unique_ptr<T>, Allocator>& vec,
                 const uint32_t /*version*/) {
    collection_size_type count;
    ar >> BOOST_SERIALIZATION_NVP(count);
    vec.clear();
    vec.reserve(count);
    while( count-- > 0 ) {
        std::unique_ptr<T> i;
        ar >> boost::serialization::make_nvp("item", i);
        vec.push_back(std::move(i)); // use std::move
    }
}

template<class Archive, class T, class Allocator>
inline void serialize(Archive& ar,
                      std::vector<std::unique_ptr<T>, Allocator>& t,
                      const uint32_t file_version) {
    boost::serialization::split_free(ar, t, file_version);
}


} } // namespace boost::serialization


#endif /* BOOST_SERIALIZE_VECTOR_UNIQUEPTR_H_ */
