/*
 * boost_serialize_uniqueptr.h
 *
 *  Created on: Apr 16, 2013
 *      Author: gerlach
 */

#ifndef BOOST_SERIALIZE_UNIQUEPTR_H_
#define BOOST_SERIALIZE_UNIQUEPTR_H_

#include <memory>

//Boost serialization for std::unique_ptr<T>, following
// http://stackoverflow.com/questions/12915267/boost-serialization-stdunique-ptr-support

namespace boost { namespace serialization {

template<class Archive, class T>
inline void save(Archive& ar, const std::unique_ptr<T>& t,
                 const uint32_t file_version) {
    // only the raw pointer has to be saved
    const T* const tx = t.get();
    ar << tx;
}

template<class Archive, class T>
inline void load(Archive & ar, std::unique_ptr<T>& t,
                 const uint32_t file_version) {
    T* pTarget;
    ar >> pTarget;

    #if BOOST_WORKAROUND(BOOST_DINKUMWARE_STDLIB, == 1)
        t.release();
        t = std::unique_ptr<T>(pTarget);
    #else
        t.reset(pTarget);
    #endif
}

// split non-intrusive serialization function into separate
// non intrusive save/load functions
template<class Archive, class T>
inline void serialize(Archive& ar, std::unique_ptr<T>& t,
                      const uint32_t file_version) {
    boost::serialization::split_free(ar, t, file_version);
}

} } // namespace boost::serialization



#endif /* BOOST_SERIALIZE_UNIQUEPTR_H_ */
