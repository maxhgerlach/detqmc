/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0.  See the enclosed file LICENSE for a copy or if
 * that was not distributed with this file, You can obtain one at
 * http://mozilla.org/MPL/2.0/.
 *
 * Copyright 2017 Max H. Gerlach
 * 
 * */

/*
 * metadata.h
 *
 *  Created on: Apr 20, 2011
 *      Author: gerlach
 */

#ifndef METADATA_H_
#define METADATA_H_

#include <map>
#include <string>
#include <sstream>
#include <iterator>             // std::iterator_traits
#include <type_traits>          // std::is_same
#include "exceptions.h"

typedef std::map<std::string, std::string> MetadataMap;

template<typename Value>
void getMeta(const MetadataMap& meta, const std::string& key, Value& valOut) {
    if (meta.count(key) == 0) {
        throw_KeyUndefined(key);
    } else {
        std::stringstream parseStream(meta.at(key));
        parseStream >> valOut;
    }
}

//parse metadata from a configfile-like section of lines of the form
//key = value
//into a map string->string
MetadataMap parseMetadataBlock(std::string& lines);

//open the file and read only its leading comment block, parse that and
//return the metadata
MetadataMap readOnlyMetadata(const std::string& filename);

//write out metadata to a string in the form
//[leading] key = value
//with [leading] = "#" for instance
std::string metadataToString(const MetadataMap& meta,
        const char* leading = "");

//write out metadata to file, with leading "#" characters;
//additional block of comment lines gets written out first with "##" added
//before each line
void writeOnlyMetaData(const std::string& filename, const MetadataMap& meta,
        const std::string& leadingCommentsBlock = "",
        bool appendToEndOfFile = false);

// return a map with the (key, value) pairs that are equal in meta1 and meta2
MetadataMap getCommonMetadata(const MetadataMap& meta1, const MetadataMap& meta2);

// return a map with the (key, value) pairs that are equal in all the MetadataMaps
// in the range pointed to by the Iter iterators
template<typename Iter>
typename std::enable_if<std::is_same<typename std::iterator_traits<Iter>::value_type,
                                     MetadataMap>::value,
                        MetadataMap>::type
getCommonMetadata(Iter curr, Iter end) {
    MetadataMap commonMetadata = *curr;
    ++curr;
    while (curr != end) {
        commonMetadata = getCommonMetadata(commonMetadata, *curr);
        ++curr;
    }
    return commonMetadata;
}

#endif /* METADATA_H_ */
