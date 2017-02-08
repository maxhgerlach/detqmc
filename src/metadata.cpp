/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0.  See the enclosed file LICENSE for a copy or if
 * that was not distributed with this file, You can obtain one at
 * http://mozilla.org/MPL/2.0/.
 *
 * Copyright 2017 Max H. Gerlach
 * 
 * */

//Code by Max Henner Gerlach, 2010--2012,
//used for the diploma thesis "Directional Ordering in the Classical Compass Model in Two and Three Dimensions"
//contact: maxgerlach@gmail.com

/*
 * metadata.cpp
 *
 *  Created on: Apr 20, 2011
 *      Author: gerlach
 */

#include "metadata.h"
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wconversion"
#pragma GCC diagnostic ignored "-Wshadow"
#include "boost/property_tree/ptree.hpp"      //config file parsing
#include "boost/property_tree/ini_parser.hpp"
#include "boost/algorithm/string.hpp"   //trimming
#pragma GCC diagnostic pop
#include <sstream>
#include <fstream>
#include <iostream>
#include <vector>

//parse metadata from a configfile-like section of lines of the form
//key = value
//into a map string->string
MetadataMap parseMetadataBlock(std::string& lines) {
    namespace bpt = boost::property_tree;
    bpt::ptree pt;
    std::istringstream linesStream(lines);
    bpt::ini_parser::read_ini(linesStream, pt);
    MetadataMap meta;
    for (bpt::ptree::const_iterator p = pt.begin(); p != pt.end(); ++p) {
        meta[p->first] = p->second.data();
    }
    return meta;
}

std::string metadataToString(const MetadataMap& meta, const char* leading) {
    std::string result;
    MetadataMap::const_iterator p;
    for (p = meta.begin(); p != meta.end(); ++p) {
        result += leading + p->first + " = " + p->second + "\n";
    }
    return result;
}

MetadataMap readOnlyMetadata(const std::string & filename) {
    std::ifstream inFile(filename.c_str());
    std::string collectedLines;
    std::string line;
    bool done = false;
    if (inFile) {
    	do {
            std::getline(inFile, line);
            boost::algorithm::trim_left(line);
            if (line[0] == '#') {
                line[0] = ' ';
                collectedLines += line + '\n';
            } else {
                done = true;
            }
    	} while (inFile and not done);
    } else {
    	std::cerr << "Could not open file " << filename << " for reading.\n";
    	std::cerr << "Error code: " << strerror(errno) << "\n";
    }
    return parseMetadataBlock(collectedLines);
}

void writeOnlyMetaData(const std::string& filename, const MetadataMap& meta,
        const std::string& leadingCommentsBlock,
        bool appendToEndOfFile) {
    std::string comments = leadingCommentsBlock;
    //add "##" characters
    comments.insert(0, "##");
    size_t found = comments.find_first_of('\n');
    while (found != std::string::npos) {
        comments.insert(found + 1, "##");
        found = comments.find_first_of('\n', found + 1);
//        if (found + 1 == comments.length()) {
//            //don't add "##" after the last line
//            break;
//        }
    }
    if (*comments.rbegin() != '\n') {
        comments.append("\n");
    }
//    auto openMode = std::ios::out;
//    if (appendToEndOfFile) {
//      openMode |= std::ios::ate;
//    }
//    std::ofstream outFile(filename.c_str(), openMode);
    std::ofstream outFile;
    if (appendToEndOfFile) {
        outFile.open(filename.c_str(), std::ios::app);
    } else {
        outFile.open(filename.c_str(), std::ios::out);
    }
    if (outFile) {
    	outFile << comments;
    	outFile << metadataToString(meta, "# ");
    	outFile.flush();
    } else {
    	std::cerr << "Could not open file " << filename << " for writing.\n";
    	std::cerr << "Error code: " << strerror(errno) << "\n";
    }
}


MetadataMap getCommonMetadata(const MetadataMap& meta1, const MetadataMap& meta2) {
    MetadataMap result;
        
    typedef MetadataMap::const_iterator Iter;
    Iter it1 = meta1.cbegin();
    Iter it2 = meta2.cbegin();
    while (it1 != meta1.end() and it2 != meta2.end()) {
        // make use of the lexicographic ordering of std::map keys.
        // skip unequal keys:
        if (it1->first < it2->first) {
            ++it1;
        } else if (it2->first < it1->first) {
            ++it2;
        } else {
            // equal keys
            if (it1->second == it2->second) {
                // equal values
                result[it1->first] = it2->second;
            }
            ++it1;
            ++it2;
        }            
    }

    return result;
}
