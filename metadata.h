//Code by Max Henner Gerlach, 2010--2012,
//used for the diploma thesis "Directional Ordering in the Classical Compass Model in Two and Three Dimensions"
//contact: maxgerlach@gmail.com

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

typedef std::map<std::string, std::string> MetadataMap;

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



#endif /* METADATA_H_ */
