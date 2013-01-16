#ifndef GIT_REVISON_H_
#define GIT_REVISON_H_

extern const char GIT_REVISION_HASH[];
extern const char HOST_NAME[];
extern const char BUILD_DATE[];
extern const char BUILD_TIME[];
extern const char CPPFLAGS[];
extern const char CXXFLAGS[];

#include "metadata.h"
#include "tools.h"
inline MetadataMap collectVersionInfo() {
	MetadataMap meta;
	meta["gitRevisionHash"] = numToString(GIT_REVISION_HASH);
	meta["buildHost"] = numToString(HOST_NAME);
	meta["buildDate"] = numToString(BUILD_DATE);
	meta["buildTime"] = numToString(BUILD_TIME);
	meta["cppflags"] = numToString(CPPFLAGS);
	meta["cxxflags"] = numToString(CXXFLAGS);
	return meta;
}


#endif //GIT_REVISON_H_
