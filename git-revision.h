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
	meta["git revision hash"] = numToString(GIT_REVISION_HASH);
	meta["build host"] = numToString(HOST_NAME);
	meta["build date"] = numToString(BUILD_DATE);
	meta["build time"] = numToString(BUILD_TIME);
	meta["cppflags"] = numToString(CPPFLAGS);
	meta["cxxflags"] = numToString(CXXFLAGS);
	return meta;
}


#endif //GIT_REVISON_H_
