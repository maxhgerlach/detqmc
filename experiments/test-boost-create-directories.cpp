#include "boost/filesystem.hpp"

int main() {
    namespace fs = boost::filesystem;

    // the following also works repeatedly, i.e. if the path exists already
    fs::create_directories(fs::path("test"));

    // this fails at run time with an exception
    fs::create_directories(fs::path("/test"));

    return 0;
}
