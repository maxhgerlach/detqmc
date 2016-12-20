#if defined (MAX_DEBUG) && ! defined(DUMA_NO_DUMA)
#include "dumapp.h"
#endif

#include "mrpt-highlevel.h"

int main(int argc, char **argv) {
    initFromCommandLineBC(argc, argv);

    return 0;
}
