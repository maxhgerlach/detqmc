/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0.  See the enclosed file LICENSE for a copy or if
 * that was not distributed with this file, You can obtain one at
 * http://mozilla.org/MPL/2.0/.
 *
 * Copyright 2017 Max H. Gerlach
 * 
 * */


#ifndef TOOLSDEBUG_H
#define TOOLSDEBUG_H




// as suggested by
//   http://www.open-mpi.org/faq/?category=debugging
// this can be used to attach to one running MPI process
// via
//   ssh $hostname gdb $executable --pid $pid
#define FREEZE_FOR_DEBUGGER()                                           \
    do {                                                                \
        int i = 0;                                                      \
        char hostname[256];                                             \
        gethostname(hostname, sizeof(hostname));                        \
        std::cout << "PID " << getpid() << " on " << hostname << " ready for attach" << std::endl; \
        while (0 == i)                                                  \
            /*sleep(5)*/;                                               \
    } while(false)

#define CHECK_NAN(mat) for (uint32_t y = 0; y < mat.n_rows; ++y) {      \
        for (uint32_t x = 0; x < mat.n_cols; ++x) {                     \
            if (std::isnan(std::real(mat(y, x))) or                     \
                std::isnan(std::imag(mat(y, x)))) {                     \
                std::cout << #mat << " " << y << " " << x << "\n";      \
                    FREEZE_FOR_DEBUGGER();                              \
            }                                                           \
        }                                                               \
    }

#define CHECK_VEC_NAN(vec) for (uint32_t y = 0; y < vec.n_elem; ++y) {  \
        if (std::isnan(std::real(vec[y])) or                            \
            std::isnan(std::imag(vec[y]))) {                            \
            std::cout << __FILE__ << ":" <<  __LINE__ << " " << #vec << " " << y << "\n"; \
                FREEZE_FOR_DEBUGGER();                                  \
        }                                                               \
    }
 



#endif /* TOOLSDEBUG_H */
