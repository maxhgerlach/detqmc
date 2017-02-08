/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0.  See the enclosed file LICENSE for a copy or if
 * that was not distributed with this file, You can obtain one at
 * http://mozilla.org/MPL/2.0/.
 *
 * Copyright 2017 Max H. Gerlach
 * 
 * */

#if defined (MAX_DEBUG) && ! defined(DUMA_NO_DUMA)
#include "dumapp.h"
#endif

//Code by Max Henner Gerlach, 2010--2012,
//used for the diploma thesis "Directional Ordering in the Classical Compass Model in Two and Three Dimensions"
//contact: maxgerlach@gmail.com

// generalized for SDW DQMC (2015-02-06 - )

/*
 * main-mrpt.cpp
 *
 *  Created on: Jul 11, 2011
 *      Author: gerlach
 */

#include "mrpt-highlevel.h"

int main(int argc, char **argv) {
    initFromCommandLine(argc, argv);

    return 0;
}
