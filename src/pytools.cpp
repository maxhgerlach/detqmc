/*
 * pytools.cpp
 *
 *  Created on: Mar 27, 2014
 *      Author: max
 */


#include <cstdio>
#include <string>
#include <python2.7/Python.h>
#include <armadillo>

void python_matshow(const arma::mat& mat, const std::string& title) {
	std::string tempfile = std::tmpnam(0);
	tempfile += ".csv";
	mat.save(tempfile, arma::csv_ascii);

	static bool initialized = false;

	if (not initialized) {
		Py_Initialize();
	}
    PyRun_SimpleString("import matplotlib.pyplot as plt");
    PyRun_SimpleString("import numpy as np");
    std::string loadcommand = "mat = np.loadtxt(\"" + tempfile + "\", delimiter=',')";
    PyRun_SimpleString(loadcommand.c_str());
    PyRun_SimpleString("plt.matshow(mat); plt.colorbar()");
    PyRun_SimpleString(("plt.title(r\"" + title + "\")").c_str());
    PyRun_SimpleString("plt.show()");
    std::remove(tempfile.c_str());
}

void python_matshow2(const arma::mat& mat1, const std::string& title1,
		const arma::mat& mat2, const std::string& title2) {
	std::string tempfile1 = std::tmpnam(0);
	tempfile1 += ".csv";
	mat1.save(tempfile1, arma::csv_ascii);
	std::string tempfile2 = std::tmpnam(0);
	tempfile2 += ".csv";
	mat2.save(tempfile2, arma::csv_ascii);

	static bool initialized = false;

	if (not initialized) {
		Py_Initialize();
	}
    PyRun_SimpleString("import matplotlib.pyplot as plt");
    PyRun_SimpleString("import numpy as np");

    std::string loadcommand1 = "mat1 = np.loadtxt(\"" + tempfile1 + "\", delimiter=',')";
    PyRun_SimpleString(loadcommand1.c_str());
    PyRun_SimpleString("plt.matshow(mat1); plt.colorbar()");
    PyRun_SimpleString(("plt.title(r\"" + title1 + "\")").c_str());

    std::string loadcommand2 = "mat2 = np.loadtxt(\"" + tempfile2 + "\", delimiter=',')";
    PyRun_SimpleString(loadcommand2.c_str());
    PyRun_SimpleString("plt.matshow(mat2); plt.colorbar()");
    PyRun_SimpleString(("plt.title(r\"" + title2 + "\")").c_str());

    PyRun_SimpleString("plt.show()");

    std::remove(tempfile1.c_str());
    std::remove(tempfile2.c_str());
}



