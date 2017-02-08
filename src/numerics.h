/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0.  See the enclosed file LICENSE for a copy or if
 * that was not distributed with this file, You can obtain one at
 * http://mozilla.org/MPL/2.0/.
 *
 * Copyright 2017 Max H. Gerlach
 * 
 * */

/*
 * numerics.h
 *
 *  Created on: Oct 25, 2011
 *      Author: gerlach
 */

// update for SDW DQMC (2015-02-06 - )


#ifndef NUMERICS_H_
#define NUMERICS_H_

#include <iostream>
#include <cmath>
#include <limits>
#include "tools.h"
#include "exceptions.h"

//find brentMinimize<f> finds the minimum of a function (callable
//object) f. It searches intervalStart <= x <= intervalEnd.
//Afterwards: @minLocation = x_min, @minValue = f(x_min).
//eps -- accuracy of evaluation
//The Brent algorithm is used to keep the number of evaluations of f low.
template <typename Function>
void brentMinimize(double& minLocation, double& minValue, Function& f,
        double intervalStart, double intervalEnd,
        double tolerance = 3.0e-8) {
    using namespace std;
    const double eps = std::numeric_limits<double>::epsilon() * 1.0e-3;

    const double goldenRatio = 0.3819660;       //(3 - sqrt(5))/2
    //minimum between a, b:
    double a = intervalStart;
    double b = intervalEnd;
    double x;   //current minimum
    double v;   //second best point
    double w;   //previous value of w
    v = w = x = a + goldenRatio * (b - a);
    double u;   //most recent point of evaluation
    //function evaluated at x,v,w
    double fx, fv, fw, fu;
    fv = fw = fx = f(x);
    double tol1, tol2;         //minimal relative movement in x
    double m;       //mid point of a and b
    double d = 0;       //distance moved in last step
    double e = 0;       //distance moved in before last step
    do {
        m = 0.5 * (a + b);
//        tol = eps * abs(x) + eps / 4;       //Brent: parameter t instead of eps/4
        tol1 = tolerance * abs(x) + eps;
        tol2 = 2 * tol1;
        //Check stopping criterion max(x-a, b-x) <= 2*tol
        if (abs(x - m) <= tol2 - (b - a) / 2) break;
        if (abs(e) > tol1) {
            //fit parabola in p, q, r
            double r = (x - w) * (fx - fv);
            double q = (x - v) * (fx - fw);
            double p = (x - v) * q - (x - w) * r;
            q = 2 * (q - r);
            if (q > 0) {
                p = -p;
            } else {
                q = -q;
            }
            double eLast = e;           //assigned to r in Brent
            e = d;
            //changed relation signs in following line:
            if (abs(p) >= abs(0.5 * q * eLast) or p <= q * (a - x) or p >= q * (b - x)) {      //TODO:check signs
                //can't do parabolic fit, do a golden section step instead:
                e = (x < m) ? (b - x) : (a - x);
                d = goldenRatio * e;
            } else {
                //parabolic interpolation step
                d = p / q;
                u = x + d;
                //don't evaluate f to close to a or b
                if (u - a < tol2 or b - u < tol2) {
                    d = (x < m) ? tol1 : -tol1;
                }
            }
        } else {
            //golden section
            e = (x < m) ? (b - x) : (a - x);
            d = goldenRatio * e;
        }
        //update current position, don't evaluate f too close to x
        u = (abs(d) > tol1) ? (x + d) : ((d > 0) ? (x + tol1) : (x - tol1));
        fu = f(u);
        //update a,b,v,w,x
        if (fu <= fx) {
            //accept new point!
            if (u < x) {
                b = x;
            } else {
                a = x;
            }
            v = w;
            fv = fw;
            w = x;
            fw = fx;
            x = u;
            fx = fu;
        } else {
            //new point is worse, but must be better than one of a or b
            if (u < x) {
                a = u;
            } else {
                b = u;
            }
            if (fu <= fw or w == x) {
                v = w;
                fv = fw;
                w = u;
                fw = fu;
            } else if (fu <= fv or v ==x or v == w) {
                v = u;
                fv = fu;
            }
        }
    } while (true);         //TODO: enforce max number of iterations?
    //finished:
    minLocation = x;
    minValue = fx;
}


//use a combination of bisection and Newton-Raphson to find a root x0 of
//func(x), uses the derivative dfunc as well
//pass an interval [xMin, xMax] within which x0 is to be found
//convergence condition: abs(func(x0)) < acc (this is different
//from normal root finding -- adapted to the CEI problem)
//this is similar to Numerical Recipes 3rd editon, Ch 9.4, rtsafe()
template<typename Function, typename Derivative>
double findRoot(Function& func, Derivative& dfunc,
        double xMin, double xMax, double acc,
        int maxIterations = 100, bool verbose = true) {
    double xh, xl;          //"high" and "low" limits of bracket
    double fl = func(xMin);
    if (fl == .0) return xMin;
    double fh = func(xMax);
    if (fh == .0) return xMax;
    if ((fl > .0 and fh > .0) or (fl < .0 and fh < .0)) {
        throw GeneralError("Root is not bracketed between "
                           + numToString(xMin) + " and " + numToString(xMax));
    }
    //find direction of search:
    if (fl <.0) {
        xl = xMin;
        xh = xMax;
    } else {
        xh = xMin;
        xl = xMax;
    }
    double root = .5 * (xMin + xMax);       //initial guess
    double dxOld = std::fabs(xMax - xMin);
    double dx = dxOld;
    double f = func(root);
    double df = dfunc(root);

    for (int j = 0; j < maxIterations; ++j) {
        bool finished = false;
        if ( (((root-xh)*df - f) * ((root-xl)*df - f) > .0)
              or (std::fabs(2.0*f) > std::fabs(dxOld*df)) ) {
            //Newton is out of range or not decreasing fast enough
            // -> bisection step
            dxOld = dx;
            dx = .5 * (xh - xl);
            root = xl + dx;
//            if (xl == root) finished = true;    //negligible change -> finished
        } else {
            //Newton step can be taken
            dxOld = dx;
            dx = f / df;
            double temp = root;
            root -= dx;
//            if (temp == root) finished = true;  //negligible change -> finished
        }
        //new function evaluation:
        f = func(root);
        df = dfunc(root);

        if ((j > 0) and (std::fabs(f) < acc)) finished = true;    //converged
        if (verbose) {
            std::cout << " [" << xl << ", " << xh << "]: "
                    "f(" << root << ")=" << f << std::endl;
        }
        if (finished) {
            return root;
        }

        //update bracket:
        if (f < .0) {
            xl = root;
        } else {
            xh = root;
        }
    }

    throw GeneralError(numToString(maxIterations) + " iterations did"
                       "not suffice to find root with relative accuracy " +
                       numToString(acc) + ", last estimate: " + numToString(root)
                       + ", accuracy " + numToString(std::fabs(dx / root)));
}


#endif /* NUMERICS_H_ */
