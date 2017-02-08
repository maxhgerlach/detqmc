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
 * logval.h
 *
 *  Created on: Jan 12, 2011
 *      Author: gerlach
 */

#ifndef LOGVAL_H_
#define LOGVAL_H_

#include <cmath>
#include <iostream>
#include <vector>

//calculations with LogVal's internally use the natural logarithm of values
//to avoid overflows
class LogVal {
public:
    static constexpr double LogZero = 1e-3; //ln(LogZero) << ln(1)
    double lnx;

    LogVal();
    explicit LogVal(double x);
    LogVal(const LogVal& lv);

    LogVal& addSmaller(LogVal smaller); //add a smaller number
    LogVal& addLarger(LogVal larger); //add a larger number
    LogVal& subtractSmaller(LogVal smaller);
    LogVal& subtractLarger(LogVal larger);

    LogVal& operator *=(LogVal rhs);
    LogVal& operator /=(LogVal rhs);
    LogVal& operator +=(LogVal rhs); //these use an additional check as compared to addSmaller and addLarger
    LogVal& operator -=(LogVal rhs); //of course only for rhs < lhs!

    LogVal& operator *=(double rhs);
    LogVal& operator /=(double rhs);

    bool operator==(LogVal rhs) const;
    bool operator!=(LogVal rhs) const;
    bool operator<=(LogVal rhs) const;
    bool operator>=(LogVal rhs) const;
    bool operator<(LogVal rhs) const;
    bool operator>(LogVal rhs) const;
};

/*
double toDouble(LogVal val);        //returns exp(lnx)
double toDouble(double val);        //does nothing to val

LogVal toLogValExp(double exponent);    //returns a LogVal representing exp(exponent), i.e. set lnx to exponent

LogVal operator+(LogVal a, LogVal b);
LogVal operator-(LogVal a, LogVal b);
LogVal operator*(LogVal a, LogVal b);
LogVal operator/(LogVal a, LogVal b);
LogVal operator*(LogVal a, double b);
LogVal operator/(LogVal a, double b);

LogVal pow(LogVal base, double exponent);

std::ostream& operator<<(std::ostream& stream, LogVal val);
*/

inline bool LogVal::operator==(LogVal rhs) const {
    return lnx == rhs.lnx;
}
inline bool LogVal::operator!=(LogVal rhs) const {
    return lnx != rhs.lnx;
}
inline bool LogVal::operator<=(LogVal rhs) const {
    return lnx <= rhs.lnx;
}
inline bool LogVal::operator>=(LogVal rhs) const {
    return lnx >= rhs.lnx;
}
inline bool LogVal::operator<(LogVal rhs) const {
    return lnx < rhs.lnx;
}
inline bool LogVal::operator>(LogVal rhs) const {
    return lnx > rhs.lnx;
}

inline LogVal::LogVal() {
//  lnx = LogZero;
}

inline LogVal::LogVal(double x) {
    lnx = std::log(x);
}

inline LogVal::LogVal(const LogVal& lv) {
    lnx = lv.lnx;
}

inline LogVal& LogVal::addSmaller(LogVal smaller) {
    lnx += log1p(std::exp(smaller.lnx - lnx));
    return *this;
}

inline LogVal& LogVal::subtractSmaller(LogVal smaller) {
    lnx += log1p(-std::exp(smaller.lnx - lnx));
    return *this;
}

inline LogVal& LogVal::addLarger(LogVal larger) {
    lnx = larger.lnx + log1p(std::exp(lnx - larger.lnx));
    return *this;
}

inline LogVal& LogVal::subtractLarger(LogVal larger) {
    lnx = larger.lnx + log1p(-std::exp(lnx - larger.lnx));
    return *this;
}


inline LogVal& LogVal::operator *=(LogVal rhs) {
    lnx += rhs.lnx;
    return *this;
}

inline LogVal& LogVal::operator /=(LogVal rhs) {
    lnx -= rhs.lnx;
    return *this;
}

inline LogVal& LogVal::operator +=(LogVal rhs) {
//  if (lnx == std::log(0)) {
//      lnx = rhs.lnx;
//  } else
    if (rhs.lnx <= lnx) {
        addSmaller(rhs);
    } else {
        addLarger(rhs);
    }
    return *this;
}

inline LogVal& LogVal::operator -=(LogVal rhs) {
    if (rhs.lnx <= lnx) {
        subtractSmaller(rhs);
    } else {
        subtractLarger(rhs);
    }
    return *this;
}

inline LogVal& LogVal::operator *=(double rhs) {
    lnx += std::log(rhs);
    return *this;
}

inline LogVal& LogVal::operator /=(double rhs) {
    lnx -= std::log(rhs);
    return *this;
}

inline LogVal operator+(LogVal a, LogVal b) {
    LogVal t = a;
    return t += b;
}
inline LogVal operator-(LogVal a, LogVal b) {
    LogVal t = a;
    return t -= b;
}
inline LogVal operator*(LogVal a, LogVal b) {
    LogVal t = a;
    return t *= b;
}
inline LogVal operator/(LogVal a, LogVal b) {
    LogVal t = a;
    return t /= b;
}
inline LogVal operator*(LogVal a, double b) {
    LogVal t = a;
    return t *= b;
}
inline LogVal operator/(LogVal a, double b) {
    LogVal t = a;
    return t /= b;
}

inline LogVal pow(LogVal base, double exponent) {
    LogVal r;
    r.lnx = base.lnx * exponent;
    return r;
}

inline std::ostream& operator<<(std::ostream& stream, LogVal val) {
//  return stream << std::exp(val.lnx);
//  return stream << "exp(" << val.lnx << ")";
    return stream << val.lnx;
}

inline std::istream& operator>>(std::istream& stream, LogVal& val) {
//  return stream << std::exp(val.lnx);
//  return stream << "exp(" << val.lnx << ")";

    return stream >> val.lnx;
}

inline double toDouble(LogVal val) {
    return std::exp(val.lnx);
}

inline double toDouble(double val) {
    return val;
}

inline LogVal toLogValExp(double exponent) {
    LogVal temp;
    temp.lnx = exponent;
    return temp;
}


//numerically stable addition of a list of logarithmic values:
inline LogVal logSum(const std::vector<LogVal>& logvals) {
    LogVal result = logvals[0];
    for (unsigned n = 1; n < logvals.size(); ++n) {
        LogVal l = logvals[n];
        if (l < result) {
            result.addSmaller(l);
        } else {
            result.addLarger(l);
        }
    }
    return result;
}

#endif /* LOGVAL_H_ */
