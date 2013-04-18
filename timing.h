/*
 * timing.h
 *
 *  Created on: Mar 28, 2013
 *      Author: gerlach
 */

#ifndef TIMING_H_
#define TIMING_H_

#ifdef TIMING

#include <string>
#include <map>
#include <iostream>
#include "boost/timer/timer.hpp"


class Timing {
public:
	Timing() : timers() {
	}
	//on destruction print timing summary
	~Timing() {
		std::cout << "\nTimings:\n";
		for (const auto& kv : timers) {
			std::cout << kv.first << ": " << kv.second.format();
		}
	}
	//start or resume a timer:
	void start(const std::string& timerKey) {
		if (timers.count(timerKey) == 0) {
			timers[timerKey] = boost::timer::cpu_timer();
		}
		timers[timerKey].resume();
	}
	void stop(const std::string& timerKey) {
		timers.at(timerKey).stop();
	}
private:
	std::map<std::string,boost::timer::cpu_timer> timers;
};

#else
//version that does nothing:
#include <string>
class Timing {
public:
	Timing() {
	}
	~Timing() {
	}
	void start(const std::string& timerKey) {
		(void)timerKey;
	}
	void stop(const std::string& timerKey) {
		(void)timerKey;
	}
};

#endif //TIMING

extern Timing timing;		//one global timing object, defined in timing.cpp



#endif /* TIMING_H_ */
