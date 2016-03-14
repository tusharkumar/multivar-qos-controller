#include <iostream>
#include <cassert>

#include "debuglog.h"
#include "timing.h"

namespace fctrl {

timeval get_curr_timeval() {
	struct timeval tv;

	gettimeofday(&tv, NULL);
	return tv;
}

double diff_time(timeval start, timeval end) {
	timeval diff;

#ifdef FC_TIMING_DEBUG
	std::cout << "diff_time: "
		<< "start = (" << start.tv_sec << ", " << start.tv_usec << ") "
		<< "end = (" << end.tv_sec << ", " << end.tv_usec << ") ";
#endif //FC_TIMING_DEBUG

	assert(end.tv_sec >= start.tv_sec);

	diff.tv_sec = end.tv_sec - start.tv_sec;
	if(end.tv_usec >= start.tv_usec) {
		diff.tv_usec = end.tv_usec - start.tv_usec;
	}
	else {
		diff.tv_sec--;
		diff.tv_usec = end.tv_usec + (1000000 - start.tv_usec);
	}

	double diff_val = double(diff.tv_sec) + double(diff.tv_usec)/1000000;
#ifdef FC_TIMING_DEBUG
	std::cout << " diff = " << diff_val << std::endl;
#endif //FC_TIMING_DEBUG

	return diff_val; //in seconds
}

} //namespace fctrl
