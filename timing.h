#ifndef FCTRL_TIMING_H
#define FCTRL_TIMING_H

#include <sys/time.h> 

namespace fctrl {

	timeval get_curr_timeval();

	double diff_time(timeval start, timeval end);
}

#endif //FCTRL_TIMING_H
