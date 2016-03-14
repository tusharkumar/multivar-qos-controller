#ifndef __SAMPLE_STATS_H__
#define __SAMPLE_STATS_H__

#include "samplestats_internal.h"

namespace fctrl {

void add_sample(
		const std::vector<double>& v_x, // sample to add
		std::vector<double>& mean,      // input: mean, std, maxswing, wmin and wmax before adding sample
		std::vector<double>& std,       // output: mean, std, maxswing, wmin and wmax after adding sample
		std::vector<double>& maxswing,
		std::vector<double>& wmin,
		std::vector<double>& wmax,
		std::size_t len,                // length of sample sequence prior to adding v_x
		double f                        // forget rate
);

template<typename HistoryIndexer>
void remove_sample(
		const HistoryIndexer& h,        // h[i] must give std::vector<double> for i(th) past frame (i=0, 1, ..., h.size()-1)
		std::vector<double>& mean,      // input: mean, std, maxswing, wmin, wmax before removing weighted sample
		std::vector<double>& std,       // output: mean, std, maxswing, wmin, wmax after removing weighted sample
		std::vector<double>& maxswing,
		std::vector<double>& wmin,
		std::vector<double>& wmax,
		std::size_t len,                // effective length of history prior to removing oldest weighted sample
		double f                        // forget rate
);


} //namespace fctrl

#endif //__SAMPLE_STATS_H__
