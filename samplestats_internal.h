#ifndef __SAMPLE_STATS__INTERNAL_H__
#define __SAMPLE_STATS__INTERNAL_H__

#include <cassert>
#include <cmath>
#include <vector>

namespace fctrl {


void remove_sample_std(
		const std::vector<double>& v_x, // sample to remove
		std::vector<double>& mean,      // input: mean, std before removing weighted sample
		std::vector<double>& std,       // output: mean, std after removing weighted sample
		std::size_t len,                // length of sample sequence prior to removing weighted sample
		double f                        // forget rate
);


static inline bool approx_eq(double val1, double val2) {
	double diff = std::abs(val1 - val2);
	double max  = std::max( std::abs(val1), std::abs(val2) );
	return (diff < 0.01 * max);
}

template<typename HistoryIndexer, typename ReduceOp>
static inline
double weighted_reduce(
		const HistoryIndexer& h,
		size_t                len,
		size_t                index,
		ReduceOp              op,
		double                f
)
{
	double r = (h[0])[index];
	double w = f;
	for(size_t i=1; i<len; i++) {
		r = op(r, w * (h[i])[index]);
		w *= f;
	}
	return r;
}

template<typename HistoryIndexer>
void remove_sample_maxswing(
		const HistoryIndexer& h,
		std::vector<double>& maxswing,  // input: maxswing, wmin, wmax before removing weighted sample
		std::vector<double>& wmin,      // output: maxswing, wmin, wmax after removing weighted sample
		std::vector<double>& wmax,
		std::size_t len,                // effective length of history prior to removing oldest weighted sample
		double f                        // forget rate
)
{
	assert(len > 0);
	assert(len <= h.size());

	const std::vector<double>& v_x = h[len-1]; // sample to remove

	assert(maxswing.size() == v_x.size());
	assert(wmin.size() == v_x.size());
	assert(wmax.size() == v_x.size());

	if(len > 1) {
		double f_to_len_minus_1 = std::pow(f, len-1);
		std::vector<double> wv_x = v_x;
		for(size_t i=0; i<wv_x.size(); i++)
			wv_x[i] *= f_to_len_minus_1;

		for(size_t i=0; i<maxswing.size(); i++) {
			if( approx_eq(wmin[i], wv_x[i]) )
				wmin[i] = weighted_reduce(h, len-1, i, std::min<double>, f);
			if( approx_eq(wmax[i], wv_x[i]) )
				wmax[i] = weighted_reduce(h, len-1, i, std::max<double>, f);
			maxswing[i] = wmax[i] - wmin[i];
		}
	}
	else {
		for(std::size_t i=0; i<maxswing.size(); i++) {
			maxswing[i] = 0.0;
			wmin[i] = 0.0;
			wmax[i] = 0.0;
		}
	}
}

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
)
{
	assert(len > 0);
	assert(len <= h.size());

	const std::vector<double>& v_x = h[len-1]; // sample to remove

	remove_sample_std(v_x, mean, std, len, f);
	remove_sample_maxswing(h, maxswing, wmin, wmax, len, f);
}

} //namespace fctrl

#endif //__SAMPLE_STATS__INTERNAL_H__
