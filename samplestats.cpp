#include <iostream>
#include <deque>

#include "debuglog.h"
#include "samplestats.h"

namespace fctrl {

void add_sample_std(
		const std::vector<double>& v_x, // sample to add
		std::vector<double>& mean,      // input: mean, std before adding sample
		std::vector<double>& std,       // output: mean, std after adding sample
		std::size_t len,                // length of sample sequence prior to adding v_x
		double f                        // forget rate
)
{
	assert(mean.size() == v_x.size());
	assert(std.size() == v_x.size());

	if(len == 0) {
		for(std::size_t i=0; i<mean.size(); i++)
			mean[i] = 0.0;
		for(std::size_t i=0; i<std.size(); i++)
			std[i] = 0.0;
	}

	std::vector<double> var(std.size());
	for(std::size_t i=0; i<std.size(); i++)
		var[i] = std[i] * std[i];

	double f_to_len        = std::pow(f, len);
	double f_to_len_plus_1 = f_to_len * f;

	std::vector<double> new_mean(mean.size());
	for(std::size_t i=0; i<mean.size(); i++)
		new_mean[i] = ((1 - f) * v_x[i] + f * (1 - f_to_len) * mean[i]) /
		              (1 - f_to_len_plus_1);

	std::vector<double> new_var(var.size());
	for(std::size_t i=0; i<var.size(); i++)
		new_var[i] = (1 - f) / (1 - f_to_len_plus_1) *
		             ( f * (1 - f_to_len) / (1 - f) * var[i] + (v_x[i] - mean[i]) * (v_x[i] - new_mean[i]) );

	for(std::size_t i=0; i<mean.size(); i++)
		mean[i] = new_mean[i];

	for(std::size_t i=0; i<std.size(); i++)
		std[i] = std::sqrt(new_var[i] >= 0.0 ? new_var[i] : 0.0);
}

void remove_sample_std(
		const std::vector<double>& v_x, // sample to remove
		std::vector<double>& mean,      // input: mean, std before removing weighted sample
		std::vector<double>& std,       // output: mean, std after removing weighted sample
		std::size_t len,                // length of sample sequence prior to removing weighted sample
		double f                        // forget rate
)
{
	assert(mean.size() == v_x.size());
	assert(std.size() == v_x.size());

	if(len > 1) {
		double f_to_len_minus_1 = std::pow(f, len-1);
		double f_to_len         = f_to_len_minus_1 * f;

		std::vector<double> new_mean(mean.size());
		for(std::size_t i=0; i<mean.size(); i++)
			new_mean[i] = ((1 - f_to_len) * mean[i] - (1 - f) * f_to_len_minus_1 * v_x[i]) /
			              (1 - f_to_len_minus_1);

		std::vector<double> var(std.size());
		for(std::size_t i=0; i<std.size(); i++)
			var[i] = std[i] * std[i];

		std::vector<double> new_var(var.size());
		for(std::size_t i=0; i<var.size(); i++)
			new_var[i] = (1 - f) / (1 - f_to_len_minus_1) *
			             ((1 - f_to_len) / (1 - f) * (var[i] + (mean[i] - new_mean[i])*(mean[i] - new_mean[i])) -
			              f_to_len_minus_1 * (v_x[i] - new_mean[i]) * (v_x[i] - new_mean[i]));

		for(std::size_t i=0; i<mean.size(); i++)
			mean[i] = new_mean[i];

		for(std::size_t i=0; i<std.size(); i++)
			std[i] = std::sqrt(new_var[i] >= 0.0 ? new_var[i] : 0.0);
	}
	else {
		for(std::size_t i=0; i<mean.size(); i++)
			mean[i] = 0.0;
		for(std::size_t i=0; i<std.size(); i++)
			std[i] = 0.0;
	}
}

/////////////////////////////////////////////////////

void add_sample_maxswing(
		const std::vector<double>& v_x, // sample to add
		std::vector<double>& maxswing,  // input: maxswing, wmin and wmax before adding sample
		std::vector<double>& wmin,      // output: maxswing, wmin and wmax after adding sample
		std::vector<double>& wmax,
		std::size_t len,                // length of sample sequence prior to adding v_x
		double f                        // forget rate
)
{
	assert(maxswing.size() == v_x.size());
	assert(wmin.size() == v_x.size());
	assert(wmax.size() == v_x.size());

	if(len == 0) {
		for(std::size_t i=0; i<maxswing.size(); i++) {
			maxswing[i] = 0.0;
			wmin[i] = v_x[i];
			wmax[i] = v_x[i];
		}
		return;
	}

	std::vector<double> gwmin = wmin;
	std::vector<double> gwmax = wmax;

	for(size_t i=0; i<gwmin.size(); i++) {
		gwmin[i] *= f;
		gwmax[i] *= f;
	}

	for(size_t i=0; i<maxswing.size(); i++) {
		maxswing[i] = std::max( f * maxswing[i], std::max(std::abs(v_x[i] - gwmin[i]), std::abs(v_x[i] - gwmax[i])) );
		wmin[i] = std::min(v_x[i], gwmin[i]);
		wmax[i] = std::max(v_x[i], gwmax[i]);
	}
}



/////////////////////////////////////////////////////

void add_sample(
		const std::vector<double>& v_x, // sample to add
		std::vector<double>& mean,      // input: mean, std, maxswing, wmin and wmax before adding sample
		std::vector<double>& std,       // output: mean, std, maxswing, wmin and wmax after adding sample
		std::vector<double>& maxswing,
		std::vector<double>& wmin,
		std::vector<double>& wmax,
		std::size_t len,                // length of sample sequence prior to adding v_x
		double f                        // forget rate
)
{
	add_sample_std(v_x, mean, std, len, f);
	add_sample_maxswing(v_x, maxswing, wmin, wmax, len, f);
}

} //namespace fctrl

//#define UNITTESTS_SAMPLESTATS
#ifdef UNITTESTS_SAMPLESTATS


int main() {

	// 1-d test
	std::vector<double> mean, std;
	mean.push_back(0.0);
	std.push_back(0.0);

	std::vector<double> maxswing, wmin, wmax;
	maxswing.push_back(0);
	wmin.push_back(0);
	wmax.push_back(0);

	const double f = 0.9;

	std::deque< std::vector<double> > h; //h[0] is newest, h[h.size()-1] is oldest sample
	for(int i=1; i<10; i++) {
		std::vector<double> v_x;
		v_x.push_back(i);
		fctrl::add_sample(v_x, mean, std, maxswing, wmin, wmax, h.size(), f);
#ifdef FC_SAMPLESTATS_DEBUG
		std::cout << "added " << i << ": mean[0]=" << mean[0] << " std[0]=" << std[0] << std::endl;
		std::cout << "         maxswing[0]=" << maxswing[0] << " wmin[0]=" << wmin[0] << " wmax[0]=" << wmax[0] << std::endl;
#endif //FC_SAMPLESTATS_DEBUG
		h.push_front(v_x);
	}

#if 1 //in-order removal
	for(int i=1; i<10; i++) {
		std::vector<double> v_x;
		v_x.push_back(i);
		fctrl::remove_sample(h, mean, std, maxswing, wmin, wmax, 10-i, f);
#ifdef FC_SAMPLESTATS_DEBUG
		std::cout << "removed " << i << ": mean[0]=" << mean[0] << " std[0]=" << std[0] << std::endl;
		std::cout << "          maxswing[0]=" << maxswing[0] << " wmin[0]=" << wmin[0] << " wmax[0]=" << wmax[0] << std::endl;
#endif //FC_SAMPLESTATS_DEBUG
	}

#else //reverse-order removal
	for(int i=9; i>0; i--) {
		std::vector<double> v_x;
		v_x.push_back(i);
		fctrl::remove_sample_std(v_x, mean, std,i,f);
#ifdef FC_SAMPLESTATS_DEBUG
		std::cout << "removed " << i << ": mean[0]=" << mean[0] << " std[0]=" << std[0] << std::endl;
#endif //FC_SAMPLESTATS_DEBUG
	}
#endif
	
	return 0;
}

#endif //UNITTESTS_SAMPLESTATS

