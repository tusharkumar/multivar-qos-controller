#ifndef __FAILURE_MODE_H__
#define __FAILURE_MODE_H__

#include <cassert>

namespace fctrl {

	class OscillationFailureMode {
	public:
		int L;
			//length of current half-cycle

		double halfcycle_max;
			//maximum value reached by current halfcycle
		

		double halfcycle_min;
			//minimum value reached by current halfcycle

		int direction;
			//+1 ==> rising halfcycle  (halfcycle_min final)
			// 0 ==> not in a halfcycle (only initially)
			//-1 ==> falling halfcycle (halfcycle_max final)

		double eta;
			//converging sum

		double d;
			//0 < d < 1: determines "forget rate" in converging sum

		double y_max;
		double y_min;
			//[y_max, y_min] defines objective-range

		int W;
			//the perception window

		OscillationFailureMode(double y_max = 0.0, double y_min = 0.0, int W = 1, double d = 0.6)
			: L(0), halfcycle_max(0.0), halfcycle_min(0.0), direction(0), eta(0.0),
				d(d), y_max(y_max), y_min(y_min), W(W)
		{ }

		double get_threshold() const {
			double threshold = 1.0 / (1 - d) * (y_max - y_min) * 1.0;
			return threshold;
		}

		double test_failure(double y);
		//Invoked at each frame, determines if Oscillation failure is detected at current frame given y.
		//[y_min, y_max] is objective range for y.
		//
		//Return multiplicative correction to the error-scaling of y, if Oscillation failure detected.
		//Return -1 if no Oscillation failure detected

		void reset() {
		//reset metrics to start state
			L = 0;
			halfcycle_max = 0.0;
			halfcycle_min = 0.0;
			direction = 0;
			eta = 0.0;
		}
	};

	class SluggishnessFailureMode {
	public:

		long int K;

		int side;
			//+1 ==> stayed above objective continuously for K frames
			//0  ==> currently within objective
			//-1 ==> stayed below objective continuously for K frames

		double mu;
			//converging sum

		double d;
			//0 < d < 1: determines "forget rate" in converging sum

		double y_max;
		double y_min;
			//[y_max, y_min] defines objective-range

		int W;
			//the perception window

		SluggishnessFailureMode(double y_max = 0.0, double y_min = 0.0, int W = 1, double d = 0.6)
			: K(0), side(0), mu(0.0), d(d), y_max(y_max), y_min(y_min), W(W)
		{ }

		double get_threshold() const {
			double threshold = 1.0 / (1 - d) * 1;
			return threshold;
		}

		double test_failure(double y);
		//Invoked at each frame, determines if Sluggishness failure is detected at current frame given y.
		//[y_min, y_max] is objective range for y.
		//
		//Return multiplicative correction to the error-scaling of y, if Sluggishness failure detected.
		//Return -1 if no Sluggishness failure detected

		void reset() {
		//reset metrics to start state
			K = 0;
			side = 0;
			mu = 0.0;
		}
	};
} //namespace fctrl

#endif //__FAILURE_MODE_H__
