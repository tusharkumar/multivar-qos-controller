#include "failuremodes.h"

namespace fctrl {

///////////////////////////////////////////////
/////////// class OscillationFailureMode
///////////////////////////////////////////////

double OscillationFailureMode::test_failure(double y) {
	assert(y_min <= y_max);

	switch(direction) {

	case 0: { //initially
		if(y < y_min) {
			halfcycle_min = y;
			halfcycle_max = y_min; //i.e., not seen yet
			direction = +1; //since starting below objective
			L = 1;
		}
		else if(y > y_max) {
			halfcycle_max = y;
			halfcycle_min = y_max; //i.e., not seen yet
			direction = -1; //since starting above objective
			L = 1;
		}
		//else, direction stays 0
		L = 0;
	}
	break;

	case +1: { //rising
		assert(halfcycle_min < y_min); //pre-condition
		L++;
		if(halfcycle_max > y_max) { //previously risen above objective
			if(y < halfcycle_max) { //dipped below peak, indicates start of new halfcycle
				eta = d * eta + (halfcycle_max - halfcycle_min) * W / (L-1);

				direction = -1;

				halfcycle_min = y_max;
				if(y < y_min) //already below objective
					halfcycle_min = y;

				L = 1;
				if(eta > get_threshold()) { //failure test
					double correction = get_threshold() / eta;
						//reduce error-scaling of y, so control-system reduces magnitude of variations in x
					eta = 0.0;
					return correction;
				}
			}
			else if(y > halfcycle_max) {
				halfcycle_max = y; //new max
			}
		}
		else { //previously not risen above objective
			if(y < y_min)
				halfcycle_min = y; //new min
			else if(y > y_max)
				halfcycle_max = y; //just risen above objective
		}
	}
	break;

	case -1: { //falling
		assert(halfcycle_max > y_max); //pre-condition
		L++;
		if(halfcycle_min < y_min) { //previously fallen below objective
			if(y > halfcycle_min) { //risen above trough, indicates start of new halfcycle
				eta = d * eta + (halfcycle_max - halfcycle_min) * W / (L-1);

				direction = +1;

				halfcycle_max = y_min;
				if(y > y_max) //already above objective
					halfcycle_max = y;

				L = 1;
				if(eta > get_threshold()) { //failure test
					double correction = get_threshold() / eta;
						//reduce error-scaling of y, so control-system reduces magnitude of variations in x
					eta = 0.0;
					return correction;
				}
			}
			else if(y < halfcycle_min) {
				halfcycle_min = y; //new min
			}
		}
		else { //previously not fallen below objective
			if(y > y_max)
				halfcycle_max = y; //new max
			else if(y < y_min)
				halfcycle_min = y; //just fallen below objective
		}
	}
	break;

	default: //invalid direction
	{ assert(0); }
	break;

	} //end switch(direction)

	return -1; //no failure detected
}



///////////////////////////////////////////////
/////////// class SluggishnessFailureMode
///////////////////////////////////////////////

double SluggishnessFailureMode::test_failure(double y) {
	assert(y_min <= y_max);

	switch(side) {
	case 0: {
		if(y > y_max) {
			side = +1;
			K = 1;
		}
		else if(y < y_min) {
			side = -1;
			K = 1;
		}
		//continue inside objective
		K = 0;
	}
	break;

	case +1: {
		if(y > y_max) {
			K++;

			double intermediate_mu = d * mu + ((double)(K-1)) / W;
			if(intermediate_mu > get_threshold()) { //long-running, interrupt failure sequence
				double correction = intermediate_mu / get_threshold();
					//increase error-scaling of y, so control-system increases magnitude of variations in x
				mu = 0.0;
				K = 1;
				return correction;
			}
		}
		else { //end of failure sequence
			mu = d * mu + ((double)(K-1)) / W;

			if(y < y_min) {
				side = -1;
				K = 1;
			}
			else { //within objective
				side = 0;
				K = 0;
			}

			if(mu > get_threshold()) { //failure test
				double correction = mu / get_threshold();
					//increase error-scaling of y, so control-system increases magnitude of variations in x
				mu = 0.0;
				return correction;
			}
		}
	}
	break;

	case -1: {
		if(y < y_min) {
			K++;

			double intermediate_mu = d * mu + ((double)(K-1)) / W;
			if(intermediate_mu > get_threshold()) { //long-running, interrupt failure sequence
				double correction = intermediate_mu / get_threshold();
					//increase error-scaling of y, so control-system increases magnitude of variations in x
				mu = 0.0;
				K = 1;
				return correction;
			}
		}
		else { //end of failure sequence
			mu = d * mu + ((double)(K-1)) / W;

			if(y > y_max) {
				side = +1;
				K = 1;
			}
			else { //within objective
				side = 0;
				K = 0;
			}

			if(mu > get_threshold()) { //failure test
				double correction = mu / get_threshold();
					//increase error-scaling of y, so control-system increases magnitude of variations in x
				mu = 0.0;
				return correction;
			}
		}
	}
	break;

	default: //invalid side
	{ assert(0); }
	break;

	} //end switch(side)

	return -1; //no failure detected
}

} //namespace fctrl
