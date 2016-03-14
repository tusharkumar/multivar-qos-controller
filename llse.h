#ifndef __LLSE_H__
#define __LLSE_H__

#include <armadillo>
#include <iostream>
#include <deque>

#include "fctrl.h"
#include "frameobservation.h"
#include "modelstructure.h"


namespace fctrl {

	bool solve_LLSE(
		ModelStructure * ms,
		const std::deque<FrameObservation>& deqFrameObservation,
		const arma::colvec X_normalizer,
		const arma::colvec Y_normalizer,
		const arma::colvec Y_importances,
		const double lambda,
		arma::colvec& q, //results
		double& objective_fit_sqerror,
		double& q_magnitude_sqerror
	);
	//Uses given sequence of FrameObservations to solve for q, as directed by ms.
	//Before LLSE is applied, the ControlParameter observations are normalized by X_normalizer,
	//  and the Objective observations are normalized by Y_normalizer.
	// Also, the relative importance of fitting the i'th Objective is given by Y_importances(i).
	//
	//Returns true if model-estimation was a success, false otherwise.
}

#endif //__LLSE_H__
