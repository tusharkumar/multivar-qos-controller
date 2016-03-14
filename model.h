#ifndef __MODEL_H__
#define __MODEL_H__

#include <armadillo>
#include <iostream>
#include <deque>


#include "fctrl.h"
#include "frameobservation.h"
#include "modelstructure.h"
#include "linearpredictionmodel.h"


namespace fctrl {
	class Model {
	public:
		ModelStructure * ms;	

		double lambda;
		arma::colvec q;
			//The learned linear model coefficients

		//Following defined if q is not empty: set by llse_estimate()
		double objective_fit_sqerror;
		double q_magnitude_sqerror;

		LinearPredictionModel lp;
		double mfe;
			//model-fit-error: the model-tracking-error of the model using the
			//  history used for llse_estimate

		arma::colvec model_Y_operating; //WARNING: raw, not normalized!!

			//the following reduce computation for compute_mte()
			//through use of optimized_compute_mte()
		double cached_mte;
		LDS_State cached_hist_state;
		long long int last_frame_number_of_compute_mte;
		double cached_gamma;



		Model(ModelStructure * ms)
			: ms(ms), lambda(0.000001 * 0.000001), lp(ms)
		{ reset(); }

		void reset() {
			q.reset();
			objective_fit_sqerror = -1;
			q_magnitude_sqerror = -1;
			lp.reset();
			mfe = -1;
			model_Y_operating.reset();

			cached_mte = -1;
			last_frame_number_of_compute_mte = -1;
			cached_gamma = -1;
		}

		bool isDefined() const
		{ return !q.is_empty(); }

		void llse_estimate(
			const std::deque<FrameObservation>& deqFrameObservation,
			double lambda,
			double gamma
		);
		//Applies the given lambda and uses the given frame-history to
		//  - estimate a linear fit: q
		//  - compute fit metrics: objective_fit_sqerror, q_magnitude_sqerror
		//  - construct a linear prediction model: lp
		//  - fit error on given history data and provided lambda: mfe
		//isDefined() == true once llse_estimate() is called

		bool isBalanced() const;
			//checks whether opposing trade-off metrics are sufficiently
			//  balanced in magnitude to suggest that LLSE has taken all
			//  constraints into account during optimization

		double refine_lambda() const;
			//suggests an improved value for lambda based on current metrics

		std::pair<double, LDS_State> compute_mte(
			const std::deque<FrameObservation>& deqFrameObservation,
			double gamma
		);
		//Compute the Model-Tracking-Error of the current model w.r.t. to a given
		//  history (deqFrameObservation), and a history forget-rate parameter (gamma)
		//Error if isDefined() == false
		//Return the mte and the last state (if any)

		double optimized_compute_mte(
			const std::deque<FrameObservation>& deqFrameObservation,
			double gamma
		);
			//minimizes recomputation by using cached information from previous
			//  mte computation.
			//Note: the mte computed by optimized_compute_mte() can be more
			//   accurate than compute_mte() as maintaining cached_hist_state
			//   essentially allows mte computation over an infinite history
			//   rather than the necessarily finite history for compute_mte()

		double compute_mte_of_history_subrange(
			const std::deque<FrameObservation>& deqFrameObservation,
			long long int offset_from_latest_frame, // >= 0
			long long int num_frames, // not counting frames to setup state
			double gamma
		);
		//Returns 0 if offset_from_latest_frame, num_frames and state setup frames
		// exceed available history in deqFrameObservation.
		// compute_mte() used on suitably extracted subrange of deqFrameObservation
		//  to compute mte.
	};

} //namespace fctrl

std::ostream& operator<<(std::ostream& os, const fctrl::Model& m);


#endif //__MODEL_H__
