#ifndef __LINEAR_PREDICTION_MODEL_H__
#define __LINEAR_PREDICTION_MODEL_H__

#include <cassert>
#include <armadillo>

#include "fctrl.h"
#include "modelstructure.h"

namespace fctrl {

	class LinearPredictionModel;

	class LDS_State {
	public:
		ModelStructure * ms;

		arma::colvec sy;
			//Part of system state capturing the past *outputs* relevant for predicting next output y_t
			//sy = [y_(t-1), y_(t-2), ..., y_(t-ms->y_order)

		arma::colvec sx;
			//Part of system state capturing the past *inputs* relevant for predicting next output y_t
			//sx = [x_(t-1), x_(t-2), ..., x_(t-ms->x_order)

		arma::colvec offset_y;

		int affine_size;

		LDS_State()
			: ms(0), affine_size(0), saturated_num_transitions_applied(0)
		{}

		LDS_State(const LinearPredictionModel& lp, const arma::colvec& modelest_offset_y);

		bool isStateFullySetup() const {
			assert(ms != 0);
			return (saturated_num_transitions_applied >= ms->x_order and saturated_num_transitions_applied >= ms->y_order);
		}

		void state_transition(const arma::colvec& x_tm1, const arma::colvec& y_tm1);
			//input-output observed for completed frame

		arma::colvec get_full_LQR_state_repr() const;
			//state representation as a single vector, suitable for using with controller
			//  designed for corresponding LDS

	private:
		int saturated_num_transitions_applied;
	};



	class LinearPredictionModel {
	public:
		ModelStructure * ms;

		arma::mat L1;
		arma::mat L2;
		arma::mat L3;
		arma::colvec L4;
			//Define *affine* transform of past outputs, current input and past inputs for estimating y_t
			//y_t <- L1 * sy + L2 * xt + L3 * sx + L4,
			//   where xt is the input for the current time-step
		
		LinearPredictionModel(ModelStructure * ms = 0)
			: ms(ms)
		{ reset(); }

		void reset() {
			L1.reset();
			L2.reset();
			L3.reset();
			L4.reset();
		}

		void construct_linearpredictionmodel(const arma::colvec& q);

		arma::colvec get_prediction(const LDS_State& st, const arma::colvec& xt) const;
			//predicted current output yt under a given input xt (and past state)

	};


} //namespace fctrl

#endif //__LINEAR_PREDICTION_MODEL_H__
