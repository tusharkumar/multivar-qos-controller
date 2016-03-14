#ifndef __LQR_H__
#define __LQR_H__

#include <list>

#include <armadillo>

#include "linearpredictionmodel.h"

//TODO: Operations to support
//  - Extend solution for one additional time-step, with possibly new R
//  - Simulate LQR+Linearmodel for N time-steps to generate all control-inputs
//       + incrementally, if LQR model is getting intermittently re-evaluated

namespace fctrl {

	double matrix_convergence(const arma::mat& matrix1, const arma::mat& matrix2);

	class LQR_Controller {
	public:
		//Dynamical system model
		arma::mat A;

		arma::mat B;

		arma::mat C;


		//Cost matrices
		arma::mat Q;
		arma::mat R;

		int N_steps;
			//horizon used for the LQR solution below


		//LQR solution: set by solve_lqr(), modified by increment_lqr_solution()
		std::list<arma::mat> list_K_t;
		std::list<arma::mat> list_W_t;
		std::list<arma::colvec> list_v_t;
		std::list<arma::mat> list_Kv_t;
			//ordered for time-steps 0 to (N_steps-1)

		bool isConverged;
			//indicates whether the LQR solution has converged


		//Intermediate solution values: these are maintained for increment_lqr_solution()
		std::list<arma::colvec> list_r_t;
			//list_r_t provides trajectory r_t for (N_steps+1) time-steps: 0 to N_steps
		arma::mat trans_A;
		arma::mat trans_B;
		arma::mat trans_C;

		arma::mat trans_C_x_Q;
		arma::mat trans_C_x_Q_x_C;


		////////////////

		LQR_Controller()
		{ reset(); }

		bool isDefined() const
		{ return (N_steps > 0); }

		void reset();

		void design_controller(LinearPredictionModel& lp, const int N_steps, const arma::colvec& diagR);
			//Sets up the dynamical system model, cost matrices and horizon, then invokes solve_lqr()


		void increment_controller_solution(const arma::colvec& diagR);
			//Increments N_steps, and extends previous LQR solution by one time-step,
			//  using a potentially modified R

		void solve_lqr();
			// Dynamical system model, cost matrices and horizon must be set up prior to invocation.

		arma::colvec get_control_input(const arma::colvec s_t0, const int t0) const;
			//Time-step t0 must be between 0 and (N_steps-1), both inclusive.
			//Here N_steps refers to the N_steps specified in a prior invocation of solve_lqr().
			//s_t0 provides the system-state at time-step t0.

		bool projected_inputs_exceed_bounds() const;
	};
} //namespace fctrl

#endif //__LQR_H__
