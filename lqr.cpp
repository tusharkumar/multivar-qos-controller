#include <cassert>

#include "debuglog.h"
#include "utils.h"
#include "lqr.h"

namespace fctrl {

double matrix_convergence(const arma::mat& matrix1, const arma::mat& matrix2)
{
	arma::mat abs_diff = arma::abs(matrix1 - matrix2);
	double sum_abs_diff = arma::sum(arma::sum(abs_diff));

	double abs_sum1 = arma::sum(arma::sum(arma::abs(matrix1)));
	double abs_sum2 = arma::sum(arma::sum(arma::abs(matrix2)));

	double avg_abs_sum = (abs_sum1 + abs_sum2)/2;

	double convergence_error = 0.0;
	if(sum_abs_diff > 0)
		convergence_error = sum_abs_diff / avg_abs_sum;

	//std::cout << "matrix1 = " << matrix1 << std::endl;
	//std::cout << "matrix2 = " << matrix2 << std::endl;
	//std::cout << "convergence_error=" << convergence_error << std::endl;

	return convergence_error;
}


void LQR_Controller::reset()
{
	A.reset();
	B.reset();
	C.reset();

	N_steps = 0;

	list_K_t.clear();
	list_W_t.clear();
	list_v_t.clear();
	list_Kv_t.clear();

	isConverged = false;

	trans_A.reset();
	trans_B.reset();
	trans_C.reset();

	trans_C_x_Q.reset();
	trans_C_x_Q_x_C.reset();
}


void LQR_Controller::design_controller(LinearPredictionModel& lp, const int N_steps, const arma::colvec& diagR)
{
	assert(N_steps > 0);
	this->N_steps = N_steps;

	const arma::mat& L1 = lp.L1;
	const arma::mat& L2 = lp.L2;
	const arma::mat& L3 = lp.L3;
	const arma::colvec& L4 = lp.L4;

	int numXs = lp.ms->vActive_CP_IDs.size();
	int numYs = lp.ms->vActive_OB_IDs.size();

	int x_order = lp.ms->x_order;
	int y_order = lp.ms->y_order;
	
	int past_inputs_size  = numXs * x_order;
	int past_outputs_size = numYs * y_order;
	int affine_size = (lp.ms->flag_use_affine_model ? 1 : 0);

	int state_size = past_inputs_size + past_outputs_size + affine_size;


	assert(state_size > 0);
	A.set_size(state_size, state_size);

	assert((int)L1.n_cols == past_outputs_size); //both could be 0
	if(!L1.is_empty()) {
		assert((int)L1.n_rows == numYs);
		A( arma::span(0, numYs-1), arma::span(0, past_outputs_size-1) ) = L1;
	}

	assert((int)L3.n_cols == past_inputs_size); //both could be 0
	if(!L3.is_empty()) {
		assert((int)L3.n_rows == numYs);
		A( arma::span(0, numYs-1), arma::span(past_outputs_size, state_size - affine_size - 1) ) = L3;
	}

	if(!L4.is_empty()) {
		assert(affine_size == 1);
		assert((int)L4.n_rows == numYs);
		A( arma::span(0, numYs-1), A.n_cols-1 ) = L4;
	}

	if(state_size > numYs) {
		arma::mat Bs = arma::zeros<arma::mat>(state_size - numYs, state_size);
		if(past_outputs_size - numYs > 0)
			Bs( arma::span(0, past_outputs_size - numYs - 1), arma::span(0, past_outputs_size - numYs - 1) )
								= arma::eye<arma::mat>(past_outputs_size - numYs, past_outputs_size - numYs);
		if(past_inputs_size > numXs) //i.e., at least two previous time-steps involved in state
			Bs(
				arma::span(Bs.n_rows - affine_size - past_inputs_size + numXs, Bs.n_rows - affine_size - 1),
				arma::span(past_outputs_size, Bs.n_cols - affine_size - numXs - 1) )
								= arma::eye<arma::mat>(past_inputs_size - numXs, past_inputs_size - numXs);
		if(affine_size == 1)
			Bs(Bs.n_rows-1, Bs.n_cols-1) = 1; //to generate 1 at the bottom of the next state vector
		else
			assert(affine_size == 0);

		A( arma::span(numYs, A.n_rows-1), arma::span::all ) = Bs;
	}


	B.zeros(state_size, numXs); //use .zeros() explicitly instead of .set_size() as the affine rows are not overwritten
	assert(lp.ms->x_order < 0 or (int)L2.n_cols == numXs);
	if(!L2.is_empty()) {
		assert((int)L2.n_rows == numYs);
		B( arma::span(0, L2.n_rows-1), arma::span::all) = L2;
	}
	if(past_inputs_size > 0) { //i.e. next state needs current input, so S should be used
		arma::mat S = arma::zeros<arma::mat>(state_size - numYs, numXs);
		S( arma::span(past_outputs_size - numYs, past_outputs_size - numYs + numXs - 1), arma::span::all ) = arma::eye<arma::mat>(numXs, numXs);

		B( arma::span(L2.n_rows, B.n_rows-1), arma::span::all ) = S;
	}

	C = arma::zeros<arma::mat>(numYs, state_size);
	C( arma::span::all, arma::span(0, numYs-1) ) = arma::eye<arma::mat>(numYs, numYs);


	Q = arma::diagmat( lp.ms->get_OB_importance() );
		//dealing with a normalized system model, so deltas = 1.0
		//ASSUMPTION: ModelStructure ms is unchanged since LinearEstimator was run and LinearPredictionModel was constructed
		//   (so that the normalizing factors, arising from Y_deltas, are consistent)

	assert((int)diagR.n_rows == numXs);
	R = arma::diagmat(diagR);
#ifdef FC_LQR_DEBUG
	std::cout << "R: " << R << std::endl;
#endif //FC_LQR_DEBUG

	arma::colvec r_t = lp.ms->get_OB_center() / lp.ms->get_OB_delta();
	list_r_t.clear();
	for(int t=0; t<this->N_steps+1; t++) //r_t for t=0 to N_steps+1
		list_r_t.push_back(r_t);

	solve_lqr();

}




void LQR_Controller::solve_lqr()
{
	assert(N_steps > 0);

	assert(!Q.is_empty());
	assert(Q.n_rows == Q.n_cols);

	assert(!R.is_empty());
	assert(R.n_rows == R.n_cols);

	assert(!A.is_empty());
	assert(A.n_rows == A.n_cols);

	assert(!B.is_empty());
	assert(B.n_rows == A.n_rows);

	assert(!C.is_empty());
	assert(C.n_cols == A.n_rows);

	list_K_t.clear();
	list_W_t.clear();
	list_v_t.clear();
	list_Kv_t.clear();
	isConverged = false;

	//setup intermediate values
	trans_A = arma::trans(A);
	trans_B = arma::trans(B);
	trans_C = arma::trans(C);

	trans_C_x_Q = trans_C * Q;
	trans_C_x_Q_x_C = trans_C_x_Q * C;


	arma::mat W_tp1 = trans_C_x_Q_x_C;
	list_W_t.push_front(W_tp1);

	std::list< arma::colvec >::const_reverse_iterator ptr_r_t = list_r_t.rbegin();
	arma::colvec v_tp1 = trans_C * Q * (*ptr_r_t);
	list_v_t.push_front(v_tp1);

	//std::cout << "solve_lqr(): N_steps=" << N_steps << " list_r_t=" << list_r_t << std::endl;
	//std::cout << "A = " << A << std::endl;
	//std::cout << "B = " << B << std::endl;
	//std::cout << "C = " << C << std::endl;

	for(int t=N_steps-1; t>=0; t--) {
		ptr_r_t++;
		assert(ptr_r_t != list_r_t.rend());

		arma::mat inv_term_t = arma::inv(trans_B * W_tp1 * B + R);

		arma::mat K_t    = inv_term_t * trans_B * W_tp1 * A;

		arma::mat A_minus_B_x_K_t = (A - B * K_t);

		arma::mat W_t    = trans_A * W_tp1 * A_minus_B_x_K_t + trans_C_x_Q_x_C;
		arma::colvec v_t = arma::trans(A_minus_B_x_K_t) * v_tp1 + trans_C_x_Q * (*ptr_r_t);
		arma::mat Kv_t   = inv_term_t * trans_B;

		list_K_t.push_front(K_t);
		list_W_t.push_front(W_t);
		list_v_t.push_front(v_t);
		list_Kv_t.push_front(Kv_t);

		W_tp1 = W_t;
		v_tp1 = v_t;
	}
	assert((++ptr_r_t) == list_r_t.rend());

	if(N_steps > 1) {
		const arma::mat& K_t = *(list_K_t.begin());
		const arma::mat& K_tp1 = *((list_K_t.begin())++);

		double K_t_convergence_error = matrix_convergence(K_t, K_tp1);

		const arma::mat& Kv_t = *(list_Kv_t.begin());
		const arma::mat& Kv_tp1 = *((list_Kv_t.begin())++);

		double Kv_t_convergence_error = matrix_convergence(Kv_t, Kv_tp1);

		if(K_t_convergence_error < 0.01 and Kv_t_convergence_error < 0.01)
			isConverged = true;
	}

	//std::cout << "---------- LQR SOLUTION ------------" << std::endl;
	//std::cout << "list_W_t =" << list_W_t << std::endl;
	//std::cout << "list_K_t =" << list_K_t << std::endl;
	//std::cout << "list_Kv_t=" << list_Kv_t << std::endl;
	//std::cout << "list_v_t=" << list_v_t << std::endl;
}




void LQR_Controller::increment_controller_solution(const arma::colvec& diagR)
{
	assert(isDefined());
	N_steps++;

	assert(diagR.n_rows == R.n_rows);
	R = arma::diagmat(diagR);

	arma::colvec r_t = *(list_r_t.begin());
	list_r_t.push_front(r_t);

	////

	const arma::mat& W_tp1 = *(list_W_t.begin());
	const arma::colvec& v_tp1 = *(list_v_t.begin());

	//One more time-step
	{
		arma::mat inv_term_t = arma::inv(trans_B * W_tp1 * B + R);

		arma::mat K_t    = inv_term_t * trans_B * W_tp1 * A;

		arma::mat A_minus_B_x_K_t = (A - B * K_t);

		arma::mat W_t    = trans_A * W_tp1 * A_minus_B_x_K_t + trans_C_x_Q_x_C;
		arma::colvec v_t = arma::trans(A_minus_B_x_K_t) * v_tp1 + trans_C_x_Q * r_t;
		arma::mat Kv_t   = inv_term_t * trans_B;

		list_K_t.push_front(K_t);
		list_W_t.push_front(W_t);
		list_v_t.push_front(v_t);
		list_Kv_t.push_front(Kv_t);
	}

	assert(N_steps > 1);
	{
		const arma::mat& K_t = *(list_K_t.begin());
		const arma::mat& K_tp1 = *((list_K_t.begin())++);

		double K_t_convergence_error = matrix_convergence(K_t, K_tp1);

		const arma::mat& Kv_t = *(list_Kv_t.begin());
		const arma::mat& Kv_tp1 = *((list_Kv_t.begin())++);

		double Kv_t_convergence_error = matrix_convergence(Kv_t, Kv_tp1);

		if(K_t_convergence_error < 0.01 and Kv_t_convergence_error < 0.01)
			isConverged = true;
	}

	//std::cout << "---------- INCREMENTED LQR SOLUTION ------------" << std::endl;
	//std::cout << "list_W_t =" << list_W_t << std::endl;
	//std::cout << "list_K_t =" << list_K_t << std::endl;
	//std::cout << "list_Kv_t=" << list_Kv_t << std::endl;
	//std::cout << "list_v_t=" << list_v_t << std::endl;
}






arma::colvec LQR_Controller::get_control_input(const arma::colvec s_t0, const int t0) const
{
	assert(t0 >= 0);

	std::list<arma::mat>::const_iterator ptr_K_t      = list_K_t.begin();
	std::list<arma::mat>::const_iterator ptr_Kv_t     = list_Kv_t.begin();
	std::list<arma::colvec>::const_iterator ptr_v_tp1 = list_v_t.begin();

	assert(ptr_K_t   != list_K_t.end());
	assert(ptr_Kv_t  != list_Kv_t.end());
	assert(ptr_v_tp1 != list_v_t.end());

	ptr_v_tp1++;
	assert(ptr_v_tp1 != list_v_t.end());

	for(int t=0; t<t0; t++) {
		ptr_K_t++;
		ptr_Kv_t++;
		ptr_v_tp1++;

		assert(ptr_K_t   != list_K_t.end());
		assert(ptr_Kv_t  != list_Kv_t.end());
		assert(ptr_v_tp1 != list_v_t.end());
	}

#ifdef FC_LQR_DEBUG
	std::cout << "get_control_input(): s_t = " << s_t0 << std::endl;
	std::cout << "  K_t = " << *ptr_K_t << std::endl;
	std::cout << "  v_(t+1) = " << *ptr_v_tp1 << std::endl;
	std::cout << "  Kv_t = " << *ptr_Kv_t << std::endl;
#endif //FC_LQR_DEBUG

	arma::colvec u_t0 = - (*ptr_K_t) * s_t0 + (*ptr_Kv_t) * (*ptr_v_tp1);

	return u_t0;
}

} //namespace fctrl
