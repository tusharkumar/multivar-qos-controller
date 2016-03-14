#include <cassert>

#include "linearpredictionmodel.h"

namespace fctrl {

///////////////////////////////
//class LDS_State
///////////////////////////////


LDS_State::LDS_State(const LinearPredictionModel& lp, const arma::colvec& modelest_offset_y)
	: ms(lp.ms)
{
	assert(ms != 0);
	sx.reset();
	sy.reset();
	saturated_num_transitions_applied = 0;

	if(ms->x_order >= 1)
		sx.set_size(lp.L3.n_cols);

	if(ms->y_order >= 1)
		sy.set_size(lp.L1.n_cols);

	if(ms->flag_subtract_yoperating_in_estimation_for_DC_offset_correction)
		offset_y = modelest_offset_y / ms->get_OB_delta();
	else if(ms->flag_subtract_yobjective_in_estimation_for_DC_offset_correction)
		offset_y = ms->get_OB_center() / ms->get_OB_delta();
	else
		offset_y = arma::zeros<arma::colvec>((int)ms->vActive_OB_IDs.size());

	affine_size = (ms->flag_use_affine_model ? 1 : 0);
}

void LDS_State::state_transition(const arma::colvec& x_tm1, const arma::colvec& y_tm1)
{
	assert(ms != 0);
	int numXs = ms->vActive_CP_IDs.size();
	int numYs = ms->vActive_OB_IDs.size();

	assert((int)x_tm1.n_rows == numXs or ms->x_order < 0);
	assert((int)y_tm1.n_rows == numYs or ms->y_order < 1);

	if(ms->y_order >= 2) { //space to shift down
		sy.subvec(numYs, sy.n_rows-1) = sy.subvec(0, sy.n_rows - numYs - 1);
		sy.subvec(0, numYs-1) = y_tm1 - offset_y;
	}
	else if(ms->y_order == 1) {
		sy = y_tm1 - offset_y;
	}

	if(ms->x_order >= 2) { //space to shift down
		sx.subvec(numXs, sx.n_rows-1) = sx.subvec(0, sx.n_rows - numXs - 1);
		sx.subvec(0, numXs-1) = x_tm1;
	}
	else if(ms->x_order == 1) {
		sx = x_tm1;
	}

	if(ms->x_order == 0 and ms->y_order == 0) { //zero-order special case
		sy = y_tm1 - offset_y;
		sx = x_tm1;
	}

	saturated_num_transitions_applied++;
	if(saturated_num_transitions_applied > ms->x_order and saturated_num_transitions_applied > ms->y_order)
		saturated_num_transitions_applied--;
}


arma::colvec LDS_State::get_full_LQR_state_repr() const
{
	arma::colvec full_LQR_state(sy.n_rows + sx.n_rows + affine_size);
	if(sy.n_rows > 0)
		full_LQR_state.subvec(0, sy.n_rows-1) = sy;
	if(sx.n_rows > 0)
		full_LQR_state.subvec(sy.n_rows, full_LQR_state.n_rows - affine_size - 1) = sx;

	if(affine_size == 1)
		full_LQR_state(full_LQR_state.n_rows-1) = 1;
	else
		assert(affine_size == 0);
	
	return full_LQR_state;
}

///////////////////////////////
//class LinearPredictionModel
///////////////////////////////


void LinearPredictionModel::construct_linearpredictionmodel(const arma::colvec& q)
{
	assert(ms != 0);
	//sy.reset();
	//sx.reset();
	L1.reset();
	L2.reset();
	L3.reset();
	L4.reset();
	//saturated_num_transitions_applied = 0;

	int numXs = ms->vActive_CP_IDs.size();
	int numYs = ms->vActive_OB_IDs.size();
	assert(numYs > 0);

	int x_order_size = ms->x_order+1;
	int y_order_size = ms->y_order;
	int affine_size  = (ms->flag_use_affine_model ? 1 : 0);

	int x_block = numXs * x_order_size;
	int y_block = numYs * y_order_size;
	int block_size = x_block + y_block + affine_size;

	int num_coeffs = (numXs * x_order_size + numYs * y_order_size + affine_size) * numYs;
	
	assert((int)q.n_rows == num_coeffs);

	if(numXs > 0) {
		if(ms->x_order >= 0) { //xt used, need L2
			L2.set_size(numYs, numXs);
			for(int i=0; i<numYs; i++) {
				for(int k=0; k<numXs; k++) {
					int j = k;
					int r = 0;

					int index_into_q = i*block_size + j*x_order_size + r;
					L2(i,k) = q(index_into_q);
				}
			}
		}

		if(ms->x_order >= 1) { //sx used, need L3
			L3.set_size(numYs, numXs * ms->x_order);
			for(int i=0; i<numYs; i++) {
				for(int k=0; k<(int)L3.n_cols; k++) {
					int j = k % numXs;
					int r = k / numXs + 1;

					int index_into_q = i*block_size + j*x_order_size + r;
					L3(i,k) = q(index_into_q);
				}
			}
			//sx.set_size(L3.n_cols);
		}
	}

	if(ms->y_order >= 1) { //sy used, need L1
		L1.set_size(numYs, numYs * ms->y_order);
		for(int i=0; i<numYs; i++) {
			for(int k=0; k<(int)L1.n_cols; k++) {
				int j = k % numYs;
				int r = k / numYs + 1;

				int index_into_q = i*block_size + x_block + j*y_order_size + r-1;
				L1(i,k) = q(index_into_q);
			}
		}
		//sy.set_size(L1.n_cols);
	}

	if(affine_size == 1) { //affine model used, need L4
		L4.set_size(numYs);
		for(int i=0; i<numYs; i++) {
			int index_into_q = i*block_size + x_block + y_block;
			L4(i) = q(index_into_q);
		}
	}
	else
	{ assert(affine_size == 0); }
}


arma::colvec LinearPredictionModel::get_prediction(const LDS_State& st, const arma::colvec& xt) const
{
	assert(ms != 0);
	assert(xt.n_rows == ms->vActive_CP_IDs.size() or ms->x_order < 0);
	assert(st.isStateFullySetup() == true);

	arma::colvec yt = arma::zeros<arma::colvec>(ms->vActive_OB_IDs.size());
	if(!L1.is_empty())
		yt += L1 * st.sy;
	if(!L2.is_empty())
		yt += L2 * xt;
	if(!L3.is_empty())
		yt += L3 * st.sx;
	if(!L4.is_empty())
		yt += L4;

	return yt;
}


} //namespace fctrl
