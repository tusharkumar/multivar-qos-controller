#include <cstdlib>
#include <cassert>

#include "debuglog.h"
#include "utils.h"
#include "llse.h"

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
)
{
	assert(ms != 0);
	assert(ms->x_order >= 0);
	assert(ms->y_order >= 0);
	assert(X_normalizer.n_rows == ms->vActive_CP_IDs.size());
	assert(Y_normalizer.n_rows == ms->vActive_OB_IDs.size());

	int numXs = ms->vActive_CP_IDs.size();
	int numYs = ms->vActive_OB_IDs.size();

	int x_order_size = ms->x_order+1;
	int y_order_size = ms->y_order;
	int affine_size = (ms->flag_use_affine_model ? 1 : 0);

	int num_rows = numYs;
	int num_cols = (numXs * x_order_size + numYs * y_order_size + affine_size) * numYs;

	for(int f=0; f<(int)deqFrameObservation.size(); f++) {
		const FrameObservation& fo = deqFrameObservation[f];

		assert(fo.vCP_observations.size() == ms->vActive_CP_IDs.size());
		assert(fo.vOB_observations.size() == ms->vActive_OB_IDs.size());
	}

	int x_history_length = deqFrameObservation.size() - ms->x_order;
	int y_history_length = deqFrameObservation.size() - ms->y_order;
	int usable_history_length = std::min(x_history_length, y_history_length);

#ifdef FC_LLSE_DEBUG
	std::cout << "LinearEstimator:: usable_history_length = " << usable_history_length << std::endl;
	std::cout << "LinearEstimator:: HISTORY:" << std::endl;
	std::cout << deqFrameObservation << std::endl;
#endif //FC_LLSE_DEBUG
	if(usable_history_length <= 0)
		return false; //FAILURE: insufficient history to estimate model, based on structuring parameters provided for order

	assert(
		   (ms->flag_subtract_ytm1_in_estimation_for_DC_offset_correction == true
		and ms->flag_subtract_yoperating_in_estimation_for_DC_offset_correction == false
		and ms->flag_subtract_yobjective_in_estimation_for_DC_offset_correction == false)

	  	or (ms->flag_subtract_ytm1_in_estimation_for_DC_offset_correction == false
		and ms->flag_subtract_yoperating_in_estimation_for_DC_offset_correction == true
		and ms->flag_subtract_yobjective_in_estimation_for_DC_offset_correction == false)

		or (ms->flag_subtract_ytm1_in_estimation_for_DC_offset_correction == false
		and ms->flag_subtract_yoperating_in_estimation_for_DC_offset_correction == false
		and ms->flag_subtract_yobjective_in_estimation_for_DC_offset_correction == true)

		or (ms->flag_subtract_ytm1_in_estimation_for_DC_offset_correction == false
		and ms->flag_subtract_yoperating_in_estimation_for_DC_offset_correction == false
		and ms->flag_subtract_yobjective_in_estimation_for_DC_offset_correction == false)
	);

	assert(ms->flag_subtract_ytm1_in_estimation_for_DC_offset_correction == false or ms->y_order >= 1);

	arma::mat A(num_rows * usable_history_length, num_cols);
	for(int f=0; f<usable_history_length; f++) {

		std::vector<double> offset_y(numYs, 0.0); //zero vector
		if(ms->flag_subtract_ytm1_in_estimation_for_DC_offset_correction == true) {
			int dep_ytm1 = f + 1;
			const FrameObservation& dep_ytm1_fo = deqFrameObservation.at(dep_ytm1);
			offset_y = dep_ytm1_fo.vOB_observations;
		}
		else if(ms->flag_subtract_yoperating_in_estimation_for_DC_offset_correction == true) {
			offset_y = get_std_for_arma_colvec( ms->Y_operating );
		}
		else if(ms->flag_subtract_yobjective_in_estimation_for_DC_offset_correction == true) {
			offset_y = get_std_for_arma_colvec( ms->get_OB_center() );
		}

		arma::mat P_f;
		P_f.zeros(num_rows, num_cols);
		for(int i=0; i<numYs; i++) { //iterate over y_i's

			arma::rowvec p_i_f_xpart(numXs * x_order_size);
			for(int j=0; j<numXs; j++) { //iterate over the x_j's
				for(int o = 0; o <= ms->x_order; o++) { //iterate over dependent time-steps
					int dep = f + o;

					const FrameObservation& dep_fo = deqFrameObservation.at(dep);
					p_i_f_xpart(j * x_order_size + o) = dep_fo.vCP_observations.at( j ) / X_normalizer(j);
								//observation corresponding to ms->vActive_CP_IDs[j] from the (f+o)'th frame
				}
			}
			p_i_f_xpart *= Y_importances(i);

			arma::rowvec p_i_f_ypart(numYs * y_order_size);
			for(int j=0; j<numYs; j++) { //iterate over the y_j's
				for(int o = 1; o <= ms->y_order; o++) { //iterate over dependent time-steps
					int dep = f + o;

					const FrameObservation& dep_fo = deqFrameObservation.at(dep);
					p_i_f_ypart(j * y_order_size + o-1) = (dep_fo.vOB_observations.at( j ) - offset_y.at( j )) / Y_normalizer(j);
								//observation corresponding to ms->vActive_OB_IDs[j] from the (f+o)'th frame
				}
			}
			p_i_f_ypart *= Y_importances(i);

			arma::rowvec p_i_f_cpart(affine_size);
			if(p_i_f_cpart.n_cols == 1)
				p_i_f_cpart(0) = 1;
			else
				assert(affine_size == 0);

			int p_i_f_size = p_i_f_xpart.n_cols + p_i_f_ypart.n_cols + p_i_f_cpart.n_cols;

			if(p_i_f_xpart.is_empty() == false)
				P_f(i, arma::span(
							i * p_i_f_size,
							i * p_i_f_size + p_i_f_xpart.n_cols - 1
					) ) = p_i_f_xpart;

			if(p_i_f_ypart.is_empty() == false)
				P_f(i, arma::span(
							i * p_i_f_size + p_i_f_xpart.n_cols,
							i * p_i_f_size + p_i_f_xpart.n_cols + p_i_f_ypart.n_cols - 1
					) ) = p_i_f_ypart;

			if(p_i_f_cpart.is_empty() == false)
				P_f(i, arma::span(
							i * p_i_f_size + p_i_f_xpart.n_cols + p_i_f_ypart.n_cols,
							i * p_i_f_size + p_i_f_xpart.n_cols + p_i_f_ypart.n_cols + p_i_f_cpart.n_cols - 1
					) ) = p_i_f_cpart;

		}

		A(arma::span(f * P_f.n_rows, (f+1) * P_f.n_rows - 1), arma::span::all) = P_f;
	}
	//Now: matrix A fully constructed


	arma::colvec vdes(A.n_rows);
	for(int f=0; f<usable_history_length; f++) {

		std::vector<double> offset_y(numYs, 0.0); //zero vector
		if(ms->flag_subtract_ytm1_in_estimation_for_DC_offset_correction == true) {
			int dep_ytm1 = f + 1;
			const FrameObservation& dep_ytm1_fo = deqFrameObservation.at(dep_ytm1);
			offset_y = dep_ytm1_fo.vOB_observations;
		}
		else if(ms->flag_subtract_yoperating_in_estimation_for_DC_offset_correction == true) {
			offset_y = get_std_for_arma_colvec( ms->Y_operating );
		}
		else if(ms->flag_subtract_yobjective_in_estimation_for_DC_offset_correction == true) {
			offset_y = get_std_for_arma_colvec( ms->get_OB_center() );
		}

		const FrameObservation& fo = deqFrameObservation[f];

		arma::colvec y_f(num_rows);
		for(int i=0; i<numYs; i++) { //iterate over y_i's
			y_f(i) = (fo.vOB_observations.at( i ) - offset_y.at( i )) / Y_normalizer(i) * Y_importances(i);
								//observation corresponding to ms->vActive_OB_IDs[i] from the f'th frame
		}

		vdes.subvec(f * y_f.n_rows, (f+1) * y_f.n_rows -1) = y_f;
	}
	//Now: matrix vdes fully constructed
	
#ifdef FC_LLSE_DEBUG
	std::cout << "----LLSE:" << std::endl;
	std::cout << "A=" << A << std::endl;
	std::cout << "vdes=" << vdes << std::endl;
#endif //FC_LLSE_DEBUG

	assert(lambda >= 0.0);
	double sqrt_lambda = sqrt(lambda);

	//Now construct and solve the regularized version of the problem
	arma::mat A_reg(A.n_rows + A.n_cols, A.n_cols);
	A_reg(arma::span(0, A.n_rows-1), arma::span(0, A.n_cols-1)) = A;
	A_reg(arma::span(A.n_rows, A_reg.n_rows-1), arma::span(0, A.n_cols-1))
		= sqrt_lambda * arma::eye<arma::mat>(A.n_cols, A.n_cols);

	arma::colvec v_reg = arma::zeros<arma::colvec>(A_reg.n_rows);
	v_reg.subvec(0, vdes.n_rows-1) = vdes;

	//throws std::runtime_error exception on failure to find solution
	q = arma::solve(A_reg, v_reg);
	if(q.n_rows == 0) //older versions of armadillo don't throw exception here
		throw std::runtime_error("arma::solve() failed, throwing on its behalf");
	assert((int)q.n_rows == num_cols);

#ifdef FC_LLSE_DEBUG
	std::cout << "LinearEstimator:: q = " << q << std::endl;
#endif //FC_LLSE_DEBUG


	//Fit metrics
	objective_fit_sqerror = squared( arma::norm(A * q - vdes, 2) );
	q_magnitude_sqerror   = squared( arma::norm(q, 2) );


	if(ms->flag_subtract_ytm1_in_estimation_for_DC_offset_correction == true) {
		int q_x_coeffs_size = numXs * x_order_size;
		int q_each_y_i_coeffs_size = (numXs * x_order_size) + (numYs * y_order_size);

		for(int i=0; i<numYs; i++) { //iterate over y_i's
			for(int j=0; j<numYs; j++) { //iterate over the y_j's
				double coeff_sum_ytm1_j = 0.0;
				for(int o = 1; o <= ms->y_order; o++) { //iterate over dependent time-steps
					int index = i * q_each_y_i_coeffs_size + (q_x_coeffs_size + j * y_order_size + o-1);
					coeff_sum_ytm1_j += q(index);
				}
				int replace_j_index = i * q_each_y_i_coeffs_size + (q_x_coeffs_size + j * y_order_size + 1-1);
				double replace_ytm1_j_coeff = q(replace_j_index) - coeff_sum_ytm1_j;
				if(i == j)
					replace_ytm1_j_coeff += 1.0;
				q(replace_j_index) = replace_ytm1_j_coeff;
			}
		}

#ifdef FC_LLSE_DEBUG
		std::cout << "LinearEstimator:: AFTER adjusting for ytm1 DC offset: q = " << q << std::endl;
#endif //FC_LLSE_DEBUG
	}
	//else if(ms->flag_subtract_yoperating_in_estimation_for_DC_offset_correction == true) {
	//		nothing needed here: just need to make sure that state updates using the estimated model
	//		subtract out the same Y_operating value from any subsequent y-readings as well   
	//}

	return true; //SUCCESS
}

} //namespace fctrl
