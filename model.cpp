#include <cassert>

#include "debuglog.h"
#include "utils.h"
#include "llse.h"
#include "model.h"

namespace fctrl {

void Model::llse_estimate(
	const std::deque<FrameObservation>& deqFrameObservation,
	double lambda,
	double gamma
)
{
	assert(lambda >= 0.0);
	this->lambda = lambda;

	assert(ms != 0);
	ms->estimate_operating_point_from_data_samples(deqFrameObservation, gamma);
	model_Y_operating = ms->Y_operating;
		//model_Y_operating: will be needed during state-updates if ms->flag_subtract_yoperating_in_estimation_for_DC_offset_correction == true
	
#ifdef FC_MODEL_DEBUG
	std::cout << "model_Y_operating=" << model_Y_operating << std::endl;
#endif //FC_MODEL_DEBUG

	arma::colvec X_normalizer = ms->get_CP_N();
	arma::colvec Y_normalizer = ms->get_OB_delta();
	arma::colvec Y_importances = ms->get_OB_importance();
		//not control-importance, since the linear-estimation needs to be true to specific relative weightings of the Objectives

	bool wasLLSE_successful
		= solve_LLSE(ms, deqFrameObservation, X_normalizer, Y_normalizer, Y_importances, lambda,
						q, objective_fit_sqerror, q_magnitude_sqerror);
	assert(wasLLSE_successful);

	lp.construct_linearpredictionmodel(q);

	mfe = optimized_compute_mte(deqFrameObservation, gamma);
#ifdef FC_MODEL_DEBUG
	std::cout << "mfe=" << mfe << std::endl;
#endif //FC_MODEL_DEBUG
}


bool Model::isBalanced() const
{
	assert(isDefined());
	return (objective_fit_sqerror > 1/10.0 * lambda * q_magnitude_sqerror)
		and (objective_fit_sqerror < 10.0 * lambda * q_magnitude_sqerror);
}


double Model::refine_lambda() const
{
	assert(isDefined());
	double lambda_next;
#ifdef FC_MODEL_DEBUG
	std::cout << "lambda=" << lambda << " objective_fit_sqerror=" << objective_fit_sqerror << " q_magnitude_sqerror=" << q_magnitude_sqerror << std::endl;
#endif //FC_MODEL_DEBUG
  assert(lambda > 0.0);
	if(q_magnitude_sqerror > 0.0) {
		//regular update
		lambda_next = sqrt(objective_fit_sqerror * lambda / q_magnitude_sqerror);
	}

	if(q_magnitude_sqerror == 0.0 || 1.0/lambda_next == 0.0) {
		//handle pathological condition
		if(lambda > 1)
			lambda_next = sqrt(lambda);
		else
			lambda_next = lambda * lambda;
	}

  if(lambda_next == 0) {
		//handle pathological condition
		lambda_next = 10E-6; //heuristic
  }

	return lambda_next;
}


std::pair<double, LDS_State> Model::compute_mte(
	const std::deque<FrameObservation>& deqFrameObservation,
	double gamma
)
{
	assert(ms != 0);
	assert(isDefined());
	for(int f=0; f<(int)deqFrameObservation.size(); f++) {
		const FrameObservation& fo = deqFrameObservation[f];

		assert(fo.vCP_observations.size() == ms->vActive_CP_IDs.size());
		assert(fo.vOB_observations.size() == ms->vActive_OB_IDs.size());
	}
	
	int x_history_length = deqFrameObservation.size() - ms->x_order;
	int y_history_length = deqFrameObservation.size() - ms->y_order;
	int usable_history_length = std::min(x_history_length, y_history_length);

	double mte = 0.0;
	if(usable_history_length <= 0)
		return std::make_pair(mte, LDS_State());

	LDS_State hist_state(lp, model_Y_operating);
	int num_startup_frames = std::max(ms->x_order, ms->y_order);
	for(int s=1; s<=num_startup_frames; s++) {
		assert(!hist_state.isStateFullySetup());
		const FrameObservation& fo = deqFrameObservation.at(deqFrameObservation.size() - s);

		arma::colvec x = get_arma_colvec_for_std( fo.vCP_observations ) / ms->get_CP_N();
		arma::colvec y = get_arma_colvec_for_std( fo.vOB_observations ) / ms->get_OB_delta();

		hist_state.state_transition(x, y);
	}
	assert(hist_state.isStateFullySetup());

	for(int f=usable_history_length-1; f>=0; f--) {
		const FrameObservation& fo = deqFrameObservation[f];

		arma::colvec xt = get_arma_colvec_for_std( fo.vCP_observations ) / ms->get_CP_N();
		arma::colvec yt = get_arma_colvec_for_std( fo.vOB_observations ) / ms->get_OB_delta();

		arma::colvec yt_predicted = lp.get_prediction(hist_state, xt);

		arma::colvec importance_error = (yt_predicted - yt) % ms->get_OB_importance();

		mte = gamma * mte + squared( arma::norm(importance_error, 2) );

		hist_state.state_transition(xt, yt);
	}

	mte *= (1 - gamma) / ((int)ms->vActive_OB_IDs.size());

	return std::make_pair(mte, hist_state);
}


double Model::optimized_compute_mte(
	const std::deque<FrameObservation>& deqFrameObservation,
	double gamma
)
{
	long long int additional_history_length = -1; //initialize to invalid value
	if(cached_mte < 0.0 or cached_gamma != gamma or (int)deqFrameObservation.size() == 0
			or cached_hist_state.isStateFullySetup() == false
			or (additional_history_length = deqFrameObservation.at(0).frame_number - last_frame_number_of_compute_mte) < 0
			or additional_history_length > (int)deqFrameObservation.size()
	)
	{
		std::pair<double, LDS_State> cached_mte_hist_state = compute_mte(deqFrameObservation, gamma);
		cached_mte        = cached_mte_hist_state.first;
		cached_hist_state = cached_mte_hist_state.second;

		cached_gamma = gamma;
		last_frame_number_of_compute_mte = deqFrameObservation.at(0).frame_number;
		return cached_mte;
	}

	assert(last_frame_number_of_compute_mte >= 0
			and last_frame_number_of_compute_mte <= deqFrameObservation.at(0).frame_number);

	assert(gamma == cached_gamma);
	assert(0.0 < gamma and gamma < 1.0);
	cached_mte /= (1 - gamma) / ((int)ms->vActive_OB_IDs.size());

	assert(cached_hist_state.isStateFullySetup());

	assert(additional_history_length >= 0);
	assert(additional_history_length <= (int)deqFrameObservation.size());
	for(int f=additional_history_length-1; f>=0; f--) {
		const FrameObservation& fo = deqFrameObservation[f];

		arma::colvec xt = get_arma_colvec_for_std( fo.vCP_observations ) / ms->get_CP_N();
		arma::colvec yt = get_arma_colvec_for_std( fo.vOB_observations ) / ms->get_OB_delta();

		arma::colvec yt_predicted = lp.get_prediction(cached_hist_state, xt);

		arma::colvec importance_error = (yt_predicted - yt) % ms->get_OB_importance();

		cached_mte = gamma * cached_mte + squared( arma::norm(importance_error, 2) );

		cached_hist_state.state_transition(xt, yt);
	}

	cached_mte *= (1 - gamma) / ((int)ms->vActive_OB_IDs.size());

	last_frame_number_of_compute_mte = deqFrameObservation.at(0).frame_number;
	return cached_mte;
}


double Model::compute_mte_of_history_subrange(
	const std::deque<FrameObservation>& deqFrameObservation,
	long long int offset_from_latest_frame, // >= 0
	long long int num_frames, // not counting frames to setup state
	double gamma
)
{
	assert(ms != 0);
	assert(isDefined());

	assert(offset_from_latest_frame >= 0);
	assert(num_frames >= 0);

	int num_startup_frames = std::max(ms->x_order, ms->y_order);
	long long int num_frames_to_copy = num_frames + num_startup_frames;
	if(offset_from_latest_frame + num_frames_to_copy > (long long int)deqFrameObservation.size())
		return 0.0;

	std::deque<FrameObservation> subrange_deqFrameObservation;
	for(long long int i=0; i<num_frames_to_copy; i++)
		subrange_deqFrameObservation.push_back( deqFrameObservation.at(offset_from_latest_frame + i) );

	std::pair<double, LDS_State> mte_hist_state = compute_mte(subrange_deqFrameObservation, gamma);
	return mte_hist_state.first;
}


} //namespace fctrl


std::ostream& operator<<(std::ostream& os, const fctrl::Model& m)
{
	os << m.q;

	return os;
}
