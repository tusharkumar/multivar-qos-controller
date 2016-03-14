#include <algorithm>
#include <sstream>
#include <iomanip>

#include "debuglog.h"
#include "debug_control.h"
#include "utils.h"
#include "controlparameter.h"
#include "objective.h"
#include "featurecontroller.h"
#include "fctrl_version.h"

namespace fctrl {

std::vector<FCdata *>& access_vFCdata() {
	static std::vector<FCdata *> vFCdata;
	return vFCdata;
}

template<class T>
	void sorted_vector_insert(std::vector<T>& vec, T elem)
{
	//Assumption: vec is sorted in ascending order
	
	for(int i=0; i<(int)vec.size(); i++) {
		if(elem < vec[i]) {
			vec.insert(vec.begin() + i, elem);
			return;
		}
		if(elem == vec[i]) {
			std::cerr << "sorted_vector_insert(): ERROR: Attempting to insert duplicate value in vector: elem=" << elem << std::endl;
			exit(1);
		}
	}

	//elem larger than all vector elements
	vec.push_back(elem);
}


template<class T>
	void vector_delete(std::vector<T>& vec, T elem)
{
	typename std::vector<T>::iterator vi = std::find(vec.begin(), vec.end(), elem);

	if(vi == vec.end()) {
		std::cerr << "vector_delete(): ERROR: Attempting to delete non-existent value in vector: elem=" << elem << std::endl;
		exit(1);
	}

	vec.erase(vi);
}




////////////////// class EntityStatistics ///////////////////

std::string EntityStatistics::print_string(bool bPrintHistograms) const
{
	std::ostringstream oss;
	oss << "MSEQ=" << last_frame_mseq << " satisfied=" << (last_frame_satisfied ? "true" : "false") << " cumulative_MSEQ=" << cumulative_frame_mseq << " cumulative_satisfaction_ratio=" << cumulative_satisfaction_ratio << "\n";

	if(bPrintHistograms) {
		oss << std::setw(10) << " ";
		for(int b=0; b<(int)buckets.size(); b++)
			oss << std::setw(6) << std::setprecision(4) << buckets[b];
		oss << "\n";

		for(std::map<InputChoice, Histogram>::const_iterator cit = mCP_histograms.begin();
				cit != mCP_histograms.end(); cit++) {
			const InputChoice& ic = cit->first;
			const Histogram& h = cit->second;
			assert(h.size() == buckets.size());

			std::string ic_string = "[";
			for(int j=0; j<(int)ic.size(); j++) {
				if(j != 0)
					ic_string += ",";

				char s[100];
				sprintf(s, "%d", ic[j]);
				ic_string += s;
			}
			ic_string += "]";

			oss << std::setw(10) << ic_string;
			for(int b=0; b<(int)h.size(); b++)
				oss << std::setw(6) << h[b];
			oss << "\n";
		}
	}
	return oss.str();
}


void reset_entity_statistics(EntityStatistics& es)
{
	es.num_frames = 0;
	es.last_frame_mseq = -1.0;
	es.last_frame_satisfied = false;
	es.cumulative_frame_mseq = 0.0;
	es.cumulative_satisfaction_ratio = 0.0;
	//es.buckets.clear();
	es.mCP_histograms.clear();
}

void update_entity_statistics(EntityStatistics& es, double latest_mseq, const EntityStatistics::InputChoice& vChoice)
{
	es.last_frame_mseq = latest_mseq;

	es.last_frame_satisfied = (es.last_frame_mseq <= 1.0);

	es.cumulative_frame_mseq = (es.cumulative_frame_mseq * es.num_frames + es.last_frame_mseq) / (es.num_frames + 1);
	es.cumulative_satisfaction_ratio = (es.cumulative_satisfaction_ratio * es.num_frames + (es.last_frame_satisfied ? 1 : 0)) / (es.num_frames + 1);


	if(es.buckets.empty()) {
		const double table[] = {0.2, 0.5, 1.0, 2.0, 3.0, 5.0, 10.0, 20.0, 40.0, 100.0};
		for(int i=0; i<10; i++)
			es.buckets.push_back(table[i]);
	}

	int b_loc=0;
	while(b_loc<(int)es.buckets.size()-1) {
		if(es.last_frame_mseq <= es.buckets[b_loc])
			break;
		b_loc++;
	}

	if(es.mCP_histograms.count(vChoice) == 0)
		es.mCP_histograms[vChoice] = EntityStatistics::Histogram(es.buckets.size(), 0);
	es.mCP_histograms[vChoice].at(b_loc)++;

	es.num_frames++;
}

////////////////// class AllStatistics ///////////////////

std::string AllStatistics::print_all_string(bool bPrintHistograms) const
{
	std::ostringstream oss;
	oss << "frame: " << frame_stats.print_string(bPrintHistograms);

	std::map<OB_ID, EntityStatistics>::const_iterator cit = objective_stats.begin();
	std::map<OB_ID, EntityStatistics>::const_iterator second_cit;

	if(cit == objective_stats.end() or ++(second_cit = cit) == objective_stats.end())
		return oss.str();
			//want at least two objectives when printing behaviors for individual objectives

	while(cit != objective_stats.end()) {
		oss << "\n";
		oss << "objective[" << cit->first << "]: " << (cit->second).print_string(bPrintHistograms);
		cit++;
	}

	return oss.str();
}

////////////////// class FeatureController ///////////////////

FeatureController::FeatureController()
	: fcid(access_vFCdata().size())
{
	access_vFCdata().push_back( new FCdata() );
	std::cout << "Created FeatureController " << fctrl_version << std::endl;
}


void FeatureController::add_control_parameter(const ControlParameter& cp)
{
	if(cp.cpid == -1) {
		std::cerr << "FeatureController::add_control_parameter(): ERROR: Adding undefined ControlParameter (cpid=-1)" << std::endl;
		exit(1);
	}

	FCdata& fcdata = *(access_vFCdata().at(fcid));
	sorted_vector_insert(fcdata.vCP_IDs, cp.cpid);
}


void FeatureController::delete_control_parameter(const ControlParameter& cp)
{
	if(cp.cpid == -1) {
		std::cerr << "FeatureController::delete_control_parameter(): ERROR: Deleting undefined ControlParameter (cpid=-1)" << std::endl;
		exit(1);
	}

	FCdata& fcdata = *(access_vFCdata().at(fcid));
	vector_delete(fcdata.vCP_IDs, cp.cpid);
}


void FeatureController::add_objective(const Objective& ob)
{
	if(ob.obid == -1) {
		std::cerr << "FeatureController::add_objective(): ERROR: Adding undefined Objective (obid=-1)" << std::endl;
		exit(1);
	}

	FCdata& fcdata = *(access_vFCdata().at(fcid));
	sorted_vector_insert(fcdata.vOB_IDs, ob.obid);
}


void FeatureController::delete_objective(const Objective& ob)
{
	if(ob.obid == -1) {
		std::cerr << "FeatureController::delete_objective(): ERROR: Deleting undefined Objective (obid=-1)" << std::endl;
		exit(1);
	}

	FCdata& fcdata = *(access_vFCdata().at(fcid));
	vector_delete(fcdata.vOB_IDs, ob.obid);
}


void FeatureController::set_perception_window_length(int W)
{
	if(W <= 0) {
		std::cerr << "FeatureController::set_perception_window_length(): ERROR: need W > 0. Provided W=" << W << std::endl;
		exit(1);
	}
	FCdata& fcdata = *(access_vFCdata().at(fcid));
	fcdata.W = W;
}

void FeatureController::set_model_x_order(int x_order)
{
	if(x_order < 0) {
		std::cerr << "FeatureController::set_model_x_order(): ERROR: need x_order >= 0. Provided x_order=" << x_order << std::endl;
		exit(1);
	}
	FCdata& fcdata = *(access_vFCdata().at(fcid));
	fcdata.requested_x_order = x_order;
	std::cout << "FeatureController:: Setting x_order=" << x_order << std::endl;
}

void FeatureController::set_model_y_order(int y_order)
{
	if(y_order < 0) {
		std::cerr << "FeatureController::set_model_y_order(): ERROR: need y_order >= 0. Provided y_order=" << y_order << std::endl;
		exit(1);
	}
	FCdata& fcdata = *(access_vFCdata().at(fcid));
	fcdata.requested_y_order = y_order;
	std::cout << "FeatureController:: Setting y_order=" << y_order << std::endl;
}

int FeatureController::get_model_x_order()
{
	FCdata& fcdata = *(access_vFCdata().at(fcid));
	return fcdata.requested_x_order;
}

int FeatureController::get_model_y_order()
{
	FCdata& fcdata = *(access_vFCdata().at(fcid));
	return fcdata.requested_y_order;
}


void FeatureController::frame_transition()
{
	timeval fc_start_time = get_curr_timeval();

	FCdata& fcdata = *(access_vFCdata().at(fcid));

	if(fcdata.frame_number == -1) //startup
		fcdata.reinitialize_setup();

	//Steps:
	// 1. Read Objective values for terminating frame into a FrameObservation (if not first invocation)
	//
	// 2. Reinitialize if required.
	//     - expand/contract prior FrameObservations based on prior cpids and obids
	//
	// 3. Track efficacy of current scheme (metrics etc):
	//      - Determine if current scheme is failing
	//      - Determine if structuring parameters need to be changed.
	//
	//   Method:
	//     - For each y_i, detect oscillation or sluggishness failure, update failure-mode metrics and beta_i
	//     - Update metrics for Model-Estimation-Error => possibly invalidate current linear-model
	//
	// 4. If needed, re-estimate model, re-generate prediction model and re-setup lqr.
	//    a) LLSE 
	//    b) LinearEstimationModel: initial predictor from history
	//    c) Setup LQR, do not solve here
	//
	// 5. Perform state update for completed frame.
	//
	// 6. If LQR defined, apply control inputs, else generate exploratory inputs
	//    a) If LQR already Setup:
	//       - extract control-input sequence => characterize input-boundedness
	//       - re-solve LQR and repeat if necessary: increase N, and/or, adjust R
	//
	//    b) Exploratory inputs that significantly & systematically exercise application X-Y space
	//           
	//
	// 7. Write new ControlParameter values, update into next FrameObservation
	//
	// 8. Increment frame_number

	//Step 1
	std::vector<double> filtered_observations;
	if(fcdata.deqFrameObservation.size() > 0) {
		filtered_observations = fcdata.read_objectives(); //record the observed Ys
	}

	//Step 2
	if(fcdata.wereFixedSettingsAppliedPreviousFrame == false and (int)fcdata.vCP_fixedSetting_pairs.size() != 0) {
#ifdef FC_CORE_DEBUG
		std::cout << "--- SWITCHING to FIXED SETTINGS ---" << std::endl;
#endif //FC_CORE_DEBUG
		fcdata.bReinitialize = true;
	}

	if(fcdata.wereFixedSettingsAppliedPreviousFrame == true  and (int)fcdata.vCP_fixedSetting_pairs.size() == 0) {
#ifdef FC_CORE_DEBUG
		std::cout << "--- SWITCHING to CONTROLLER ---" << std::endl;
#endif //FC_CORE_DEBUG
		fcdata.bReinitialize = true;
	}

	if(fcdata.bReinitialize)
		fcdata.reinitialize_setup();

	//Step 3
	if(fcdata.deqFrameObservation.size() > 0) {
		const FrameObservation& last_frame = fcdata.deqFrameObservation[0];
		assert(filtered_observations.size() == last_frame.vOB_observations.size());

		arma::colvec Y_filtered_measured(filtered_observations.size());
		for(int i=0; i<(int)filtered_observations.size(); i++)
			Y_filtered_measured(i) = filtered_observations[i];

		double last_frame_mseq = compute_MSEQ(Y_filtered_measured, fcdata.ms.get_OB_center(), fcdata.ms.get_OB_delta(), fcdata.ms.get_OB_importance());

		fcdata.deqMSEQ.push_front(last_frame_mseq);
		assert(fcdata.deqMSEQ.size() == fcdata.deqFrameObservation.size());

		EntityStatistics::InputChoice vChoice;
		for(int j=0; j<(int)last_frame.vCP_observations.size(); j++)
			vChoice.push_back( (int)last_frame.vCP_observations[j] );

		update_entity_statistics(fcdata.stats.frame_stats, last_frame_mseq, vChoice);

		arma::colvec last_frame_seq_objectives = compute_SEQ_objectives(Y_filtered_measured, fcdata.ms.get_OB_center(), fcdata.ms.get_OB_delta(), fcdata.ms.get_OB_importance());

		for(int i=0; i<(int)fcdata.ms.vActive_OB_IDs.size(); i++) {
			OB_ID obid = fcdata.ms.vActive_OB_IDs[i];
			if(fcdata.stats.objective_stats.count(obid) == 0)
				fcdata.stats.objective_stats[obid] = EntityStatistics();

			update_entity_statistics(fcdata.stats.objective_stats[obid], last_frame_seq_objectives(i), vChoice);
		}
	}
	
  //Determine whether to run controller or apply fixed settings
	arma::colvec curr_control_input;
	if((int)fcdata.vCP_fixedSetting_pairs.size() == 0) {
		curr_control_input = fcdata.run_controller(fc_start_time);
		fcdata.wereFixedSettingsAppliedPreviousFrame = false;
	}
	else {
		curr_control_input = fcdata.run_fixed_settings();
		fcdata.wereFixedSettingsAppliedPreviousFrame = true;
	}

  //Override control inputs, if desired for this frame
  const arma::colvec CP_N = fcdata.ms.get_CP_N();
  for(size_t p=0; p<fcdata.vCP_overrideSettings_pairs.size(); p++) {
    CP_ID cpid  = fcdata.vCP_overrideSettings_pairs[p].first;
    int   value = fcdata.vCP_overrideSettings_pairs[p].second;

    bool   found = false;
    for(size_t l=0; l<fcdata.ms.vActive_CP_IDs.size(); l++) {
      if(fcdata.ms.vActive_CP_IDs[l] == cpid) {
        assert(value >= -CP_N(l) && value <= CP_N(l));
        curr_control_input(l) = ((double)value)/CP_N(l);
        found = true;
        break;
      }
    }
    assert(found);
  }
  fcdata.override_settings_requested_for_frame = (fcdata.vCP_overrideSettings_pairs.size() != 0);
    //flag for next frame to read

  fcdata.vCP_overrideSettings_pairs.clear();

	//Step 7
	fcdata.write_controlparameters(curr_control_input);
	fcdata.prev_control_input = curr_control_input;
	
	//Step 8
	fcdata.frame_number++;

	timeval fc_end_time = get_curr_timeval();
	fcdata.prev_time = diff_time(fc_start_time, fc_end_time);
#ifdef FC_CORE_DEBUG
	std::cout << "prev_time=" << fcdata.prev_time << std::endl;
#endif //FC_CORE_DEBUG
}


void FeatureController::reinitialize()
{
	FCdata& fcdata = *(access_vFCdata().at(fcid));
	fcdata.bReinitialize = true;
}

long long int FeatureController::get_frame_number() const
{
	FCdata& fcdata = *(access_vFCdata().at(fcid));
	return fcdata.frame_number;
}


void FeatureController::set_frame_budget(double time_budget)
{
	FCdata& fcdata = *(access_vFCdata().at(fcid));
	fcdata.frame_budget = time_budget;
	std::cout << "FeatureController: frame_budget=" << fcdata.frame_budget << std::endl;
}

void FeatureController::set_minimum_history_length(int min_length)
{
	assert(min_length >= 0);
	FCdata& fcdata = *(access_vFCdata().at(fcid));
	fcdata.L_min = min_length;
	std::cout << "FeatureController: L_min=" << fcdata.L_min << std::endl;
}

const AllStatistics& FeatureController::get_allstatistics() const
{
	const FCdata& fcdata = *(access_vFCdata().at(fcid));
	return fcdata.stats;
}


void FeatureController::force_fixed_control_parameters(
	const std::vector< std::pair<CP_ID, int> >& vCP_fixedSetting_pairs
)
{
	FCdata& fcdata = *(access_vFCdata().at(fcid));
	fcdata.vCP_fixedSetting_pairs = vCP_fixedSetting_pairs;
}

void FeatureController::override_control_parameter_in_next_frame(
    CP_ID controlparm_id,
    int   override_value
)
{
	FCdata& fcdata = *(access_vFCdata().at(fcid));
	fcdata.vCP_overrideSettings_pairs.push_back( std::make_pair(controlparm_id, override_value) );
}



////////////////// class FCdata ///////////////////

const double oscillation_d = 0.6;
const double sluggishness_d = 0.6;

void FCdata::reinitialize_setup()
{
#define OVERRIDE_WITH_Y_OPERATING false //FIXME
	ms.flag_subtract_ytm1_in_estimation_for_DC_offset_correction = false; //not OVERRIDE_WITH_Y_OPERATING; //FIXME
	ms.flag_subtract_yoperating_in_estimation_for_DC_offset_correction = false; //OVERRIDE_WITH_Y_OPERATING; //FIXME

	ms.flag_subtract_yobjective_in_estimation_for_DC_offset_correction = false;

#define USE_AFFINE_MODEL true //FIXME
	ms.flag_use_affine_model = USE_AFFINE_MODEL;

	ms.flag_use_adaptive_betas = true;

	ms.active_W = W;
	ms.x_order = requested_x_order;
	ms.y_order = requested_y_order;
	assert(ms.x_order >= 0);
	assert(ms.y_order >= 0);

	if(ms.flag_subtract_ytm1_in_estimation_for_DC_offset_correction == true) {
		if(ms.x_order == 0 and ms.y_order == 0) {
			ms.flag_subtract_ytm1_in_estimation_for_DC_offset_correction = false;
			ms.flag_subtract_yobjective_in_estimation_for_DC_offset_correction = true;
		}
	}

	bool changed_CP_IDs = (ms.vActive_CP_IDs != vCP_IDs);
	bool changed_OB_IDs = (ms.vActive_OB_IDs != vOB_IDs);

	ms.vActive_CP_IDs = vCP_IDs;
	ms.vActive_OB_IDs = vOB_IDs;

	if(changed_CP_IDs or changed_OB_IDs) {
		deqFrameObservation.clear(); //previous data is invalidated
		deqMSEQ.clear();
		cs.reset();
	}

	bReinitialize = false;
	coverage_met = false;
	M.reset();
	prev_newM = false;

	stability_length_metrics.reset();
	coverage_length_metrics.reset();

	prev_apply_controller = false;

	int numXs = ms.vActive_CP_IDs.size();
	diagR = 0.001 * arma::ones<arma::colvec>(numXs);
	c = arma::ones<arma::colvec>((int)ms.vActive_CP_IDs.size());
	f_barely_underconstrained = std::vector<bool>((int)ms.vActive_CP_IDs.size(), false);

	C.reset();
	C_ZO.reset();

	prev_control_input = arma::colvec();

	vOscillationFailureMode.clear();
	vSluggishnessFailureMode.clear();

	ms.estimate_operating_point_from_data_samples(deqFrameObservation, gamma);

	arma::colvec y_max = (ms.get_OB_center() + ms.get_OB_delta()) / ms.get_OB_delta();
	arma::colvec y_min = (ms.get_OB_center() - ms.get_OB_delta()) / ms.get_OB_delta();

	last_correction = arma::zeros<arma::colvec>(ms.vActive_OB_IDs.size());
	betas = arma::ones<arma::colvec>(ms.vActive_OB_IDs.size());
	for(int i=0; i<(int)ms.vActive_OB_IDs.size(); i++) {
		vOscillationFailureMode.push_back ( OscillationFailureMode (y_max(i), y_min(i), ms.active_W, oscillation_d) );
		vSluggishnessFailureMode.push_back( SluggishnessFailureMode(y_max(i), y_min(i), ms.active_W, sluggishness_d) );
	}

	last_exp_frame_number = -1;

  v_SR_level_metrics.clear();

  l_tr                    = 0;

  max_frame_MSEQ          = 1.0;

	num_frames_of_M         = 1;

	theta                   = 0.0;

	forced_length_remaining = 0;
	num_PFE_clusters        = 1;

	d_cdf.clear();


	gamma = 0.0;
	L_gamma = -1;
	candidate_L_gamma = 0;
	prev_frame_ended_PFE = false;

	debt = 0.0;


	//reset stats
	reset_entity_statistics(stats.frame_stats);

	if(changed_OB_IDs)
		stats.objective_stats.clear();

	std::map<OB_ID, EntityStatistics>::iterator it = stats.objective_stats.begin();
	while(it != stats.objective_stats.end()) {
		reset_entity_statistics(it->second);
		it++;
	}


#ifdef FC_CORE_DEBUG
	std::cout << "FCdata:reinitialize_setup() COMPLETE:" << std::endl;
	std::cout << ms << std::endl;
	std::cout << "oscillation_d=" << oscillation_d << " sluggishness_d=" << sluggishness_d << std::endl;
#endif //FC_CORE_DEBUG
}

void FCdata::reset_beta_failuremodes() {
	last_correction = arma::zeros<arma::colvec>(ms.vActive_OB_IDs.size());
	betas = arma::ones<arma::colvec>(ms.vActive_OB_IDs.size());

	assert(vOscillationFailureMode.size() == vSluggishnessFailureMode.size());
	for(int i=0; i<(int)vOscillationFailureMode.size(); i++) {
		vOscillationFailureMode[i].reset();
		vSluggishnessFailureMode[i].reset();
	}
}

std::vector<double> FCdata::read_objectives()
{
	FrameObservation& curr_fo = deqFrameObservation.at(0);
	assert(curr_fo.frame_number == frame_number);
	assert(curr_fo.vOB_observations.size() == 0);

	std::vector<double> filtered_vOB_observations;

	for(int i=0; i<(int)ms.vActive_OB_IDs.size(); i++) {
		OB_ID obid = ms.vActive_OB_IDs[i];
		OBJdata& objdata = access_vOBJdata().at(obid);

		assert(objdata.isDefined);
		objdata.perform_measurement();

		//FIXME
		//double raw_obj_value                     = objdata.read_raw_measurement();
		//curr_fo.vOB_observations.push_back(raw_obj_value); //store

		double sliding_window_averaged_obj_value = objdata.read_sliding_window_averaged_measurement();
		filtered_vOB_observations.push_back(sliding_window_averaged_obj_value); //store
		curr_fo.vOB_observations.push_back(sliding_window_averaged_obj_value); //store
	}
	arma::colvec normalized_observations = (get_arma_colvec_for_std(curr_fo.vOB_observations) - ms.get_OB_center()) / ms.get_OB_delta();
	arma::colvec normalized_filtered_observations = (get_arma_colvec_for_std(filtered_vOB_observations) - ms.get_OB_center()) / ms.get_OB_delta();
#ifdef FC_ANALYSIS_DEBUG
	std::cout << " READ frame-normalized-observations=" << get_std_for_arma_colvec(normalized_observations) << std::endl;
#endif //FC_ANALYSIS_DEBUG
#ifdef FC_CORE_DEBUG
	std::cout << "   normalized-filtered-observations=" << get_std_for_arma_colvec(normalized_filtered_observations) << std::endl;
#endif //FC_CORE_DEBUG

	return filtered_vOB_observations;
}


void FCdata::write_controlparameters(const arma::colvec& curr_control_input)
{
	const std::vector<int> vCP_choices = ms.convert_normalized_to_clipped_CP_choices(curr_control_input);
#ifdef FC_ANALYSIS_DEBUG
	std::cout << "write_controlparameters(): frame_number=" << (frame_number+1) << " curr_control_input=" << curr_control_input
											<< " APPLYING frame-choices=" << vCP_choices << std::endl;
	std::cout << "------------------------------------------------" << std::endl << std::endl;
#endif //FC_ANALYSIS_DEBUG
	assert(vCP_choices.size() == ms.vActive_CP_IDs.size());

	deqFrameObservation.push_front(FrameObservation());
	FrameObservation& next_fo = deqFrameObservation.at(0);
	next_fo.frame_number = frame_number+1;

	if((int)cs.num_dims() == 0)
		cs.reinitialize_dims(ms.vActive_CP_IDs.size());

	for(int j=0; j<(int)ms.vActive_CP_IDs.size(); j++) {
		CP_ID cpid = ms.vActive_CP_IDs[j];
		CPdata& cpdata = access_vCPdata().at(cpid);
		int choice = vCP_choices[j];

		assert(cpdata.isDefined);
		if(cpdata.vCallbacks.at(choice + cpdata.N) == 0) {
			std::cerr << "ControlParameter::ERROR: callback not provided for writing value:"
				<< " cpid = " << cpid << std::endl;
			exit(1);
		}

		cpdata.vCallbacks.at(choice + cpdata.N)(choice); //invoke callback
		next_fo.vCP_observations.push_back(choice); //store
	}

	std::vector<double> v_x(ms.vActive_CP_IDs.size());
	for(int j=0; j<(int)ms.vActive_CP_IDs.size(); j++)
		v_x[j] = next_fo.vCP_observations[j];

	add_sample(v_x, cs.v_means, cs.v_stds, cs.v_maxswings, cs.v_wmins, cs.v_wmaxs, deqFrameObservation.size()-1, gamma);
}



arma::colvec FCdata::run_fixed_settings()
{
	if(vCP_fixedSetting_pairs.size() != ms.vActive_CP_IDs.size()) {
		std::cerr << "FeatureController::ERROR: fixed settings must be provided corresponding to control-parameters"
				<< " associated with FeatureController.\n"
				<< " vCP_fixedSetting_pairs = " << vCP_fixedSetting_pairs << "\n"
				<< " vActive_CP_IDs = " << ms.vActive_CP_IDs << std::endl;
		exit(1);
	}

	arma::colvec curr_control_input = arma::zeros<arma::colvec>((int)ms.vActive_CP_IDs.size());

	arma::colvec Xassigns = arma::zeros<arma::colvec>(ms.vActive_CP_IDs.size());
	for(int p=0; p<(int)vCP_fixedSetting_pairs.size(); p++) {
		std::pair<CP_ID, int> pair = vCP_fixedSetting_pairs[p];
		CP_ID cpid = pair.first;
		int choice = pair.second;

		std::vector<CP_ID>::iterator it = std::find(ms.vActive_CP_IDs.begin(), ms.vActive_CP_IDs.end(), cpid);
		if(it == ms.vActive_CP_IDs.end()) { //not found
			std::cerr << "FeatureController::ERROR: CP_ID=" << cpid << " not associated with controller.\n"
				<< " Attempted to apply fixed settings = " << vCP_fixedSetting_pairs << std::endl;
			exit(1);
		}
		int j = it - ms.vActive_CP_IDs.begin();
		assert(ms.vActive_CP_IDs.at(j) == cpid);

		CPdata& cpdata = access_vCPdata().at(cpid);
		if(choice < -cpdata.N or choice > cpdata.N) {
			std::cerr << "FeatureController::ERROR: fixed choice=" << choice << " for CP_ID=" << cpid
				<< " exceeds bounds" << std::endl;
			exit(1);
		}
		curr_control_input(j) = choice / (double)cpdata.N;
		Xassigns(j) = 1;
	}

	double check_sum = arma::sum(Xassigns);
	if(check_sum != (double)ms.vActive_CP_IDs.size()) {
		std::cerr << "FeatureController::ERROR: fixed settings must be provided corresponding to control-parameters"
				<< " associated with FeatureController.\n"
				<< " vCP_fixedSetting_pairs = " << vCP_fixedSetting_pairs << "\n"
				<< " vActive_CP_IDs = " << ms.vActive_CP_IDs << std::endl;
		exit(1);
	}

	return curr_control_input;
}

arma::colvec FCdata::run_controller(timeval fc_start_time)
{
	adaptive_state_transition();

#ifdef FC_CORE_DEBUG
	std::cout << cs;
#endif //FC_CORE_DEBUG
	check_history_coverage(cs.v_stds, cs.v_maxswings, deqFrameObservation.size(), gamma,
							kappa, coverage_met //outputs
	);
#ifdef FC_ANALYSIS_DEBUG
	std::cout << "coverage_met=" << coverage_met << " kappa=" << kappa << std::endl;
#endif //FC_ANALYSIS_DEBUG

	bool forcedExploration = probabilistic_forced_exploration();
#ifdef FC_CORE_DEBUG
	std::cout << "forcedExploration=" << forcedExploration << std::endl;
#endif //FC_CORE_DEBUG
#ifdef FC_ANALYSIS_DEBUG
	if(forcedExploration)
		std::cout << "----FORCED EXPLORATION----" << std::endl;
#endif //FC_ANALYSIS_DEBUG

	double frame_budget_usable = budget_allocate();

	bool newM = update_model(fc_start_time, frame_budget_usable, forcedExploration);
	if(newM)
		reset_beta_failuremodes();

	if(ms.x_order > 0 or ms.y_order > 0) { //Update LQR controller
		bool newC = update_regulator(fc_start_time, frame_budget_usable);
		if(newC == true) {
			LQR_timestep = 0;
		}
	}
	else if(ms.x_order == 0 and ms.y_order == 0) { //Update zero-order controller
		update_regulatorZO(fc_start_time, frame_budget_usable);
	}

	bool controller_defined = C.isDefined() or C_ZO.isDefined();

	bool apply_controller = controller_defined == true and M.isDefined() == true
						and st.isStateFullySetup() == true and forcedExploration == false;
#ifdef FC_CORE_DEBUG
	std::cout << "apply_controller=" << apply_controller << std::endl;
#endif //FC_CORE_DEBUG
	arma::colvec curr_control_input;
	if(apply_controller) {
		//determine control inputs from LQR

		assert(deqFrameObservation.size() > 0);
		if(C.isDefined()) {
			arma::colvec full_LQR_state = st.get_full_LQR_state_repr();

			LQR_timestep = 0; //FIXME: overrides
			curr_control_input = C.get_control_input(full_LQR_state, LQR_timestep);

			LQR_timestep++;
			if(LQR_timestep > C.N_steps-1) {
				LQR_timestep = C.N_steps-1;
			}
		}
		else if(C_ZO.isDefined()) {
			curr_control_input = C_ZO.get_control_input(st);
		}
		else
		{ assert(0); }
	}
	else {
		//exercise range of control inputs to stimulate learning
		curr_control_input = input_explorer();

		last_exp_frame_number = frame_number;
	}
	assert(curr_control_input.n_rows == ms.vActive_CP_IDs.size());


	resize_history(newM, forcedExploration);


#ifdef FC_ANALYSIS_DEBUG
	std::cout << "resize_history(): history_len=" << deqFrameObservation.size()
		      << " gamma=" << gamma << " L_gamma=" << L_gamma << std::endl;
#endif //FC_ANALYSIS_DEBUG

	prev_newM             = newM;
	prev_apply_controller = apply_controller;

	return curr_control_input;
} //FCdata::run_controller()

bool FCdata::update_model(timeval fc_start_time, double frame_budget_usable, bool forcedExploration)
{
	if(M.isDefined() == true) {
		double M_mte = M.optimized_compute_mte(deqFrameObservation, gamma);
#ifdef FC_CORE_DEBUG
		std::cout << "update_model(): M_mte=" << M_mte << " M.mfe=" << M.mfe << " M_est_framenumber=" << M_est_framenumber << " L_gamma=" << L_gamma << std::endl;
#endif //FC_CORE_DEBUG
		if((M_mte > 10 * M.mfe and frame_number - M_est_framenumber > ms.active_W) ||
		   (M.isBalanced() == false and frame_number - M_est_framenumber > 10 * L_gamma)
		)
		{
			M.reset();
			C.reset();
			C_ZO.reset();
#ifdef FC_ANALYSIS_DEBUG
			std::cout << "----RESET M,C----" << std::endl;
#endif //FC_ANALYSIS_DEBUG
		}
	}

	bool newM = false;
	bool newMp = false;

	try { //llse_and_refine_lambda() produces exception if no llse solution found
		if(coverage_met)
		{
			if(M.isDefined() == false) {
				M = llse_and_refine_lambda();
				M_est_framenumber = frame_number;
#ifdef FC_ANALYSIS_DEBUG
				std::cout << "----APPLIED M(new)----" << std::endl;
#endif //FC_ANALYSIS_DEBUG

				Mp = M;
				Mp_est_framenumber = M_est_framenumber;

				newM = true;
				newMp = true;
				sameModel = true;
			}


			while( diff_time(fc_start_time, get_curr_timeval()) < frame_budget_usable/2 && Mp.isBalanced() == false)
			{
				Mp = llse_and_refine_lambda();
				Mp_est_framenumber = frame_number;
#ifdef FC_ANALYSIS_DEBUG
				std::cout << "----REFINED Mp(balance)----" << std::endl;
#endif //FC_ANALYSIS_DEBUG

				newMp = true;
				sameModel = false;
			}

			if(M.isBalanced() == false and Mp.isBalanced() == true)
			{
				M = Mp;
				M_est_framenumber = Mp_est_framenumber;
#ifdef FC_ANALYSIS_DEBUG
				std::cout << "----APPLIED M(balanced)----" << std::endl;
#endif //FC_ANALYSIS_DEBUG

				newM = true;
				sameModel = true;
			}
		}
	}
	catch(const std::runtime_error& e) { //no llse solution found
#ifdef FC_ANALYSIS_DEBUG
		std::cout << "----LLSE NO SOLUTION FOUND----" << e.what() << std::endl;
#endif //FC_ANALYSIS_DEBUG
		//M and Mp will have their valid values from prior to the excepting llse_and_refine_lambda() call.
		//Continue execution as if nothing happened.
	}

	if(M.isDefined() == true) {
		assert(Mp.isDefined() == true);

		if(newM == false and newMp == false) {
			double M_mte_last_frame = M.compute_mte_of_history_subrange(deqFrameObservation, 0, 1, gamma);
			double Mp_mte_last_frame = Mp.compute_mte_of_history_subrange(deqFrameObservation, 0, 1, gamma);
			double max_mte = std::max(M_mte_last_frame, Mp_mte_last_frame);
			if(max_mte > 0)
				advantage_Mp_over_M = (M_mte_last_frame - Mp_mte_last_frame) / max_mte + gamma * advantage_Mp_over_M;
			else
				advantage_Mp_over_M = 0.0 + gamma * advantage_Mp_over_M;

			subst_history_count++;

			bool subst_coverage_met = coverage_met;

			if(subst_history_count < deqFrameObservation.size()) {
				const FrameObservation& last_fo = deqFrameObservation.at(0);
				std::vector<double> v_x(ms.vActive_CP_IDs.size());
				for(int j=0; j<(int)ms.vActive_CP_IDs.size(); j++)
					v_x[j] = last_fo.vCP_observations[j];

				if(subst_history_count == 1) {
					subst_cs.reinitialize_dims(v_x.size());
				}
				add_sample(v_x, subst_cs.v_means, subst_cs.v_stds, subst_cs.v_maxswings, subst_cs.v_wmins, subst_cs.v_wmaxs, subst_history_count-1, gamma);
				double subst_kappa = 0.0; //dummy
				check_history_coverage(subst_cs.v_stds, subst_cs.v_maxswings, subst_history_count, gamma,
										subst_kappa, subst_coverage_met //outputs
				);
			}

			double threshold = (1.0 - std::pow(gamma, std::max((long long int)L_gamma, subst_history_count))) / (1.0 - gamma) * 0.10;

#ifdef FC_ANALYSIS_DEBUG
			std::cout << "advantage_Mp_over_M=" << advantage_Mp_over_M << " threshold=" << threshold << " subst_coverage_met=" << subst_coverage_met << std::endl;
#endif //FC_ANALYSIS_DEBUG

			if(advantage_Mp_over_M > threshold and subst_coverage_met) {
				M = Mp;
				M_est_framenumber = Mp_est_framenumber;
#ifdef FC_ANALYSIS_DEBUG
				std::cout << "----APPLIED M(better)----" << std::endl;
#endif //FC_ANALYSIS_DEBUG

				newM = true;
				sameModel = true;
			}
		}

		try { //llse_and_refine_lambda() produces exception if no llse solution found
			if( diff_time(fc_start_time, get_curr_timeval()) < frame_budget_usable && coverage_met and forcedExploration == false
					and C.isDefined() and C_refine_timestep == frame_number - 1
					and (Mp_est_framenumber <= frame_number - L_gamma or sameModel == true)
			)
			{
				Mp = llse_and_refine_lambda();
				Mp_est_framenumber = frame_number;
#ifdef FC_ANALYSIS_DEBUG
				std::cout << "----REFINED Mp(better)----" << std::endl;
#endif //FC_ANALYSIS_DEBUG

				newMp = true;
				sameModel = false;
			}
		}
		catch(const std::runtime_error& e) { //no llse solution found
#ifdef FC_ANALYSIS_DEBUG
			std::cout << "----LLSE NO SOLUTION FOUND: during Mp REFINEMENT----: " << e.what() << std::endl;
#endif //FC_ANALYSIS_DEBUG
			//M and Mp will have their valid values from prior to the excepting llse_and_refine_lambda() call.
			//Continue execution as if nothing happened.
		}
	}

	if(newM or newMp) {
		advantage_Mp_over_M = 0;
		subst_history_count = 0;
	}

	if(newM == true) {
		C.reset();
		C_ZO.reset();
#ifdef FC_ANALYSIS_DEBUG
		std::cout << "----RESET C(newM)----" << std::endl;
#endif //FC_ANALYSIS_DEBUG
	}

	return newM;
}


bool FCdata::update_regulator(timeval fc_start_time, double frame_budget_usable)
{
	bool newC = false;
	if(M.isDefined() == true and (diff_time(fc_start_time, get_curr_timeval()) < frame_budget_usable || C.isDefined() == false)) {

		bool any_C_refinement = false;

		if(C.isDefined() == false) {
			C.design_controller(M.lp, 1, diagR);
			C_refine_timestep = frame_number;
			newC = true;
			any_C_refinement = true;
			f_barely_underconstrained = std::vector<bool>((int)ms.vActive_CP_IDs.size(), false);


			st = LDS_State(M.lp, M.model_Y_operating); //reset state
			//initialize state from available history
			int num_init_steps = std::min( std::max(ms.x_order, ms.y_order), (int)deqFrameObservation.size() );
			arma::colvec CP_N = ms.get_CP_N();
			arma::colvec OB_delta = ms.get_OB_delta();
			for(int k=num_init_steps-1; k>=0; k--) {
				arma::colvec x_tmk = get_arma_colvec_for_std( deqFrameObservation[k].vCP_observations ) / CP_N;
				arma::colvec y_tmk = get_arma_colvec_for_std( deqFrameObservation[k].vOB_observations ) / OB_delta;

				st.state_transition(x_tmk, y_tmk); //betas reset to unity, so no adaptive transition needed
			}
		}


		assert(M.isDefined());
		assert(C.isDefined());

		arma::colvec CP_N = ms.get_CP_N();
		assert(CP_N.n_rows == ms.vActive_CP_IDs.size());

		if( st.isStateFullySetup() and diff_time(fc_start_time, get_curr_timeval()) < frame_budget_usable ) {

			arma::colvec full_LQR_state = st.get_full_LQR_state_repr();
			assert(full_LQR_state.n_rows > 0); //state is fully setup now


			arma::colvec diagR_prev = diagR;
			arma::colvec x_prev = C.get_control_input(full_LQR_state, 0);


			for(int j=0; j<(int)ms.vActive_CP_IDs.size(); j++) {
				diagR(j) = initial_refinement( f_barely_underconstrained[j], std::abs(x_prev(j)), diagR_prev(j), c(j) );
			}

			//C.increment_controller_solution(diagR); //changed v23
			C.design_controller(M.lp, 1, diagR);
			C_refine_timestep = frame_number;

			any_C_refinement = true;

			arma::colvec diagR_refined = diagR;
			arma::colvec x_refined     = C.get_control_input(full_LQR_state, 0);

#ifdef FC_CORE_DEBUG
			std::cout << "  -- after initial_refinement()" << std::endl;
			std::cout << "diagR_prev=" << diagR_prev;
			std::cout << "x_prev=" << x_prev;
			std::cout << "diagR_refined=" << diagR_refined;
			std::cout << "x_refined=" << x_refined;
			std::cout << std::endl;
#endif //FC_CORE_DEBUG

			
			std::vector<bool> f_term((int)ms.vActive_CP_IDs.size(), false);

			bool all_term;
			do {
				arma::colvec proj_norm_y_t = M.lp.get_prediction(st, x_refined);
				arma::colvec y_norm_center = ms.get_OB_center() / ms.get_OB_delta();
				arma::colvec proj_y_t_error = proj_norm_y_t - y_norm_center;
				double trajectory_tracking_error = arma::dot(proj_y_t_error,  C.Q * proj_y_t_error);


				for(int j=0; j<(int)ms.vActive_CP_IDs.size(); j++) {
					double diagR_j                     = diagR(j);
					bool   f_term_j                    = f_term[j];
					bool   f_barely_underconstrained_j = f_barely_underconstrained[j];
					double c_j                         = c(j);

					refine_input_cost(
						diagR_j, f_term_j, f_barely_underconstrained_j, c_j, //Outputs, except f_barely_underconstrained_j and c_j are also inputs
						std::abs(x_prev(j)), std::abs(x_refined(j)), diagR_prev(j), diagR_refined(j), trajectory_tracking_error //Inputs
					);

					diagR(j)                     = diagR_j;
					f_term[j]                    = f_term_j;
					f_barely_underconstrained[j] = f_barely_underconstrained_j;
					c(j)                         = c_j;
				}

				//C.increment_controller_solution(diagR); //changed v23
				C.design_controller(M.lp, 1, diagR);
				C_refine_timestep = frame_number;
				
				diagR_prev = diagR_refined;
				x_prev     = x_refined;

				diagR_refined = diagR;
				x_refined     = C.get_control_input(full_LQR_state, 0);

				all_term = true;
				for(int j=0; j<(int)ms.vActive_CP_IDs.size(); j++)
					all_term = all_term && f_term[j];

#ifdef FC_CORE_DEBUG
				std::cout << "  -- after refine_input_cost()" << std::endl;
				std::cout << "f_barely_underconstrained=" << f_barely_underconstrained << " f_term=" << f_term << std::endl;
				std::cout << "  -- shift -- " << std::endl;
				std::cout << "diagR_prev=" << diagR_prev;
				std::cout << "x_prev=" << x_prev;
				std::cout << "diagR_refined=" << diagR_refined;
				std::cout << "x_refined=" << x_refined;
				std::cout << std::endl;
#endif //FC_CORE_DEBUG

			} while( all_term == false && diff_time(fc_start_time, get_curr_timeval()) < frame_budget_usable );

		}


		if(any_C_refinement) {
#ifdef FC_CORE_DEBUG
			std::cout << "solve_lqr(): N_steps=" << C.N_steps << " list_r_t=" << C.list_r_t << std::endl;
			std::cout << "A = " << C.A << std::endl;
			std::cout << "B = " << C.B << std::endl;
			std::cout << "C = " << C.C << std::endl;
#endif //FC_CORE_DEBUG

#ifdef FC_ANALYSIS_DEBUG
			std::cout << "---------- LQR SOLUTION ------------" << std::endl;
#endif //FC_ANALYSIS_DEBUG
#ifdef FC_CORE_DEBUG
			std::cout << "list_W_t =" << C.list_W_t << std::endl;
			std::cout << "list_K_t =" << C.list_K_t << std::endl;
			std::cout << "list_Kv_t=" << C.list_Kv_t << std::endl;
			std::cout << "list_v_t=" << C.list_v_t << std::endl;
#endif //FC_CORE_DEBUG
		}

	}

	return newC;
}


const double required_std_coverage_fraction = 0.5;

double get_required_maxswing_coverage_fraction(double history_gamma)
{ return (1.0 + history_gamma) / 2.0; }

void FCdata::check_history_coverage(
		const std::vector<double>& v_history_stds,
		const std::vector<double>& v_history_maxswings,
		std::size_t                history_len,
		double                     history_gamma,
		double&                    history_kappa,       //output
		bool&                      history_coverage_met //output
)
{
	if(history_len <= (std::size_t)ms.x_order or history_len <= (std::size_t)ms.y_order
		or history_len == 0
		or history_len <= (std::size_t)ms.active_W //readings possibly not stable yet
	)
	{
		history_coverage_met = false;
		history_kappa = 0.0;
		return;
	}

	arma::colvec CP_N = ms.get_CP_N();

	assert(v_history_stds.size() == ms.vActive_CP_IDs.size());
	assert(v_history_maxswings.size() == ms.vActive_CP_IDs.size());

	assert(0 < history_gamma && history_gamma < 1);

	double required_maxswing_coverage_fraction = get_required_maxswing_coverage_fraction(history_gamma);

#ifdef FC_CORE_DEBUG
	std::cout << "required_std_coverage_fraction=" << required_std_coverage_fraction
	          << " required_maxswing_coverage_fraction=" << required_maxswing_coverage_fraction << std::endl;
#endif //FC_CORE_DEBUG

	int num_spanning_dims = 0;
	for(int j=0; j<(int)v_history_stds.size(); j++) {
		double normalized_range_j = (v_history_stds[j] / CP_N(j));
		double normalized_maxswing_j = (v_history_maxswings[j] / CP_N(j));

#ifdef FC_CORE_DEBUG
		std::cout << "normalized_range_" << j << "=" << normalized_range_j
		          << " normalized_maxswing_" << j << "=" << normalized_maxswing_j << std::endl;
#endif //FC_CORE_DEBUG
		if(normalized_range_j >= required_std_coverage_fraction &&
		   normalized_maxswing_j >= required_maxswing_coverage_fraction)
		{ num_spanning_dims++; }
	}

	history_coverage_met = (num_spanning_dims == (int)ms.vActive_CP_IDs.size());
	history_kappa        = ((double)num_spanning_dims) / ms.vActive_CP_IDs.size();
}

Model FCdata::llse_and_refine_lambda()
{
	assert(0.0 < gamma and gamma < 1.0);
	Model m(&ms);
	m.llse_estimate(deqFrameObservation, lambda_next, gamma);
	lambda_next = m.refine_lambda();

#ifdef FC_CORE_DEBUG
	std::cout << "lambda_next=" << lambda_next << std::endl;
#endif //FC_CORE_DEBUG

	return m;
}


bool FCdata::force_model_estimation(long long int last_model_estimation_timestamp)
{
	assert(last_model_estimation_timestamp >= 0);
	long long int age = frame_number - last_model_estimation_timestamp;
	assert(age >= 0);

	if(age < ms.active_W)
		return false;

	double threshold = ms.active_W / (double)age;

	double r = (double)std::rand()/(double)RAND_MAX; //random number uniformly between 0.0 and 1.0
	return (r >= threshold);
}


double FCdata::initial_refinement(
	double f_barely_underconstrained_j,
	double abs_x_prev_j,
	double diagR_prev_j,
	double c_j
)
{
	assert(abs_x_prev_j >= 0);
	if(f_barely_underconstrained_j == false)
		c_j = 1.0;

	double diagR_j;
	if( abs_x_prev_j > 1.0 )
		diagR_j = increase_input_cost(diagR_prev_j, c_j);
	else
		diagR_j = decrease_input_cost(diagR_prev_j, c_j);

	return diagR_j;
}



void FCdata::refine_input_cost(
		//Outputs
	double& diagR_j,
	bool  & f_term_j,
	bool  & f_barely_underconstrained_j, //also an input
	double& c_j,                         //also an input
		//Inputs
	double  abs_x_prev_j,
	double  abs_x_refined_j,
	double  diagR_prev_j,
	double  diagR_refined_j,
	double  trajectory_tracking_error
)
{
	if(f_barely_underconstrained_j == false)
		c_j = 1;

	assert(abs_x_prev_j >= 0);
	assert(abs_x_refined_j >= 0);

	if(f_term_j == true) {
		assert(f_barely_underconstrained_j == false);
		evaluate_terminated_refinement(diagR_j, f_term_j, abs_x_prev_j, abs_x_refined_j, diagR_prev_j, diagR_refined_j, c_j);
		return;
	}

	assert(f_term_j == false);
	assert(diagR_prev_j != diagR_refined_j);

	double prod = (diagR_prev_j - diagR_refined_j) * (abs_x_prev_j - abs_x_refined_j);

	if(abs_x_prev_j <= 1.0 and abs_x_refined_j <= 1.0) {
		if( practically_unchanged(abs_x_prev_j, abs_x_refined_j) ) {
#ifdef FC_CORE_DEBUG
			std::cout << "trajectory_tracking_error=" << trajectory_tracking_error << " ";
#endif //FC_CORE_DEBUG
			if( practically_unchanged(abs_x_refined_j, 0.0) and 0.5*0.5*diagR_refined_j > 10 * trajectory_tracking_error ) {
				diagR_j = decrease_input_cost(diagR_refined_j, c_j);
#ifdef FC_CORE_DEBUG
				std::cout << "R down" << std::endl;
#endif //FC_CORE_DEBUG
			}
			else {
				diagR_j = diagR_prev_j;
				f_barely_underconstrained_j = false;
				f_term_j = true;
#ifdef FC_CORE_DEBUG
				std::cout << "R revert term" << std::endl;
#endif //FC_CORE_DEBUG
			}
		}
		else if(prod > 0) {
			diagR_j = diagR_prev_j;
			f_barely_underconstrained_j = false;
		}
		else if(prod < 0) {
			diagR_j = decrease_input_cost(diagR_refined_j, c_j);
		}
		else
		{
			std::cout << "Unexpected prod=" << prod << " diagR_prev_j=" << diagR_prev_j
				<< " diagR_refined_j=" << diagR_refined_j << " abs_x_prev_j=" << abs_x_prev_j
				<< " abs_x_refined_j=" << abs_x_refined_j << std::endl;
			assert(0);
		}
	}

	else if(abs_x_prev_j <= 1.0 and abs_x_refined_j > 1.0) {
		if(diagR_prev_j > diagR_refined_j) {
			f_barely_underconstrained_j = true;
			c_j++;
			diagR_j = increase_input_cost(diagR_refined_j, c_j);
		}
		else if(diagR_prev_j < diagR_refined_j) {
			diagR_j = diagR_prev_j;
			f_barely_underconstrained_j = false;
		}
		else
		{ assert(0); }
	}

	else if(abs_x_prev_j > 1.0 and abs_x_refined_j <= 1.0) {
		if(diagR_prev_j < diagR_refined_j) {
			f_barely_underconstrained_j = true;
			c_j++;
			diagR_j = decrease_input_cost(diagR_refined_j, c_j);
		}
		else if(diagR_prev_j > diagR_refined_j) {
			diagR_j = diagR_prev_j;
			f_barely_underconstrained_j = false;
		}
		else
		{ assert(0); }
	}

	else if(abs_x_prev_j > 1.0 and abs_x_refined_j > 1.0) {
		if(prod <= 0) {
			diagR_j = increase_input_cost(diagR_refined_j, c_j);
		}
		else {
			diagR_j = diagR_prev_j;
			f_barely_underconstrained_j = false;
		}
	}

	else
	{ assert(0); }
}


void FCdata::evaluate_terminated_refinement(
		//Outputs
	double& diagR_j,
	bool  & f_term_j, //also an input
		//Inputs
	double  abs_x_prev_j,
	double  abs_x_refined_j,
	double  diagR_prev_j,
	double  diagR_refined_j,
	double  c_j
)
{
	assert(abs_x_prev_j >= 0);
	assert(abs_x_refined_j >= 0);
	assert(f_term_j == true);

	if(diagR_prev_j != diagR_refined_j) {
		diagR_j = diagR_refined_j;
	}
	else {
		if(practically_unchanged(abs_x_prev_j, abs_x_refined_j) == true) {
			diagR_j = diagR_refined_j;
		}
		else {
			diagR_j = decrease_input_cost(diagR_refined_j, c_j);
			f_term_j = false;
		}
	}
}


double FCdata::increase_input_cost(double diagR_j_before, double c_j)
{
	double diagR_j_after = diagR_j_before * (1 + 1.0/c_j);
	return diagR_j_after;
}

double FCdata::decrease_input_cost(double diagR_j_before, double c_j)
{
	double diagR_j_after = diagR_j_before / (1 + 1.0/c_j);
	return diagR_j_after;
}

bool FCdata::practically_unchanged(double abs_x_prev_j, double abs_x_refined_j) {
	bool result = std::abs(abs_x_prev_j - abs_x_refined_j) * 10 * ms.active_W < 0.5;
#ifdef FC_CORE_DEBUG
	std::cout << "practically_unchanged=" << result << " for " << " abs_x_prev_j=" << abs_x_prev_j << " abs_x_refined_j=" << abs_x_refined_j << std::endl;
#endif //FC_CORE_DEBUG
	
	return result;
}

bool FCdata::probabilistic_forced_exploration() {
	if(frame_number == -1) //startup
		return false;

  assert(deqMSEQ.size() > 0);
  max_frame_MSEQ = std::max(max_frame_MSEQ, deqMSEQ[0]);
  assert(max_frame_MSEQ >= 1.0);

  static const double log_level = log(FCdata::SR_level);

  double log_max_frame_MSEQ = log(max_frame_MSEQ);
  int new_l_max = int(log_max_frame_MSEQ / log_level + 1.00); //conservative approx ceil
  assert(new_l_max >= 0);

  if(new_l_max >= (int)v_SR_level_metrics.size())
    v_SR_level_metrics.resize(new_l_max + 1);

	if(prev_newM and num_frames_of_M > 1) {
			//save any unsaved metrics for the last M
		if(num_frames_of_M % (10 * ms.active_W) != 0) {
			update_achievability();
		}

			//initially, assume the best possible performance for the new M
    for(size_t l=0; l<v_SR_level_metrics.size(); l++)
      v_SR_level_metrics[l].reset_for_new_M();

		num_frames_of_M = 1;
	}

	if(prev_apply_controller) {
		if(frame_number - M_est_framenumber >= ms.active_W) {
      double k = 1.0;
      for(size_t l=0; l<v_SR_level_metrics.size(); l++) {
        v_SR_level_metrics[l].add_SR_sample(num_frames_of_M,
                                            (deqMSEQ[0] <= k),
                                            10 * ms.active_W);
        k *= FCdata::SR_level;
      }
			num_frames_of_M++;
		}

		//update achievability stats every 10 * W model-driven frames
		if(num_frames_of_M % (10 * ms.active_W) == 0) {
			update_achievability();
		}

    // k_tr : k that maximizes benefit.
    // l_tr : the index of k_tr
/*    if(l_tr > 0 &&
       (1 - v_SR_level_metrics[l_tr-1].achievable_SR) <= (1 - v_SR_level_metrics[l_tr].achievable_SR) * 1.50)
    { l_tr--; }
    else if(l_tr+1 < v_SR_level_metrics.size() && 
       (1 - v_SR_level_metrics[l_tr].achievable_SR) > (1 - v_SR_level_metrics[l_tr+1].achievable_SR) * 1.50)
    { l_tr++; }
    */

    l_tr = 0;
    double improvement_k_tr = 0.0;
    for(size_t l=0; l+1<v_SR_level_metrics.size(); l++) {
      double improvement_k = (v_SR_level_metrics[l+1].achievable_SR - v_SR_level_metrics[l].achievable_SR) / (l + 1);
      if(improvement_k > improvement_k_tr) {
        l_tr = l;
        improvement_k_tr = improvement_k;
      }
    }

    double achievable_SR_k  = v_SR_level_metrics[l_tr].achievable_SR;
    double estimated_SR_M_k = v_SR_level_metrics[l_tr].estimated_SR_M;
    double current_SR_M_k   = v_SR_level_metrics[l_tr].current_SR_M;

		theta = 0.2 /* FIXME 0.5*/ * (1 - current_SR_M_k) *
									std::max( (achievable_SR_k - estimated_SR_M_k) / (achievable_SR_k + 0.01), 0.0 ) + 0.01;
	}

#ifdef FC_ANALYSIS_DEBUG
  double achievable_SR  = v_SR_level_metrics[0].achievable_SR;
  double estimated_SR_M = v_SR_level_metrics[0].estimated_SR_M;
  double current_SR_M   = v_SR_level_metrics[0].current_SR_M;

	std::cout << "l_tr=" << l_tr << " achievable_SR=" << achievable_SR
			  << " estimated_SR_M=" << estimated_SR_M << " current_SR_M=" << current_SR_M << std::endl;

  std::cout << "l    =";
  for(size_t l=0; l<v_SR_level_metrics.size(); l++)
    std::cout << "  " << std::setw(8) << l;
  std::cout << std::endl;

  std::cout << "aSR  =";
  for(size_t l=0; l<v_SR_level_metrics.size(); l++)
    std::cout << "  " << std::setw(8) << v_SR_level_metrics[l].achievable_SR;
  std::cout << std::endl;

  std::cout << "eSRM =";
  for(size_t l=0; l<v_SR_level_metrics.size(); l++)
    std::cout << "  " << std::setw(8) << v_SR_level_metrics[l].estimated_SR_M;
  std::cout << std::endl;

  std::cout << "cSRM =";
  for(size_t l=0; l<v_SR_level_metrics.size(); l++)
    std::cout << "  " << std::setw(8) << v_SR_level_metrics[l].current_SR_M;
  std::cout << std::endl;
  std::cout << std::endl;
#endif //FC_ANALYSIS_DEBUG

	if(forced_length_remaining > 0) { //ongoing exploration cluster
#ifdef FC_CORE_DEBUG
		std::cout << "forced_length_remaining=" << forced_length_remaining << std::endl;
#endif //FC_CORE_DEBUG

    //indefinitely extend only the last PFE cluster in a PFE group until coverage achieved
    if(cluster_index_in_PFE_group < num_PFE_clusters ||
       forced_length_remaining > 1 ||
       coverage_met) //extend PFE until coverage met
    {
			forced_length_remaining--;
    }

    assert(cluster_index_in_PFE_group <= num_PFE_clusters);
    //reset the cluster index to indicate that a PFE group has completed
    if(cluster_index_in_PFE_group == num_PFE_clusters && forced_length_remaining == 0)
      cluster_index_in_PFE_group = 0;

		return true;
	}

	if(coverage_met) //forced exploration not needed
		return false;

	size_t L_PFE = estimate_shortest_PFE_cluster_length();

	double prev_q = -1; //initialize to invalid value
	double q = -1;
	const double threshold = theta / 5.0;

	if(num_PFE_clusters > L_PFE)
		num_PFE_clusters = L_PFE;

	double expected_d;
	int sampled_d;
	assert(num_PFE_clusters >= 1);
	assert(L_PFE >= 1);
	while(1)
	{
		const int d_peak = L_PFE / num_PFE_clusters;

		const int d_high = std::min(3 * d_peak, L_gamma);

		adjust_cluster_length_distribution_and_sample(d_peak, d_high, expected_d, sampled_d);
		
		q = theta / ( expected_d * (1 - theta) + theta);

#ifdef FC_ANALYSIS_DEBUG
		std::cout << "Tried: num_PFE_clusters=" << num_PFE_clusters << " d_peak=" << d_peak
					<< " d_high=" << d_high << " threshold=" << threshold << " q=" << q << std::endl;
#endif //FC_ANALYSIS_DEBUG

		if(prev_q >= 0 &&
		   ((prev_q < threshold && q >= threshold) || (prev_q >= threshold && q < threshold))
		)
		{ break; }

		if(q < threshold) {
			if(num_PFE_clusters == L_PFE)
				break;
			num_PFE_clusters++;
		}
		else {
			if(num_PFE_clusters == 1)
				break;
			num_PFE_clusters--;
		}

		prev_q = q;
	}

#ifdef FC_ANALYSIS_DEBUG
	std::cout << "L_PFE=" << L_PFE << " num_PFE_clusters=" << num_PFE_clusters
				<< " expected_d=" << expected_d << " sampled_d=" << sampled_d << " q=" << q << std::endl;
#endif //FC_ANALYSIS_DEBUG

	double r = (double)std::rand()/(double)RAND_MAX; //random number uniformly between 0.0 and 1.0
	if(q <= r)
		return false;

	assert(sampled_d >= 1);

  //start the next PFE cluster
	forced_length_remaining = sampled_d - 1;
  cluster_index_in_PFE_group = std::min(cluster_index_in_PFE_group+1, num_PFE_clusters);
#ifdef FC_ANALYSIS_DEBUG
  std::cout << "Starting cluster_index_in_PFE_group=" << cluster_index_in_PFE_group << std::endl;
#endif //FC_ANALYSIS_DEBUG
	return true;
}

void FCdata::update_achievability()
{
  double k = 1.0;
  for(size_t l=0; l<v_SR_level_metrics.size(); l++) {
#ifdef FC_CORE_DEBUG
    std::cout << "--- Level l=" << l << " k=" << k << std::endl;
#endif //FC_CORE_DEBUG
    v_SR_level_metrics[l].update_achievability();
    k *= FCdata::SR_level;
  }
}

size_t FCdata::estimate_shortest_PFE_cluster_length()
{
	arma::colvec CP_N = ms.get_CP_N();
	assert(CP_N.n_rows == cs.v_stds.size());

	size_t j_o = 0;
	double min_norm_std = cs.v_stds.at(j_o) / CP_N(j_o);

	for(size_t j=1; j<cs.v_stds.size(); j++) {
		double norm_std_j = cs.v_stds[j] / CP_N(j);
		if(norm_std_j < min_norm_std) {
			j_o = j;
			min_norm_std = norm_std_j;
		}
	}
	double min_norm_mean = cs.v_means.at(j_o) / CP_N(j_o);

	const int histLength = deqFrameObservation.size();

	int lower = 1;
	int upper = int(L_gamma / 2);
	while(lower < upper) {
		int L = int((lower + upper) / 2);
		int R = 0;
		if(L + histLength <= L_gamma)
			R = histLength;
		else
			R = L_gamma - L;

		const double one_by_one_minus_gamma_R_plus_L = 1.0 / (1.0 - std::pow(gamma, R+L));
		const double gamma_L = std::pow(gamma, L);
		const double one_minus_gamma_L = 1.0 - gamma_L;
		const double one_minus_gamma_R = 1.0 - std::pow(gamma, R);
		const double one_minus_gamma_L_by_one_plus_gamma = one_minus_gamma_L / (1.0 + gamma);
		const double one_minus_gamma = 1.0 - gamma;

		double norm_mean_estimate =
			one_by_one_minus_gamma_R_plus_L * (one_minus_gamma_L_by_one_plus_gamma * one_minus_gamma * 1.0 +
					                           gamma_L * one_minus_gamma_R * min_norm_mean);

		double one_minus_norm_mean_estimate = 1.0 - norm_mean_estimate;
		double one_plus_norm_mean_estimate  = 1.0 + norm_mean_estimate;
		double min_norm_mean_minus_norm_mean_estimate = min_norm_mean - norm_mean_estimate;
		double min_norm_mean_minus_norm_mean_estimate_squared = min_norm_mean_minus_norm_mean_estimate *
																min_norm_mean_minus_norm_mean_estimate;
		double norm_var_estimate =
			one_by_one_minus_gamma_R_plus_L * (one_minus_gamma_L_by_one_plus_gamma *
											   (one_minus_norm_mean_estimate * one_minus_norm_mean_estimate +
												gamma * one_plus_norm_mean_estimate * one_plus_norm_mean_estimate)
											   + gamma_L * one_minus_gamma_R * (min_norm_std * min_norm_std +
											                                    gamma_L * min_norm_mean_minus_norm_mean_estimate_squared));

		if(norm_var_estimate >= 0.5 * 0.5)
			upper = L;
		else if(lower + 1 == upper)
			lower = upper;
		else
			lower = L;
	}
	size_t L_PFE = size_t(upper);
	return L_PFE;
}

void FCdata::adjust_cluster_length_distribution_and_sample(
		const int d_peak,
		const int d_high,
		double& expected_d,
		int& sampled_d
)
{
	assert(d_high > d_peak);
	//const int d_low  = 0;
	assert(d_peak >= 1);
	const int d_low = d_peak - 1;

	if((int)d_cdf.size() != d_high + 1)
		d_cdf.resize(d_high+1);

	double sum_exp_d  = 0;
	double sum_height = 0;
	for(int d=d_low+1; d<=d_peak; d++) {
		double height = d - d_low;
		sum_exp_d  += d * height;
		sum_height += height;

		d_cdf[d] = sum_height;
	}
	double slope = double(d_peak) / (d_high - d_peak);
	for(int d=d_peak+1; d<=d_high; d++) {
		double height = (d_high - d) * slope;
		sum_exp_d  += d * height;
		sum_height += height;

		d_cdf[d] = sum_height;
	}

	expected_d = sum_exp_d / sum_height;

	double r = (double)std::rand()/(double)RAND_MAX; //random number uniformly between 0.0 and 1.0

	double target_cdf_val = r * sum_height;

	int d=d_low;
	for(d=d_low+1; d<=d_high; d++)
		if(d_cdf[d] > target_cdf_val)
			break;

	assert(d_low < d and d < d_high);

	sampled_d = d;
}


void FCdata::adaptive_state_transition()
{
	if(last_exp_frame_number >= 0 && last_exp_frame_number < frame_number - 1) { // i.e., regulator drove inputs in the last frame
		assert(deqFrameObservation.size() > 0);
		arma::colvec x_tm1 = get_arma_colvec_for_std( deqFrameObservation[0].vCP_observations ) / ms.get_CP_N();
		arma::colvec y_tm1 = get_arma_colvec_for_std( deqFrameObservation[0].vOB_observations ) / ms.get_OB_delta();

		assert(vOscillationFailureMode.size() == vSluggishnessFailureMode.size());
		assert(vOscillationFailureMode.size() == y_tm1.n_rows);

    if(override_settings_requested_for_frame) {
      std::cout << "Skipping adaptive betas: override_settings_requested_for_frame = true" << std::endl;
    }

		bool skip_beta_update = (frame_number - last_exp_frame_number <= ms.active_W)
                            || override_settings_requested_for_frame;
			//skip for W frames after any exploration to filter out effects of exploration on y_tm1,
      //or if control inputs were overridden for the just completed frame, making input trends unreliable this frame

		if(skip_beta_update) {
#ifdef FC_CORE_DEBUG
			std::cout << "skipping beta update" << betas << std::endl;
#endif //FC_CORE_DEBUG
		}
		else {
			for(int i=0; i<(int)vOscillationFailureMode.size(); i++) {
				double correction_OscFailure = vOscillationFailureMode[i].test_failure( y_tm1(i) );
				double correction_SlgFailure = vSluggishnessFailureMode[i].test_failure( y_tm1(i) );

				//FIXME: slow down the oscillation correction, additional corrections will occur subsequently
				//if(correction_OscFailure > 0.0)
				//	correction_OscFailure = std::sqrt(correction_OscFailure);

				if(correction_SlgFailure > 0.0 and correction_OscFailure > 0.0)
          correction_OscFailure = -1;
          //if both apply, prefer sluggishness as its correction is typically much smaller,
          // subsequent oscillation can fix a mistake quickly

				assert(correction_OscFailure == -1 or correction_SlgFailure == -1);
					//atmost one failure mode should be applied at current frame for a given objective

				if(correction_OscFailure > 0.0) {
#ifdef FC_ANALYSIS_DEBUG
					std::cout << "----RESCALING for Oscillation: [" << i << "] correction=" << 
						correction_OscFailure << std::endl;
#endif //FC_ANALYSIS_DEBUG
					betas(i) *= correction_OscFailure;
					last_correction(i) = correction_OscFailure;
#ifdef FC_CORE_DEBUG
					std::cout << "betas=" << betas << std::endl;
#endif //FC_CORE_DEBUG
				}

				if(correction_SlgFailure > 0.0) {
					//indicates a need for Sluggishness correction, now need to determine factor
					double& last_corr = last_correction(i);
#if 0
					if(last_corr == 0) //no previous correction
						last_corr = 1.20; //default conservative initial factor
					else if(last_corr < 1.0) //last correction was for Oscillation
						last_corr = 1 / last_corr * 0.80; //80% smaller counter-action
					else //last correction was for Sluggishness
						last_corr *= 1.20; //ramp up gradually over multiple corrections
#else
          if(last_corr == 0) //no previous correction
            last_corr = 1.20; //default initial

          if(betas(i) < 1.0) //oscillation correction dominates
            last_corr = 1.0 / std::sqrt( betas(i) );
          else
            last_corr *= 1.20; //ramp up correction in subsequent corrections
#endif
					correction_SlgFailure = last_corr;
#ifdef FC_ANALYSIS_DEBUG
					std::cout << "----RESCALING for Sluggishness: [" << i << "] correction=" << 
						correction_SlgFailure << std::endl;
#endif //FC_ANALYSIS_DEBUG
					betas(i) *= correction_SlgFailure;
#ifdef FC_CORE_DEBUG
					std::cout << "betas=" << betas << std::endl;
#endif //FC_CORE_DEBUG
				}
			}
		}

		arma::colvec y_norm_center = ms.get_OB_center() / ms.get_OB_delta();
		arma::colvec error_y_tm1 = y_tm1 - y_norm_center;
		arma::colvec scalederror_y_tm1 = betas % error_y_tm1;
		arma::colvec scaled_y_tm1 = scalederror_y_tm1 + y_norm_center;

		if(C.isDefined())
			st.state_transition(x_tm1, scaled_y_tm1);
		else if(C_ZO.isDefined())
			st.state_transition(prev_control_input, scaled_y_tm1);
		else
		{ assert(0); }
	}
	else { // i.e., input exploration was performed in the last frame, or this is the first application frame
		// reset internal metrics, but retain the betas
		for(int i=0; i<(int)vOscillationFailureMode.size(); i++) {
			vOscillationFailureMode[i].reset();
			vSluggishnessFailureMode[i].reset();
		}
	}
}


//wrapper needed for remove_sample()
struct HistoryIndexer {
	const std::deque<FrameObservation>& deqFrameObservation;

	HistoryIndexer(const std::deque<FrameObservation>& deqFrameObservation)
		: deqFrameObservation(deqFrameObservation)
	{}

	const std::vector<double>& operator[](size_t index) const
	{ return deqFrameObservation[index].vCP_observations; }

	size_t size() const
	{ return deqFrameObservation.size(); }
};

void FCdata::resize_history(bool newM, bool doPFE)
{
	size_t min_L_gamma = std::max<size_t>(L_min, std::max<size_t>(ms.active_W, 2 * ms.vActive_CP_IDs.size() + 1));

	//Initially
	if(gamma == 0.0) {
		gamma = std::pow(0.9, 1.0/ms.active_W);
		L_gamma = log(0.1) / log(gamma);
		candidate_L_gamma = L_gamma;
		prev_frame_ended_PFE = false;
	}

	assert(L_gamma > 0);
	if( L_gamma < (int)min_L_gamma || M.isDefined() == false ) {
		L_gamma = std::max( min_L_gamma, (M.isDefined() == false ? deqFrameObservation.size() : 0) );
		candidate_L_gamma = L_gamma;
		gamma = std::pow(0.1, 1.0/L_gamma);
	}

	size_t Ls = 0;
	size_t t_bcp = 0;
	size_t Lc = 0;

	if(frame_number < 0)
		return;

	assert(deqFrameObservation.size() > 0);
	stability_length_metrics.update_stability_length(deqFrameObservation[0],
													 newM,
													 deqFrameObservation.back().frame_number,
													 Ls,   //outputs
													 t_bcp);

	Lc = coverage_length_metrics.update_coverage_length(prev_frame_ended_PFE,
														kappa,
														deqFrameObservation.size(),
														t_bcp,
														frame_number);
	assert((int)Lc >= ms.active_W); //Lc must always be defined

	size_t Lp = 0;
	if(Ls != 0) { //i.e., Ls defined
		if(Ls >= Lc) {
			Lp = Ls;
		}
		else { //Ls < Lc
			HistogramStatistics::stability_type coverageStability = stability_length_metrics.interpolated_stability(Lc);
			switch(coverageStability) {
			case HistogramStatistics::stable:
				Lp = Lc;
				break;
			case HistogramStatistics::unknown:
				stability_length_metrics.add_stability_candidate(Lc);
				Lp = std::max(Ls, candidate_L_gamma);
				break;
			case HistogramStatistics::unstable:
				Lp = (Ls + Lc) / 2;
				break;
			case HistogramStatistics::highlyunstable:
				Lp = Ls;
				break;
			default:
				assert(0);
				break;
			}
		}
	}
	else { //Ls is not defined
		Lp = Lc;
	}

	candidate_L_gamma = (Lp + candidate_L_gamma) / 2;
	if(std::abs<int>(L_gamma - candidate_L_gamma) >= ms.active_W && candidate_L_gamma >= min_L_gamma) {
		L_gamma = candidate_L_gamma;
		gamma = std::pow(0.1, 1.0/L_gamma);
	}

	if(M.isDefined() == true && (deqFrameObservation.size() > (size_t)L_gamma || t_bcp > 0)) {
		size_t L_bcp = deqFrameObservation.size();
		if(t_bcp > 0) {
			assert(frame_number > t_bcp);
			L_bcp = frame_number - t_bcp;
		}
		size_t L_reduced = std::min<size_t>(L_gamma, L_bcp);

#ifdef FC_CORE_DEBUG
		std::cout << "SHORTENING from length=" << (int)deqFrameObservation.size()
			<< " to reduced=" << L_reduced << std::endl;
#endif //FC_CORE_DEBUG

		HistoryIndexer h(deqFrameObservation);

		//First update statistics based on samples to be removed
		size_t length = deqFrameObservation.size();
		for(std::size_t f=(std::size_t)L_reduced; f<deqFrameObservation.size(); f++) {
			remove_sample(h, cs.v_means, cs.v_stds, cs.v_maxswings, cs.v_wmins, cs.v_wmaxs, length, gamma);
			length--;
		}

		//Then remove the samples
		deqFrameObservation.resize(L_reduced);
		deqMSEQ.resize(L_reduced);
	}

#ifdef FC_ANALYSIS_DEBUG
	std::cout << "Ls=" << Ls << " t_bcp=" << t_bcp << " Lc=" << Lc << std::endl;
#endif //FC_ANALYSIS_DEBUG

	prev_frame_ended_PFE = doPFE == true && (forced_length_remaining == 0);
}


double FCdata::budget_allocate()
{
	if(frame_budget <= 0.0) {
		std::cout << "FeatureController::ERROR: a positive frame-budget must be set via set_frame_budget()"
			<< " prior to first invocation of frame_transition()" << std::endl;
		exit(1);
	}

	debt += (prev_time - frame_budget);

	double frame_budget_usable = std::min(frame_budget, frame_budget - debt);

#ifdef FC_CORE_DEBUG
	std::cout << " frame_budget_usable=" << frame_budget_usable << " debt=" << debt << std::endl;
#endif //FC_CORE_DEBUG
	return frame_budget_usable;
}


int generate_random_integer_in_range(int low, int high)
{
	assert(low <= high);
	return ((std::rand() % (int)(high - low + 1)) + low);
}

arma::colvec FCdata::input_explorer()
{
	arma::colvec x_ie = arma::zeros<arma::colvec>( (int)ms.vActive_CP_IDs.size() );

	if((int)deqFrameObservation.size() == 0) {
#ifdef FC_CORE_DEBUG
		std::cout << "input_explorer(): initial x_ie = 0" << std::endl;
#endif //FC_CORE_DEBUG
		return x_ie;
	}

	//Now deqFrameObservation.size() > 0
	arma::colvec CP_N = ms.get_CP_N();
	assert(CP_N.n_rows == x_ie.n_rows);

	x_ie = get_arma_colvec_for_std( deqFrameObservation[0].vCP_observations ) / CP_N;
	assert((int)x_ie.n_rows == (int)cs.num_dims());

#ifdef FC_CORE_DEBUG
	std::cout << "input_explorer(): ";
#endif //FC_CORE_DEBUG

#define UNBIASED_ALL_DIMS_INPUT_EXPLORER
#ifdef UNBIASED_ALL_DIMS_INPUT_EXPLORER

#ifdef FC_CORE_DEBUG
	std::cout << "UNBIASED_ALL: ";
#endif //FC_CORE_DEBUG

	for(int j=0; j<(int)x_ie.n_rows; j++) {
		x_ie(j) = generate_random_integer_in_range(-CP_N[j], +CP_N[j]) / (double)CP_N[j];
#ifdef FC_CORE_DEBUG
		std::cout << "[" << j << "]=" << x_ie(j) << " ";
#endif //FC_CORE_DEBUG
	}

#else //UNBIASED_ALL_DIMS_INPUT_EXPLORER

	if(coverage_met == true) {
		int j = generate_random_integer_in_range(0, x_ie.n_rows-1);
		x_ie(j) = generate_random_integer_in_range(-CP_N[j], +CP_N[j]) / (double)CP_N[j];
#ifdef FC_CORE_DEBUG
		std::cout << "[" << j << "]=" << x_ie(j) << " with coverage_met=1" << std::endl;
#endif //FC_CORE_DEBUG
		return x_ie;
	}

	double required_maxswing_coverage_fraction = get_required_maxswing_coverage_fraction(gamma);

	for(int j=0; j<(int)x_ie.n_rows; j++) {
		//Add 20% margin (i.e., 1.20 multiplier) so that input dimensions just barely
		//having coverage also get explored
		//  ==> reduces chance that a dim will lose coverage just as an exploration cluster ends
		if(cs.v_stds[j]      < 1.20 * required_std_coverage_fraction * CP_N[j] ||
		   cs.v_maxswings[j] < 1.20 * required_maxswing_coverage_fraction * CP_N[j])
		{
			if(std::abs(cs.v_means[j]) < 0.01 * CP_N[j])
				x_ie(j) = generate_random_integer_in_range(-CP_N[j], +CP_N[j]) / (double)CP_N[j];
			else if(cs.v_means[j] > 0.0)
				x_ie(j) = generate_random_integer_in_range(-CP_N[j], lround(cs.v_means[j])) / (double)CP_N[j];
			else
				x_ie(j) = generate_random_integer_in_range(lround(cs.v_means[j]), +CP_N[j]) / (double)CP_N[j];

#ifdef FC_CORE_DEBUG
			std::cout << "[" << j << "]=" << x_ie(j) << " ";
#endif //FC_CORE_DEBUG
		}
	}
#endif //UNBIASED_ALL_DIMS_INPUT_EXPLORER

#ifdef FC_CORE_DEBUG
	std::cout << std::endl;
#endif //FC_CORE_DEBUG

	return x_ie;
}

void FCdata::update_regulatorZO(timeval fc_start_time, double frame_budget_usable)
{
	if( M.isDefined() == true and diff_time(fc_start_time, get_curr_timeval()) < frame_budget_usable
			and (C_ZO.isDefined() == false or C_ZO.is_delta_balanced() == false) )
	{
		if(C_ZO.eta > 0.0 and C_ZO.rho > 0.0)
			C_ZO.delta = std::sqrt(C_ZO.eta * C_ZO.delta / C_ZO.rho);

		const arma::mat& L = M.lp.L2;

		arma::mat Lprime(L.n_rows + L.n_cols, L.n_cols);

		Lprime( arma::span(0, L.n_rows-1)            , arma::span::all ) = L;
		Lprime( arma::span(L.n_rows, Lprime.n_rows-1), arma::span::all ) = C_ZO.delta * arma::eye<arma::mat>(L.n_cols, L.n_cols);

#ifdef FC_ANALYSIS_DEBUG
		std::cout << "---------- LQR SOLUTION (not really, zero-order) ------------" << std::endl;
#endif //FC_ANALYSIS_DEBUG
#ifdef FC_CORE_DEBUG
		std::cout << "L=" << L << std::endl;
		std::cout << "delta=" << C_ZO.delta << std::endl;
#endif //FC_CORE_DEBUG

		arma::mat Lprime_trans = arma::trans(Lprime);
		arma::mat Linv = arma::inv(Lprime_trans * Lprime) * Lprime_trans;
#ifdef FC_CORE_DEBUG
		std::cout << "Linv=" << Linv << std::endl;
#endif //FC_CORE_DEBUG

		C_ZO.reconstruct(L, Linv);

		st = LDS_State(M.lp, ms.get_OB_center()); //reset state
		if(deqFrameObservation.size() > 0) {  //re-initialize state
			arma::colvec y_tm1 = get_arma_colvec_for_std( deqFrameObservation[0].vOB_observations ) / ms.get_OB_delta();
			st.state_transition(prev_control_input, y_tm1);
		}
	}
}



////////////////// DEBUG CONTROL /////////////////////


void set_model_order(FeatureController& fc, int x_order, int y_order)
{
	assert(x_order >= 0);
	assert(y_order >= 0);

	FCdata& fcdata = *(access_vFCdata().at(fc.fcid));

	fcdata.ms.x_order = x_order;
	fcdata.ms.y_order = y_order;

	fcdata.bReinitialize = true;
}

int get_current_model_x_order(const FeatureController& fc)
{
	const FCdata& fcdata = *(access_vFCdata().at(fc.fcid));
	return fcdata.ms.x_order;
}

int get_current_model_y_order(const FeatureController& fc)
{
	const FCdata& fcdata = *(access_vFCdata().at(fc.fcid));
	return fcdata.ms.y_order;
}

} //namespace fctrl


std::ostream& operator<<(std::ostream& os, const fctrl::FCdata& fcdata)
{
	os << std::endl;

	os << fctrl::align() << "FeatureController {" << std::endl;

	fctrl::incr_align();

	os << fcdata.ms;

	for(int f=0; f<(int)fcdata.deqFrameObservation.size(); f++) {
		os << fctrl::align() << "FrameObservation[" << f << "]:" << std::endl;

		fctrl::incr_align();
		os << fcdata.deqFrameObservation[f];
		fctrl::decr_align();
	}
	os << std::endl;

	os << fcdata.M;

	fctrl::decr_align();

	os << fctrl::align() << "}" << std::endl;

	return os;
	
}






