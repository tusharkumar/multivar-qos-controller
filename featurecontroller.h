#ifndef __FEATURE_CONTROLLER_H__
#define __FEATURE_CONTROLLER_H__

#include <cstdlib>
#include <cassert>
#include <iostream>
#include <set>

#include "frameobservation.h"
#include "periodqueue.h"
#include "samplestats.h"
#include "modelstructure.h"
#include "model.h"
#include "linearpredictionmodel.h"
#include "lqr.h"
#include "failuremodes.h"
#include "zeroorder_controller.h"
#include "fctrl.h"
#include "timing.h"
#include "stabilitycoveragemetrics.h"

namespace fctrl {

	class CoverageStats {
	public:
		std::vector<double> v_means;
		std::vector<double> v_stds;
		std::vector<double> v_maxswings;
		std::vector<double> v_wmins;
		std::vector<double> v_wmaxs;
		//Elements of vector correspond to ControlParameters in vCP_IDs (in FCdata).
		//Holds statistics for x samples in deqFrameObservation (in FCdata).

		void reset() {
			v_means.clear();
			v_stds.clear();
			v_maxswings.clear();
			v_wmins.clear();
			v_wmaxs.clear();
		}

		void reinitialize_dims(size_t num_dims) {
			reset();
			v_means.resize(num_dims, 0.0);
			v_stds.resize(num_dims, 0.0);
			v_maxswings.resize(num_dims, 0.0);
			v_wmins.resize(num_dims, 0.0);
			v_wmaxs.resize(num_dims, 0.0);
		}

		size_t num_dims() const { return v_means.size(); }
	};

  class SR_level_metrics {
  public:
    std::multiset<double> achievable_set_SR;
    double achievable_SR;
    double estimated_SR_M;
    double current_SR_M;
    std::deque<double> deq_recent_SR_M;

    SR_level_metrics()
    {
      achievable_set_SR.insert(1.0);
      achievable_SR  = 1.0;
      estimated_SR_M = 1.0;
      current_SR_M   = 1.0;
      deq_recent_SR_M.push_front(1.0);
    }

    void reset_for_new_M()
    {
      estimated_SR_M = 1.0;
      current_SR_M   = 1.0;
      deq_recent_SR_M.clear();
      deq_recent_SR_M.push_front(1.0);
    }

    void update_achievability()
    {
      achievable_set_SR.insert(current_SR_M);

      achievable_SR = 0.0;
      std::size_t quarter_size_SR = (achievable_set_SR.size() + 3) / 4;
        //average over the *largest* quarter of elements
      std::multiset<double>::iterator it_SR=achievable_set_SR.begin();
      std::advance(it_SR, achievable_set_SR.size()-quarter_size_SR);
      for( ; it_SR != achievable_set_SR.end(); it_SR++)
      { achievable_SR += *it_SR; }
      achievable_SR = achievable_SR / quarter_size_SR;

#ifdef FC_CORE_DEBUG
      std::cout << "Sampled: current_SR_M=" << current_SR_M
          << " producing: achievable_SR=" << achievable_SR << std::endl;
      std::cout << "achievable_set_SR=" << achievable_set_SR << std::endl;
#endif //FC_CORE_DEBUG
    }

    void add_SR_sample(size_t num_frames_of_M, bool frame_satisfies_level, size_t max_q_length)
    {
			estimated_SR_M  = (estimated_SR_M * num_frames_of_M + (frame_satisfies_level ? 1.0 : 0.0))
									/ (num_frames_of_M + 1);

			deq_recent_SR_M.push_front( (frame_satisfies_level ? 1.0 : 0.0) );
			if(deq_recent_SR_M.size() > max_q_length)
				deq_recent_SR_M.resize(max_q_length);

			current_SR_M = 0.0;
			for(int i=0; i<(int)deq_recent_SR_M.size(); i++)
				current_SR_M += deq_recent_SR_M[i];
			current_SR_M = current_SR_M / deq_recent_SR_M.size();
    }
  };

	class FCdata {
	public:
		long long int frame_number;

		int W;
		int requested_x_order;
		int requested_y_order;
		int L_min;

		std::vector<CP_ID> vCP_IDs;
			//associated ControlParameters, identified by their cpids.
			//Semantics: Maintained in ascending order of cpids

		std::vector<OB_ID> vOB_IDs;
			//associated Objectives, identified by their obids
			//Semantics: Maintained in ascending order of obids

		bool bReinitialize;
			//if true, will force a reinitialization at next frame-transition

		ModelStructure ms;


		std::deque<FrameObservation> deqFrameObservation;
			//The available history, ordered from most recent to oldest.
			//Each FrameObservation must have observations for ControlParameters
			//  and Objectives corresponding to the elements of ms.vActive_CP_IDs and ms.vActive_OB_IDs

		CoverageStats cs;

		std::deque<double> deqMSEQ;
			//Maintained of same length as deqFrameObservation
			//Holds the per-frame MSEQ for the corresponding frames


			//kappa and coverage_met are set at start of frame
		double kappa;
			//captures the fraction of input dimensions that sufficiently span their space
		bool   coverage_met;
			//indicates sufficient coverage in history to allow model estimation
		Model M;
		long long int M_est_framenumber;
			//frame-number during which M was estimated. Defined iff M.isDefined() == true
		bool prev_newM;
			//indicates if a new M was generated in the previous frame by update_model()
		StabilityLengthMetrics stability_length_metrics;
			//stability length metrics seen over the active M's 
		CoverageLengthMetrics coverage_length_metrics;
			//metrics that estimate the length necessary to reach coverage in the observed history
		Model Mp;
			//copy used to separately refine M, and then apply changes back to M.
			//Invariant: M.isDefined() iff Mp.isDefined()
		long long int Mp_est_framenumber;
			//frame-number during which Mp was estimated. Defined iff Mp.isDefined() == true
		double lambda_next;
			//the value of lambda to use for next LLSE (persistent across resets)
		double advantage_Mp_over_M;
			//the cumulative prediction accuracy advantage of Mp over M
			//  based on history data (defined iff M.isDefined() == true)
		CoverageStats subst_cs;
		long long int subst_history_count;
			//history stats over inputs seen since last M or Mp was applied.
			//computed over same range of frames as advantage_Mp_over_M (defined iff Mp.isDefined() == true)
		bool sameModel;
			//indicates if M and Mp are the same model (defined iff M.isDefined() == true)

		bool prev_apply_controller;
			//Indicates if controller was applied last frame.

		LDS_State st;
			//Linear Dynamical System state at current step, dimensionality defined by model-structure of M
			// (defined iff C.isDefined() == true)

		arma::colvec diagR;
			//Diagonal terms of R, the diagonal matrix specifying the LQR input costs

		arma::colvec c;
			//counters for updating R

		std::vector<bool> f_barely_underconstrained;
			//indicates whether j'th input was barely underconstrained in previous time-step
			//  and in previous step of input-cost refinement


		LQR_Controller C;
		long long int C_refine_timestep;
		int LQR_timestep;

		ZeroOrder_Controller C_ZO;

		arma::colvec prev_control_input;
			//normalized representation


		std::vector<OscillationFailureMode> vOscillationFailureMode;
		std::vector<SluggishnessFailureMode> vSluggishnessFailureMode;
			//maintained for each objective in vOB_IDs

		arma::colvec last_correction;
			//last correction applied to beta of the corresponding output

		arma::colvec betas;
			//scale-factors corresponding to each objective y_i.
			//ith term beta_i magnifies/diminishes the observed error in y_i
			//   i.e. observed_error_i = y_i - objective_i
			//     gets adjusted: scaled_error_i = beta_i * observed_error_i
			//
			//  And, scaled_y_i = scaled_error_i + objective_i
			//     is used to update state


			//Metrics and State for Forced Exploration

		long long int last_exp_frame_number;
			//last frame on which C did not drive inputs

    static const double SR_level = 1.50;

    std::vector<SR_level_metrics> v_SR_level_metrics;
      // Levels l = 0, 1, 2, 3, ... index the corresponding SR metrics for each level.

    size_t l_tr;
      // Identifies the current objective being tracked:
      //
      //   Objective: tau <= k_tr
      //
      //   where k_tr = (1.50)**l_tr

    double max_frame_MSEQ;
      //The maximum MSEQ encountered on any frame so far (sliding-window averaged).
      // Kept >= 1.0

		long long int num_frames_of_M;
			//number of frames on which the active model M has been driven inputs so far

		double theta;

		int forced_length_remaining;
			//forced length remaining in current PFE cluster (if > 0)

		size_t num_PFE_clusters;
      //current setting of number of PFE clusters in a PFE group

    size_t cluster_index_in_PFE_group;
      //index within the PFE group --- of the ongoing or the last PFE cluster completed
      //index = 1 to num_PFE_clusters
      //index = 0 ==> the next PFE cluster will start a new PFE group

		std::vector<double> d_cdf;


		double gamma;
			//history forget rate
		int L_gamma;
			//the length of the significant history dictated by gamma
		size_t candidate_L_gamma;
		bool prev_frame_ended_PFE;
			//flag to signal if previous frame was last frame of a PFE


		double frame_budget;
			//total time allocated per-frame for featurecontroller to execute optimizations
		double prev_time;
			//total time consumed by FeatureController in previous frame (persistent across resets)
		double debt;
			//time consumed by FeatureController in excess of frame_budget, accumulated over all
			//prior frames (slack to frame_budget subtracts from debt, exceeding frame_budget adds to debt)


		AllStatistics stats;

		bool wereFixedSettingsAppliedPreviousFrame;
		std::vector< std::pair<CP_ID, int> > vCP_fixedSetting_pairs;
			//vCP_fixedSetting_pairs == undefined ==> normal operation via run_controller()
			//
			//vCP_fixedSetting_pairs == defined   ==> fixed control-inputs via run_fixed_settings().
			//  A fixed setting must be provided for every ControlParameter associated with this FeatureController


    std::vector< std::pair<CP_ID, int> > vCP_overrideSettings_pairs;
      //Overrides the value applied by the controller in the next frame.
      //Cleared at the end of each frame.

    bool override_settings_requested_for_frame;
      //Indicates if control input settings override was requested in the previous frame
      // for this frame.

			//Zero-Order model-structure functionality

		//////

		void reinitialize_setup();
			//Apply the latest settings requested by the application

		void reset_beta_failuremodes();

		std::vector<double> read_objectives();
			//Store raw measurements into *this* frame's FrameObservation.
			//Objective values are stored in un-normalized form.
			//Return filtered measurements.

		void write_controlparameters(const arma::colvec& curr_control_input);
			//Store a normalized and clipped version of curr_control_input into *next* frame's FrameObservation.
			//Elements of curr_control_input must correspond to elements of ms.vActive_CP_IDs


		////// Operation alternatives

		arma::colvec run_fixed_settings();
			//Apply fixed settings (specified by user) for the control-parameters in every frame
			//Returns application control inputs for the current frame

		arma::colvec run_controller(timeval fc_start_time);
			//Use observed history to estimate model and determine control inputs to apply
			//  to meet QoS objectives.
			//Returns application control inputs for the current frame



		////// Supporting functionality for run_controller()

		bool update_model(timeval fc_start_time, double frame_budget_usable, bool forcedExploration);
			//returns whether a new model was applied

		bool update_regulator(timeval fc_start_time, double frame_budget_usable);
			//returns whether a new controller was constructed

		void check_history_coverage(
				const std::vector<double>& v_history_stds,
				const std::vector<double>& v_history_maxswings,
				std::size_t                history_len,
				double                     history_gamma,
				double&                    history_kappa,       //output
				bool&                      history_coverage_met //output
		);
			//does history explore a sufficient part of input-space to allow
			// robust model estimation?
			//outputs history_kappa and history_coverage_met

		Model llse_and_refine_lambda();
			//perform LLSE on current history in deqFrameObservation
			//   using current lambda_next, and refine lambda_next for next time

		bool force_model_estimation(long long int last_model_estimation_timestamp);
			//probabilistically determines if model-estimation should be forced this frame,
			//  based on how long ago the previous model-estimation was performed

		double initial_refinement(
			double f_barely_underconstrained_j,
			double abs_x_prev_j,
			double diagR_prev_j,
			double c_j
		);
			//returns updated diagR_j

		void refine_input_cost(
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
		);

		void evaluate_terminated_refinement(
				//Outputs
			double& diagR_j,
			bool  & f_term_j, //also an input
				//Inputs
			double  abs_x_prev_j,
			double  abs_x_refined_j,
			double  diagR_prev_j,
			double  diagR_refined_j,
			double  c_j
		);

		double increase_input_cost(double diagR_j_before, double c_j);
		double decrease_input_cost(double diagR_j_before, double c_j);

		bool practically_unchanged(double abs_x_prev_j, double abs_x_refined_j);


		bool probabilistic_forced_exploration();
			//whether input-exploration should be forced in support of maintaining sufficient
			//  history for model estimation, despite presence of controller.
			//updates forced_length_remaining to reflect number of frames left in current PFE
			//  cluster if one is in progress; otherwise = 0.
			//retains num_PFE_clusters as best guess for next frame.

		void update_achievability();
      //updates v_SR_level_metrics

		size_t estimate_shortest_PFE_cluster_length();

		void adjust_cluster_length_distribution_and_sample(
				const int d_peak,
				const int d_high,
				double& expected_d,
				int& sampled_d
		);

		void adaptive_state_transition();

		void resize_history(bool newM, bool doPFE);
			//reduces or increases length of history retained in deqFrameObservation
			//  based on estimated length suitable for model-estimation

		double budget_allocate();
			//determines how much time is available for non-critical
			//model estimation (of Mp) and refinement of C

		arma::colvec input_explorer();
			
		void update_regulatorZO(timeval fc_start_time, double frame_budget_usable);

		FCdata()
			: frame_number(-1), W(1), requested_x_order(0), requested_y_order(1), L_min(0),
				bReinitialize(false), ms(), coverage_met(false),
				M(&ms), prev_newM(false), stability_length_metrics(&M), coverage_length_metrics(&M),
				Mp(0), lambda_next(0.000001 * 0.000001),
				advantage_Mp_over_M(0), sameModel(false), prev_apply_controller(false),
				C(), LQR_timestep(0),
				C_ZO(&ms),
				last_exp_frame_number(-1),
				v_SR_level_metrics(), l_tr(0), max_frame_MSEQ(1.0), 
        num_frames_of_M(0),
			  theta(), forced_length_remaining(), num_PFE_clusters(1), cluster_index_in_PFE_group(0), d_cdf(),
				gamma(0.0), L_gamma(-1),
				candidate_L_gamma(0), prev_frame_ended_PFE(false),
				frame_budget(-1), prev_time(0.0), debt(0.0),
				wereFixedSettingsAppliedPreviousFrame(false),
        override_settings_requested_for_frame(false)
        
		{}
	};

	
	std::vector<FCdata *>& access_vFCdata();
		//Semantics: vFCdata appends a new element iff a FeatureController with new fcid is created
} //namespace fctrl

std::ostream& operator<<(std::ostream& os, const fctrl::FCdata& fcdata);

template<typename OStream>
OStream& operator<<(OStream& ostream, const fctrl::CoverageStats& cs)
{
	for(std::size_t j=0; j<cs.num_dims(); j++) {
		ostream << "[" << j << "] mean=" << cs.v_means[j] << " std=" << cs.v_stds[j]
				<< " maxswing=" << cs.v_maxswings[j] << " wmins=" << cs.v_wmins[j] << " wmaxs=" << cs.v_wmaxs[j]
				<< std::endl;
	}
	ostream << std::endl;
	return ostream;
}

#endif //__FEATURE_CONTROLLER_H__
