#ifndef __MODEL_STRUCTure_H__
#define __MODEL_STRUCTure_H__

#include <iostream>
#include <deque>

#include <armadillo>

#include "frameobservation.h"
#include "fctrl.h"

namespace fctrl {

	std::vector<double> get_std_for_arma_colvec(const arma::colvec& arma_colvec);
	arma::colvec get_arma_colvec_for_std(const std::vector<double>& vec);

	arma::colvec compute_SEQ_objectives(
		const arma::colvec& Y_measured,
		const arma::colvec& Y_centers,
		const arma::colvec& Y_deltas,
		const arma::colvec& Y_importances
	);

	double compute_MSEQ(
		const arma::colvec& Y_measured,
		const arma::colvec& Y_centers,
		const arma::colvec& Y_deltas,
		const arma::colvec& Y_importances
	);

	class ModelStructure {
	public:

		int active_W;
			//perception-window length

		std::vector<CP_ID> vActive_CP_IDs;
			//The associated ControlParameters
		std::vector<OB_ID> vActive_OB_IDs;
			//The associated Objectives

		//Semantics: maintained in sorted order of cpids and obids.



		//////// Structuring Parameters ////////

		int x_order;
			//x_order == 0 => no dependence on past inputs, only current inputs;
			//x_order >= 1 indicates past dependence
		int y_order;
			//y_order == 0 => no dependence on past outputs, hence no feedback-control
			//y_order >=1 for feedback control

			//TODO: generalize idea of order to range, and then to ControlParameter-specific range
		
		bool flag_subtract_ytm1_in_estimation_for_DC_offset_correction;

		bool flag_subtract_yoperating_in_estimation_for_DC_offset_correction;

		bool flag_subtract_yobjective_in_estimation_for_DC_offset_correction;

		bool flag_use_affine_model;

		bool flag_use_adaptive_betas;

		//////// Behavior Observations ////////

		arma::colvec Y_operating;
			//The current operating-points for the corresponding Objectives in vActive_OB_IDs.
			//Only entries corresponding to eMaximization Objectives need to be valid.


		////////////////////////////////////////


		ModelStructure()
			: active_W(1), x_order(0), y_order(1),
				flag_subtract_ytm1_in_estimation_for_DC_offset_correction(true),
				flag_subtract_yoperating_in_estimation_for_DC_offset_correction(false),
				flag_subtract_yobjective_in_estimation_for_DC_offset_correction(false),
				flag_use_affine_model(false),
				flag_use_adaptive_betas(false)
		{}


		//////// Decision Making routines ////////

		void estimate_operating_point_from_data_samples(
				const std::deque<FrameObservation>& deqSamples,
				double gamma
		);
			//Sets Y_operating


		//////// Decision Access routines ////////

		arma::colvec get_CP_N();
			//These construct the information from the ControlParameter in vActive_CP_IDs,
			//  reusing any previously computed information.

		arma::colvec get_OB_center();
		arma::colvec get_OB_delta();
		arma::colvec get_OB_importance();
			//These construct the information from the Objectives in vActive_OB_IDs,
			//  reusing any previously computed information.

		//////// Conversion routines ////////

		arma::colvec     convert_CP_choices_to_normalized(const std::vector<int>& vCP_choices);
		std::vector<int> convert_normalized_to_clipped_CP_choices(const arma::colvec& norm_CP_choices);

		arma::colvec   convert_regular_Y_to_normalized(const arma::colvec& regular_Y);
		arma::colvec   convert_normalized_Y_to_regular(const arma::colvec& normalized_Y);

	private:
		arma::colvec X_Ns;
			//The range-sizes, N, for the ControlParameters in vActive_CP_IDs.

		arma::colvec Y_centers;
			//The currently chosen centers corresponding to Objectives in vActive_OB_IDs.
			//  The eMaximization objectives have centers that are not pre-determined.

		arma::colvec Y_deltas;
			//The currently chosen deltas corresponding to Objectives in vActive_OB_IDs.
			//  The eMaximization objectives have deltas that are not pre-determined.

		arma::colvec Y_importances;
			//The currently chosen importances corresponding to Objectives in vActive_OB_IDs.
			//  The eMaximization objectives have importances that are not pre-determined.

		arma::colvec Y_control_importances;
			//Same as Y_importances, expect allowing for further tweaking by the
			//  control-optimization strategy.
	};

} //namespace fctrl

std::ostream& operator<<(std::ostream& os, const fctrl::ModelStructure& ms);

#endif //__MODEL_STRUCTure_H__
