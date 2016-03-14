#include "utils.h"
#include "controlparameter.h"
#include "objective.h"
#include "modelstructure.h"

namespace fctrl {


std::vector<double> get_std_for_arma_colvec(const arma::colvec& arma_colvec)
{
	std::vector<double> std_vec(arma_colvec.n_rows);
	for(int i=0; i<(int)arma_colvec.n_rows; i++)
		std_vec[i] = arma_colvec(i);

	return std_vec;
}

arma::colvec get_arma_colvec_for_std(const std::vector<double>& vec)
{
	arma::colvec acv(vec.size());
	for(int i=0; i<(int)vec.size(); i++)
		acv(i) = vec[i];

	return acv;
}


arma::colvec compute_SEQ_objectives(
	const arma::colvec& Y_measured,
	const arma::colvec& Y_centers,
	const arma::colvec& Y_deltas,
	const arma::colvec& Y_importances
)
{
	arma::colvec normalized_error = (Y_measured - Y_centers) / Y_deltas;
	arma::colvec sne = arma::square(normalized_error);
	
	return sne;
}

double compute_MSEQ(
	const arma::colvec& Y_measured,
	const arma::colvec& Y_centers,
	const arma::colvec& Y_deltas,
	const arma::colvec& Y_importances
)
{
	arma::colvec normalized_error = (Y_measured - Y_centers) / Y_deltas;
	arma::colvec sne = arma::square(normalized_error);
	arma::colvec weighted_sne = sne % Y_importances;
	double mseq = arma::sum(weighted_sne) / weighted_sne.n_rows;
	
	return mseq;
}


////////////////// class ModelStructure ///////////////////


#if 0
void ModelStructure::estimate_operating_point_from_data_samples(
		const std::deque<FrameObservation>& deqSamples
)
{
	Y_operating = arma::zeros<arma::colvec>(vActive_OB_IDs.size());

	if((int)deqSamples.size() == 0) { //no data samples available
		//Construct from Objectives instead
		for(int i=0; i<(int)vActive_OB_IDs.size(); i++) {
			OB_ID obid = vActive_OB_IDs[i];
			const OBJdata& objdata = access_vOBJdata().at(obid);
			assert(objdata.isDefined);

			if(objdata.type == Objective::eRange) {
				Y_operating(i) = objdata.center;
			}
			else if(objdata.type == Objective::eMaximization) {
				assert(objdata.v_moterms.size() >= 2);
				double delta_to_next = objdata.v_moterms[1].center - objdata.v_moterms[0].center;
				assert(delta_to_next > 0);
				Y_operating(i) = objdata.v_moterms[0].center - delta_to_next;
					//set below most conservative maximization-objective step
			}
		}
	}
	else { //data sample available
		for(int f=0; f<(int)deqSamples.size(); f++) {
			const FrameObservation& fo = deqSamples[f];

			assert(fo.vOB_observations.size() == vActive_OB_IDs.size());
			for(int i=0; i<(int)fo.vOB_observations.size(); i++) {
				Y_operating(i) += fo.vOB_observations[i];
			}
		}
		Y_operating = (1.0 / (double)deqSamples.size()) * Y_operating;
			//average over data samples
	}
}

#else
void ModelStructure::estimate_operating_point_from_data_samples(
		const std::deque<FrameObservation>& deqSamples,
		double gamma
)
{
	Y_operating = arma::zeros<arma::colvec>(vActive_OB_IDs.size());

	if((int)deqSamples.size() == 0) { //no data samples available
		//Construct from Objectives instead
		for(int i=0; i<(int)vActive_OB_IDs.size(); i++) {
			OB_ID obid = vActive_OB_IDs[i];
			const OBJdata& objdata = access_vOBJdata().at(obid);
			assert(objdata.isDefined);

			if(objdata.type == Objective::eRange) {
				Y_operating(i) = objdata.center;
			}
			else if(objdata.type == Objective::eMaximization) {
				assert(objdata.v_moterms.size() >= 2);
				double delta_to_next = objdata.v_moterms[1].center - objdata.v_moterms[0].center;
				assert(delta_to_next > 0);
				Y_operating(i) = objdata.v_moterms[0].center - delta_to_next;
					//set below most conservative maximization-objective step
			}
		}
	}
	else { //data sample available
		assert(0 < gamma and gamma < 1);
		double multiplier = 1.0;
		for(int f=0; f<(int)deqSamples.size(); f++) {
			const FrameObservation& fo = deqSamples[f];

			assert(fo.vOB_observations.size() == vActive_OB_IDs.size());
			for(int i=0; i<(int)fo.vOB_observations.size(); i++) {
				Y_operating(i) += fo.vOB_observations[i] * multiplier;
			}
			multiplier *= gamma;
		}
		Y_operating = (1 - gamma) * (1.0 / (double)deqSamples.size()) * Y_operating;
			//average over data samples
	}
}
#endif


arma::colvec ModelStructure::get_CP_N()
{
	if(X_Ns.n_rows > 0) { //re-use pre-computed
		assert(X_Ns.n_rows == vActive_CP_IDs.size());
		return X_Ns;
	}

	//Compute X_Ns
	X_Ns.set_size(vActive_CP_IDs.size());

	for(int j=0; j<(int)vActive_CP_IDs.size(); j++) {
		CP_ID cpid = vActive_CP_IDs[j];
		const CPdata& cpdata = access_vCPdata().at(cpid);
		assert(cpdata.isDefined);

		X_Ns(j) = cpdata.N;
	}

	return X_Ns;
}

arma::colvec ModelStructure::get_OB_center()
{
	if(Y_centers.n_rows > 0) {
		assert(Y_centers.n_rows == vActive_OB_IDs.size());
		assert(Y_deltas.n_rows == vActive_OB_IDs.size());
		assert(Y_importances.n_rows == vActive_OB_IDs.size());
		return Y_centers;
	}

	assert(Y_operating.n_rows == vActive_OB_IDs.size());
	Y_centers.set_size(vActive_OB_IDs.size());
	Y_deltas.set_size(vActive_OB_IDs.size());
	Y_importances.set_size(vActive_OB_IDs.size());

	for(int i=0; i<(int)vActive_OB_IDs.size(); i++) {
		OB_ID obid = vActive_OB_IDs[i];
		const OBJdata& objdata = access_vOBJdata().at(obid);
		assert(objdata.isDefined);

		if(objdata.type == Objective::eRange) {
			Y_centers(i) = objdata.center;
			Y_deltas(i) = objdata.width;
			Y_importances(i) = objdata.importance;
		}
		else if(objdata.type == Objective::eMaximization) {
			assert(objdata.v_moterms.size() >= 2);

			bool found = false;
			double center_benefit=1.0;
			for(int l=0; l<(int)objdata.v_moterms.size(); l++) {
				if(Y_operating(i) < objdata.v_moterms[l].center) {
					Y_centers(i) = objdata.v_moterms[l].center;
					if(l==0)
					{ Y_deltas(i) = objdata.v_moterms[l+1].center - objdata.v_moterms[l].center; }
					else
					{ Y_deltas(i) = objdata.v_moterms[l].center - objdata.v_moterms[l-1].center; }
					Y_importances(i) = objdata.importance * center_benefit;

					found = true;
					break;
				}
				center_benefit -= objdata.v_moterms[l].incr_benefit;
				assert(center_benefit > 0.0);
			}
			if(found == false) { //i.e., larger than all the centers
				int last_index = objdata.v_moterms.size()-1;
				Y_centers(i) = objdata.v_moterms.at(last_index).center;
				Y_deltas(i) = objdata.v_moterms[last_index].center - objdata.v_moterms.at(last_index-1).center;
				assert(center_benefit > 0.0); //FIXME: current fctrl.h spec will force 0.0 here, instead of residual benefit
				Y_importances(i) = objdata.importance * center_benefit;
			}
		}
		else
		{ assert(0); }
	}

	return Y_centers;
}

arma::colvec ModelStructure::get_OB_delta()
{
	if(Y_deltas.n_rows == 0)
		get_OB_center(); //force re-computation

	assert(Y_centers.n_rows == vActive_OB_IDs.size());
	assert(Y_deltas.n_rows == vActive_OB_IDs.size());
	assert(Y_importances.n_rows == vActive_OB_IDs.size());

	return Y_deltas;
}

arma::colvec ModelStructure::get_OB_importance()
{
	if(Y_importances.n_rows == 0)
		get_OB_center(); //force re-computation

	assert(Y_centers.n_rows == vActive_OB_IDs.size());
	assert(Y_deltas.n_rows == vActive_OB_IDs.size());
	assert(Y_importances.n_rows == vActive_OB_IDs.size());

	return Y_importances;
}


arma::colvec ModelStructure::convert_CP_choices_to_normalized(const std::vector<int>& vCP_choices)
{
	assert(vCP_choices.size() == vActive_CP_IDs.size());
	arma::colvec norm_CP_choices(vActive_CP_IDs.size());

	for(int j=0; j<(int)vActive_CP_IDs.size(); j++) {
		CP_ID cpid = vActive_CP_IDs[j];
		const CPdata& cpdata = access_vCPdata().at(cpid);
		assert(cpdata.isDefined);

		assert(cpdata.N > 0);

		norm_CP_choices(j) = ((double)vCP_choices.at(j)) / cpdata.N;
	}

	return norm_CP_choices;
}


std::vector<int> ModelStructure::convert_normalized_to_clipped_CP_choices(const arma::colvec& norm_CP_choices)
{
	assert(norm_CP_choices.n_rows == vActive_CP_IDs.size());
	std::vector<int> vCP_choices(vActive_CP_IDs.size());

	for(int j=0; j<(int)vActive_CP_IDs.size(); j++) {
		CP_ID cpid = vActive_CP_IDs[j];
		const CPdata& cpdata = access_vCPdata().at(cpid);
		assert(cpdata.isDefined);

		int sign = 1;
		if(norm_CP_choices(j) < 0)
			sign = -1;

		int choice = int(norm_CP_choices(j) * cpdata.N + 0.5 * sign);
		if(choice > cpdata.N)
		{ vCP_choices[j] = cpdata.N; }
		else if(choice < -cpdata.N)
		{ vCP_choices[j] = -cpdata.N; }
		else
		{ vCP_choices[j] = choice; }
	}

	return vCP_choices;
}


arma::colvec ModelStructure::convert_regular_Y_to_normalized(const arma::colvec& regular_Y)
{ return regular_Y / get_OB_delta(); }

arma::colvec ModelStructure::convert_normalized_Y_to_regular(const arma::colvec& normalized_Y)
{ return normalized_Y % get_OB_delta(); }


} //namespace fctrl


std::ostream& operator<<(std::ostream& os, const fctrl::ModelStructure& ms)
{
	os << fctrl::align() << "active_W=" << ms.active_W << std::endl;
	os << fctrl::align() << "vActive_CP_IDs= " << ms.vActive_CP_IDs << std::endl;
	os << fctrl::align() << "vActive_OB_IDs= " << ms.vActive_OB_IDs << std::endl;
	os << fctrl::align() << "x_order=" << ms.x_order << std::endl;
	os << fctrl::align() << "y_order=" << ms.y_order << std::endl;
	os << fctrl::align() << "flag_subtract_ytm1_in_estimation_for_DC_offset_correction="
			<< (ms.flag_subtract_ytm1_in_estimation_for_DC_offset_correction ? "true" : "false") << std::endl;
	os << fctrl::align() << "flag_subtract_yoperating_in_estimation_for_DC_offset_correction="
			<< (ms.flag_subtract_yoperating_in_estimation_for_DC_offset_correction ? "true" : "false") << std::endl;
	os << fctrl::align() << "flag_subtract_yobjective_in_estimation_for_DC_offset_correction="
			<< (ms.flag_subtract_yobjective_in_estimation_for_DC_offset_correction ? "true" : "false") << std::endl;
	os << fctrl::align() << "flag_use_affine_model=" << (ms.flag_use_affine_model ? "true" : "false") << std::endl;
	os << fctrl::align() << "flag_use_adaptive_betas=" << (ms.flag_use_adaptive_betas ? "true" : "false") << std::endl;

	return os;
}

