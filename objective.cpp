#include "objective.h"

namespace fctrl {

std::vector<OBJdata>& access_vOBJdata() {
	static std::vector<OBJdata> vOBJdata;
	return vOBJdata;
}


void normalize_benefits(std::vector<MOTerm>& v_moterms) {
	double total_benefit = 0.0;
	for(int i=0; i<(int)v_moterms.size(); i++) {
		assert(v_moterms[i].incr_benefit >= 0);
		total_benefit += v_moterms[i].incr_benefit;
	}
	assert(total_benefit > 0.0);

	for(int i=0; i<(int)v_moterms.size(); i++)
		v_moterms[i].incr_benefit = v_moterms[i].incr_benefit / total_benefit;
}



Objective::Objective()
	: obid(-1)
{}


Objective::Objective(
	double center,
	double width,
	int sliding_window_length,
	const std::string& description,
	double importance
) : obid(access_vOBJdata().size())
{
	access_vOBJdata().push_back( OBJdata(center, width, sliding_window_length, description, importance) );
}

Objective::Objective(
	const std::vector<MOTerm>& v_moterms,
	int sliding_window_length,
	const std::string& description,
	double importance
) : obid(access_vOBJdata().size())
{
	access_vOBJdata().push_back( OBJdata(v_moterms, sliding_window_length, description, importance) );
}


Objective::Objective(const Objective& src)
	: obid(src.obid)
{}


Objective::Type Objective::get_type() const
{
	if(obid == -1) {
		std::cerr << "Objective::get_type(): ERROR: Invoked for obid=-1" << std::endl;
		exit(1);
	}

	const OBJdata& objdata = access_vOBJdata().at(obid);
	assert(objdata.isDefined);

	return objdata.type;
}

void Objective::set_measurement_callback(MeasurementCallback msc_func)
{
	if(obid == -1) {
		std::cerr << "Objective::set_measurement_callback(): ERROR: Invoked for obid=-1" << std::endl;
		exit(1);
	}

	OBJdata& objdata = access_vOBJdata().at(obid);
	assert(objdata.isDefined);
	objdata.callback = msc_func;
}


void Objective::force_characteristics(const Objective& srcObj)
{
	if(obid == -1) {
		std::cerr << "Objective::force_characteristics(): ERROR: Invoked for obid=-1" << std::endl;
		exit(1);
	}

	OBJdata& objdata = access_vOBJdata().at(obid);
	//assert(objdata.isDefined);


	if(srcObj.obid == -1) {
		std::cerr << "Objective::force_characteristics(): ERROR: Invoked for srcObj.obid=-1" << std::endl;
		exit(1);
	}

	OBJdata& src_objdata = access_vOBJdata().at(srcObj.obid);
	//assert(src_objdata.isDefined);
	
	objdata = src_objdata;
}

} //namespace fctrl
