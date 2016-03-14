#include "controlparameter.h"

namespace fctrl {

std::vector<CPdata>& access_vCPdata() {
	static std::vector<CPdata> vCPdata;
	return vCPdata;
}

ControlParameter::ControlParameter()
	: cpid(-1)
{}

ControlParameter::ControlParameter(int N, const std::string& description)
	: cpid(access_vCPdata().size())
{
	access_vCPdata().push_back( CPdata(N, description) );
}

ControlParameter::ControlParameter(const ControlParameter& src)
	: cpid(src.cpid)
{}


ControlParameter& ControlParameter::set_callback(int s, ParameterSettingCallback psc_func)
{
	if(cpid == -1) {
		std::cerr << "ControlParameter::set_callback(): ERROR: invoked on an undefined ControlParameter (cpid=-1)" << std::endl;
		exit(1);
	}

	CPdata& cpdata = access_vCPdata().at(cpid);
	assert(cpdata.isDefined);

	if(s < -cpdata.N or s > cpdata.N) {
		std::cerr << "ControlParameter::set_callback(): ERROR: value s = " << s
				<< " is outside valid range: [" << -cpdata.N << ", " << cpdata.N << "]" << std::endl;
		exit(1);
	}

	int findex = s + cpdata.N;
	cpdata.vCallbacks.at(findex) = psc_func;

	return *this;
}

ControlParameter& ControlParameter::set_allcallback(ParameterSettingCallback psc_func)
{
	if(cpid == -1) {
		std::cerr << "ControlParameter::set_allcallback(): ERROR: invoked on an undefined ControlParameter (cpid=-1)" << std::endl;
		exit(1);
	}

	CPdata& cpdata = access_vCPdata().at(cpid);
	assert(cpdata.isDefined);

	int size = cpdata.vCallbacks.size();
	for(int i=0; i<size; i++)
		cpdata.vCallbacks[i] = psc_func;

	return *this;
}


void ControlParameter::force_characteristics(int new_N, const std::string& new_description)
{
	//tricky: make sure to respect semantics w.r.t. visibility to FeatureController

	if(cpid == -1) {
		std::cerr << "ControlParameter::force_characteristics(): ERROR: invoked on an undefined ControlParameter (cpid=-1)" << std::endl;
		exit(1);
	}

	CPdata& cpdata = access_vCPdata().at(cpid);
	assert(cpdata.isDefined);

	assert((int)cpdata.vCallbacks.size() == 2 * cpdata.N + 1);
		//before

	std::vector<ParameterSettingCallback> orig_vCallbacks = cpdata.vCallbacks;

	cpdata.N           = new_N;
	cpdata.description = new_description;
	cpdata.vCallbacks.clear();
	cpdata.vCallbacks.resize(2*cpdata.N + 1, 0);

	int orig_N = orig_vCallbacks.size();
	for(int i=0; i<std::min(cpdata.N, orig_N); i++) {
		cpdata.vCallbacks[cpdata.N - i] = orig_vCallbacks[orig_N - i];
		cpdata.vCallbacks[cpdata.N + i] = orig_vCallbacks[orig_N + i];
	}
}

void ControlParameter::force_characteristics(const ControlParameter& srcCP)
{
	if(srcCP.cpid == -1) { //copying characteristics of an undefined ControlParameter
		std::cerr << "ControlParameter::force_characteristics(): ERROR: source ControlParameter is undefined (srcCP.cpid=-1)" << std::endl;
		exit(1);
	}

	CPdata& src_cpdata = access_vCPdata().at(srcCP.cpid);
	assert(src_cpdata.isDefined);
	assert((int)src_cpdata.vCallbacks.size() == 2 * src_cpdata.N + 1);

	if(cpid == -1) {
		std::cerr << "ControlParameter::force_characteristics(): ERROR: invoked on an undefined ControlParameter (cpid=-1)" << std::endl;
		exit(1);
	}

	CPdata& cpdata = access_vCPdata().at(cpid);
	assert(cpdata.isDefined);

	assert((int)cpdata.vCallbacks.size() == 2 * cpdata.N + 1);
		//before

	cpdata.N           = src_cpdata.N;
	cpdata.description = src_cpdata.description;
	cpdata.vCallbacks  = src_cpdata.vCallbacks;
}

} //namespace fctrl

