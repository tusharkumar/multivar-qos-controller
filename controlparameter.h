#ifndef __CONTROLPARAMETER_H__
#define __CONTROLPARAMETER_H__

#include <cstdlib>
#include <cassert>
#include <iostream>

#include "fctrl.h"

namespace fctrl {
	class CPdata {
	public:
		bool isDefined;

		//characteristics
		int N;
		std::string description;
		std::vector<ParameterSettingCallback> vCallbacks;
			//Semantics: vCallbacks.size() must be maintained == 2*N+1

		CPdata()
			: isDefined(false)
		{} //undefined, violates semantics. Invalid to perform any operations on such CPdata

		CPdata(int N, std::string description)
			: isDefined(true), N(N), description(description)
		{
			assert(N >= 0);
			vCallbacks.resize(2*N+1, 0);
		}

		CPdata(const CPdata& src)
			: isDefined(true), N(src.N), description(src.description), vCallbacks(src.vCallbacks)
		{}
	};

	std::vector<CPdata>& access_vCPdata();
		//Semantics: vCPdata appends a new element iff a ControlParameter with new cpid is created

} //namespace fctrl

#endif //__CONTROLPARAMETER_H__
