#ifndef __OBJECTIVE_H__
#define __OBJECTIVE_H__

#include <cstdlib>
#include <cassert>
#include <iostream>

#include "utils.h"
#include "fctrl.h"

namespace fctrl {


	void normalize_benefits(std::vector<MOTerm>& v_moterms);

	class OBJdata {
	public:
		bool isDefined;

		Objective::Type type;
		std::string description;
		double importance;
		MeasurementCallback callback;

		//Range-objective: defined iff type == Objective::eRange
		double center;
		double width;

		//Maximization-objective: defined iff type == Objective::eMaximization
		std::vector<MOTerm> v_moterms;

		CircularQueue<double> samples;
			//recorded data samples for this Objective,
			//circular queue of length sliding_window_length

		double raw_sample;

		OBJdata()
			: isDefined(false)
		{} //undefined, violates semantics. Invalid to perform any operation on such OBJdata

		OBJdata(
			double center,
			double width,
			int sliding_window_length,
			const std::string& description = "",
			double importance = 1.0
		) : isDefined(true), type(Objective::eRange), description(description),
				importance(importance), callback(0),
				center(center), width(width), samples(sliding_window_length), raw_sample(0)
		{
			assert(sliding_window_length >= 1);
			assert(importance >= 0.0);
			assert(width >= 0.0);
		}

		OBJdata(
			const std::vector<MOTerm>& v_moterms,
			int sliding_window_length,
			const std::string& description = "",
			double importance = 1.0
		) : isDefined(true), type(Objective::eMaximization), description(description),
				importance(importance), callback(0),
				v_moterms(v_moterms), samples(sliding_window_length)
		{
			assert(sliding_window_length >= 1);
			assert(importance >= 0.0);
			assert(v_moterms.size() >= 2);
			for(int i=1; i<(int)v_moterms.size(); i++)
				assert(v_moterms[i-1].center <= v_moterms[i].center);

			normalize_benefits(this->v_moterms);
		}

		void perform_measurement() {
			if(callback == 0) {
				std::cerr << "Objective::ERROR: measurement-callback not provided for reading value:"
					<< " Objective description = " << description << std::endl;
				exit(1);
			}
			double obj_value = callback(); //invoke callback
			samples.push(obj_value);
			raw_sample = obj_value;
		}

		double read_sliding_window_averaged_measurement() const
		{ return samples.read_average(); }

		double read_raw_measurement() const
		{ return raw_sample; }

	};

	std::vector<OBJdata>& access_vOBJdata();
		//Semantics: vOBJdata appends a new element iff an Objective with new obid is created


} //namespace fctrl

#endif //__OBJECTIVE_H__
