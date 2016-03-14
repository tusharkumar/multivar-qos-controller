#ifndef __FCTRL_H__
#define __FCTRL_H__

#include <string>
#include <vector>
#include <map>

namespace fctrl {

	////////////////////////////////////////////////
	//////// ControlParameter definitions 
	////////////////////////////////////////////////

	typedef long long int CP_ID;
		//Unique ID associated with each created ControlParameter

	typedef void (* ParameterSettingCallback)(int s);
		//Controller will invoke this function to apply a control setting for
		//   the associated ControlParameter
		//s: the setting being applied is passed as the argument.


	class ControlParameter {
	public:
		//Semantics: A ControlParameter gets created with an immutable set of characteristics

		const CP_ID cpid;
			//Unique ID associated with each ControlParameter (unique independent of number of controllers created)
			// cpid = -1 for an undefined ControlParameter


		ControlParameter();
			//Undefined control-parameter, will not be used for control
			// cpid = -1

		ControlParameter(
			int N,
				//range-width of parameter: valid values of parameter are between -N to N (both inclusive)
			const std::string& description = ""
				//optional description, not used by controller
		);

		ControlParameter(const ControlParameter& src);
			//copy constructor, retains same cpid and identical characteristics



		ControlParameter& set_callback(int s, ParameterSettingCallback psc_func);
			//Associates a programmer-defined callback function to be invoked when the given setting s is to be applied. 
			//Returns the same ControlParameter object for convenient chaining.

		ControlParameter& set_allcallback(ParameterSettingCallback psc_func);
			//Associates the same callback function with every setting.
			//Setting will be passed as an argument (see ParameterSettingCallback definition).
			//Returns the same ControlParameter object for convenient chaining.

		void force_characteristics(int new_N, const std::string& new_description);
			//Force this ControlParameter to take on the new given characteristics, in violation of regular semantics.
			//The new characteristics will become visible to an associated FeatureController *during* the
			//  next invocation of frame_transition().
			//Retains the callbacks relevant for the intersection of the ranges [-N, N] and [-new_N, new_N].

		void force_characteristics(const ControlParameter& srcCP);
			//Force this ControlParameter to take on the characteristics of srcCP, in violation of regular semantics.
			//The new characteristics will become visible to an associated FeatureController *during* the
			//  next invocation of frame_transition().
			//Acquires the callbacks for of srcCP.
	};




	////////////////////////////////////////////////
	//////// Objective definitions 
	////////////////////////////////////////////////

	typedef long long int OB_ID;
		//Unique ID associated with each created Objective

	class MOTerm {
	public:
		double center;
		double incr_benefit;

		MOTerm(double center, double incr_benefit)
			: center(center), incr_benefit(incr_benefit) { }
	};
	//A term in a Maximization-Objective specification

	typedef double (* MeasurementCallback)(void);
		//Function of this type is invoked by the controller in order to read an Objective's value

	class Objective {
	public:
		//Semantics: An Objective gets created with an immutable set of characteristics

		typedef enum {eRange, eMaximization} Type;

		const OB_ID obid;
			//Unique ID associated with each Objective (unique independent of number of controllers created)
			// obid = -1 for an undefined Objective

		Objective();
			//Undefined objective, will not impact control
			// obid = -1
		
		Objective(
			double center,
			double width,
			int sliding_window_length = 1,
			const std::string& description = "",
			double importance = 1.0
		);
		//Range-objective => get_type() == eRange
		//
		//sliding_window_length > 1 indicates that sliding_window_length number of
		//  past measurements must be averaged for use by a FeatureController
		//  associated with this Objective

		Objective(
			const std::vector<MOTerm>& v_moterms,
			int sliding_window_length = 1,
			const std::string& description = "",
			double importance = 1.0
		);
		//Maximization-objective => get_type() == eMaximization:
		// centers in the v_moterms must be sorted in ascending order.
		// incr_benefit fields will be normalized before use so they sum up to 1.0

		Objective(const Objective& src);
			//copy constructor, retains same obid and identical characteristics
	
		Type get_type() const;

		void set_measurement_callback(MeasurementCallback msc_func);
			//Set the programmer-defined callback function that the controller can invoke at
			//  the end of each frame to read off the Objective's value

		void force_characteristics(const Objective& srcObj);
			//Force this Objective to take on the characteristics of srcObj, in violation of regular semantics.
			//The new characteristics will become visible to an associated FeatureController *during* the
			//  next invocation of frame_transition()
	};




	////////////////////////////////////////////////
	//////// FeatureController definitions 
	////////////////////////////////////////////////

	class EntityStatistics {
	public:
		long long int num_frames;

		double last_frame_mseq;
		bool last_frame_satisfied;

		double cumulative_frame_mseq;
		double cumulative_satisfaction_ratio;

		typedef std::vector<double> Buckets;
		Buckets buckets;
			//A list of MSEQ deviations from mean-objective.
			//Defines the buckets for Histogram below

		typedef std::vector<long long int> Histogram;

		typedef std::vector<int> InputChoice;
		std::map<InputChoice, Histogram> mCP_histograms;
			//Each element is a histogram of MSEQ-deviations corresponding to a unique InputChoice

		EntityStatistics()
			: num_frames(0), last_frame_mseq(0.0), last_frame_satisfied(false),
				cumulative_frame_mseq(0.0), cumulative_satisfaction_ratio(0.0)
		{}

		std::string print_string(bool bPrintHistograms = true) const;
	};

	class AllStatistics {
	public:
		EntityStatistics frame_stats;
			//Statistics for a frame as a whole

		std::map<OB_ID, EntityStatistics> objective_stats;
			//Statistics maintained for individual objectives

		std::string print_all_string(bool bPrintHistograms = true) const;
	};


	typedef long long int FC_ID;

	class FeatureController {
	public:
		const FC_ID fcid;

		FeatureController();

		//////// Parameter Association ////////
		void add_control_parameter(const ControlParameter& cp);
			//ControlParameter cp will be used to achieve feature-control

		void delete_control_parameter(const ControlParameter& cp);
			//ControlParameter cp will no longer be used for feature-control.
			//Error if ControlParameter cp not associated with controller



		//////// Objective Association ////////
		void add_objective(const Objective& ob);
			//Objective ob will be used to define control-objective

		void delete_objective(const Objective& ob);
			//Objective ob will no longer be used to define control-objective.
			//Error if Objective ob not associated with controller



		//////// Perception Window Association ////////
		void set_perception_window_length(int W);
			//The expected size of the user-perception window (default: W=1 assumed)


		//////// Estimation-Model Structure ////////
		void set_model_x_order(int x_order);
		void set_model_y_order(int y_order);
			//Default: x_order = 1, y_order = 1
			//Requirement:
			//   x_order >= 0
			//   y_order >= 1 (to allow for feedback control)

		int get_model_x_order();
		int get_model_y_order();


		//////// Controller Operation ////////
		void frame_transition();
			//Signals a frame transition:
			//   - Controller will read out all Objectives for the finished frame, via measurement-callbacks
			//   - Controller will apply new ControlParameters for the next frame, via setting-callbacks
			//
			//First invocation of frame_transition():
			//   - Prior invocations of set_perception_window_length(), set_model_x_order(), set_model_y_order(),
			//       add/delete_control_parameter() and add/delete_objective() take effect.
			//
			//Subsequent invocations of frame_transition():
			//   - Additional invocations of set_perception_window_length(), set_model_x_order(), set_model_y_order(),
			//       add/delete_control_parameter() and add/delete_objective() DO NOT take effect.
			//   - Invoke reinitialize() subsequently to force these additional calls to take effect.

		void reinitialize();
			//See description under frame_transition()

		long long int get_frame_number() const;
			//The current frame-number, incremented by each call to frame_transition().
			// Starts at -1.

		void set_frame_budget(double time_budget);
			//Budgeted time (in seconds) each call to frame_transition() can spend on determining an
			//  appropriate control-response to observed application behavior
			//A frame-budget must be set via set_frame_budget() before the first call to frame_transition()

		void set_minimum_history_length(int min_length);
			//The minimum number of samples of history that must be retained for model estimation.
			//The featurecontroller will retain at least W samples, even if not specified by user.

		const AllStatistics& get_allstatistics() const;
	
		//////// OVERRIDING Controller Operation (Debugging, Testing) ////////
		void force_fixed_control_parameters(
			const std::vector< std::pair<CP_ID, int> >& vCP_fixedSetting_pairs
		);
			//Forces fixed settings to be applied for each ControlParameter, disables controller:
			//   FeatureController will have extremely low runtime overhead under forced settings
			//
			//vCP_fixedSetting_pairs has elements for each ControlParameter associated with this FeatureController:
			//  each element is a pair of the ControlParameter's CP_ID and the fixed setting that must be used
			//  for that ControlParameter (subject to -N to +N constraint of the ControlParameter)


    void override_control_parameter_in_next_frame(
        CP_ID controlparm_id,
        int   override_value
    );
      //In next frame, ignore the value produced by the controller for ControlParameter controlparm_id.
      //Instead apply override_value. Call has effect only for one frame.
	private:
		FeatureController(const FeatureController& fc);
			//copy-construction disallowed

		FeatureController& operator=(const FeatureController& fc);
			//assignment disallowed
	};


	////////////////////////////////////////////////
	//////// Convenient Timing-Objective definitions 
	////////////////////////////////////////////////

	void associate_timing_callback_with_objective(Objective& timingObj);
		//Internally associates a timing-measurement-callback with the given Objective

	void start_timing_measurement(OB_ID obid);
		//Starts a timing measurement for Objective with given OB_ID.
		//associate_timing_callback_with_objective() must have previously been invoked for this Objective.

	void end_timing_measurement(OB_ID obid);
		//Ends a timing measurement for Objective with given OB_ID.
		//associate_timing_callback_with_objective() must have previously been invoked for this Objective.
		//Time since previous start_timing_measurement() on same obid will get accumulated
		//  to Objective obid's measurement.
		//Objective obid's measurement will be zeroed out after each frame_transition() invocation.
		//start-end can be invoked in pairs to accumulate multiple disjoint time-segments within a frame.
		//Invoking end_timing_measurement() without corresponding start_timing_measurement() will have no effect.


} //namespace fctrl

#endif //__FCTRL_H__
