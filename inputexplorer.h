#ifndef __INPUT_EXPLORER_H__
#define __INPUT_EXPLORER_H__

#include "modelstructure.h"

namespace fctrl {
	class InputExplorer {
		//   + Deterministic exploration over input space.
		//   + initial step-size = +0.5 (in normalized units)
		//   + Each time-step, add step-size to the next input parameter.
		//   + Once all inputs exhausted, set step-step = -step-size/2. Repeat from first input.
		//   + Terminate process when at least one input doesn't change at all
		//       (due to too fine a step-size compared to its discreteness, or clipping out-of-bounds)
		//   + Then, reset step_size = -0.5. Repeat entire process (then reset step_size = +0.5, etc).
	public:
		ModelStructure * ms;
			//The ModelStructure for which input exploration is to be performed.

		int next_input_index;
		double step_size;

		bool bTerminate;
		int terminate_restart_sign;

		InputExplorer(ModelStructure * ms = 0)
			: ms(ms), next_input_index(0), step_size(0.5),
				bTerminate(false), terminate_restart_sign(-1)
		{}

		arma::colvec gen_next_exploration_input(const arma::colvec& curr_control_input)
		{
			arma::colvec next_exploration_input = curr_control_input;

			assert(ms != 0);
			assert(ms->vActive_CP_IDs.size() > 0);
			if(next_input_index >= (int)ms->vActive_CP_IDs.size()) { //finished inputs
				next_input_index = 0;
				step_size = -step_size/2;

				if(bTerminate) {
					step_size = terminate_restart_sign * 0.5;
					terminate_restart_sign = -terminate_restart_sign;
					bTerminate = false;
				}
			}

			assert(next_input_index >= 0 and next_input_index < (int)next_exploration_input.n_rows);
			next_exploration_input(next_input_index) += step_size;
			if(next_exploration_input(next_input_index) > 1.0)
				next_exploration_input(next_input_index) = 1.0;
			else if(next_exploration_input(next_input_index) < -1.0)
				next_exploration_input(next_input_index) = -1.0;

			std::vector<int> CP_choices_curr = ms->convert_normalized_to_clipped_CP_choices(curr_control_input);
			std::vector<int> CP_choices_next = ms->convert_normalized_to_clipped_CP_choices(next_exploration_input);

			if(CP_choices_curr.at(next_input_index) == CP_choices_next.at(next_input_index))
				bTerminate = true;

			next_input_index++;
			return next_exploration_input;
		}
	};
} //namespace fctrl

#endif //__INPUT_EXPLORER_H__
