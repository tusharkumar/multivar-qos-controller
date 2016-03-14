#ifndef __FRAME_OBSERVATION_H__
#define __FRAME_OBSERVATION_H__

#include "fctrl.h"

namespace fctrl {

	class FrameObservation {
	public:
		long long int frame_number;
			//observations for which frame-number

		std::vector<double> vCP_observations;
			//ControlParameter observations for frame
		std::vector<double> vOB_observations;
			//Objective observations for frame

		FrameObservation()
			: frame_number(-1)
		{}
	};

} //namespace fctrl

std::ostream& operator<<(std::ostream& os, const fctrl::FrameObservation& fo);


#endif //__FRAME_OBSERVATION_H__
