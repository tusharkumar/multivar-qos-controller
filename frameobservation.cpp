#include "utils.h"
#include "frameobservation.h"

namespace fctrl {
} //namespace fctrl


std::ostream& operator<<(std::ostream& os, const fctrl::FrameObservation& fo)
{
	os << fctrl::align() << "frame_number=" << fo.frame_number << std::endl;
	os << fctrl::align() << "vCP_observations=" << fo.vCP_observations << std::endl;
	os << fctrl::align() << "vOB_observations=" << fo.vOB_observations << std::endl;

	return os;
}

