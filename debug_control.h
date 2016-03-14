#ifndef __DEBUG_CONTROL_H__
#define __DEBUG_CONTROL_H__

#include "fctrl.h"

namespace fctrl {

	void set_model_order(FeatureController& fc, int x_order, int y_order);
	int get_current_model_x_order(const FeatureController& fc);
	int get_current_model_y_order(const FeatureController& fc);

	//TODO: something to allow/control/disallow model-exploration policy


		//FIXME: following unimplemented
	void set_LLSE_lambda(FeatureController& fc, double lambda);
	double get_LLSE_lambda(const FeatureController& fc);

	//TODO: something to allow/control/disallow lambda variation and selection


} //namespace fctrl

#endif //__DEBUG_CONTROL_H__
