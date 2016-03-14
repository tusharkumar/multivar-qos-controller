#include <iostream>
#include <cassert>
#include <cstdlib>
#include <time.h>

#include "fctrl.h"
#include "timing.h"

void f1(int x) {
	std::cout << "Inside f1" << std::endl;
	for(int i=0; i<1000000; i++)
		x = (x * 10 + 5) / 9;
	std::cout << "x = " << x << std::endl;
}

void f2(int x) {
	std::cout << "Inside f2" << std::endl;
	for(int i=0; i<500000; i++)
		x = (x * 10 + 5) / 9;
	std::cout << "x = " << x << std::endl;
}

void f3(int x) {
	std::cout << "Inside f3" << std::endl;
	for(int i=0; i<100000; i++)
		x = (x * 10 + 5) / 9;
	std::cout << "x = " << x << std::endl;
}

void f4(int x) {
	std::cout << "Inside f4" << std::endl;
	for(int i=0; i<50000; i++)
		x = (x * 10 + 5) / 9;
	std::cout << "x = " << x << std::endl;
}


//Objective#0
timeval last_starttime;
double measure_time(void) {
	timeval endtime = fctrl::get_curr_timeval();
	double diff = fctrl::diff_time(last_starttime, endtime);
	last_starttime = endtime;

	return diff;
}

//ControlParameter#0
int f_choice  = 100;
void Xset_callback(int s) {
	f_choice = s;
}

//Objective#1
double goalval = -100.0;
double get_goalval(void)
{ return goalval; }

//ControlParameter#1
void set_goalval(int s) {
	goalval = 5.0 * s;
#if 1
	f_choice -= s;
	if(f_choice > 2)
		f_choice = 2;
	if(f_choice < -2)
		f_choice = -2;
#endif
}


void run_chosen_f(int x)
{
	if (f_choice == -2)
		f1(x);
	else if(f_choice == -1)
		f2(x);
	else if(f_choice == 0)
		f3(x);
	else if(f_choice == 1 or f_choice == 2)
		f4(x);
	else
	{ assert(0); }
}


//fctrl::FeatureController fc;

void ww() {
	std::cout << "Inside f_ww" << std::endl;

	int x = 5;
	run_chosen_f(x);
}

int main() {
	std::srand(time(NULL));

	fctrl::FeatureController fc;

	double frame_time_objective = 0.005;
	fc.set_frame_budget(frame_time_objective * 0.16);
	fc.set_model_x_order(0);
	fc.set_model_y_order(1);

	//Objective#0
	fctrl::Objective ob(frame_time_objective, frame_time_objective*0.3);
	ob.set_measurement_callback(measure_time);
	fc.add_objective( ob );

	//ControlParameter#0
	fctrl::ControlParameter cp(2);
	cp.set_allcallback(Xset_callback);
	fc.add_control_parameter(cp);

	
	/*
	//Objective#1
	fctrl::Objective ob_goal(5, 5*0.2);
	ob_goal.set_measurement_callback(get_goalval);
	fc.add_objective( ob_goal );
	

	
	//ControlParameter#1
	fctrl::ControlParameter cp_goal(3);
	cp_goal.set_allcallback(set_goalval);
	fc.add_control_parameter(cp_goal);
	*/
	
	std::vector< std::pair<fctrl::CP_ID, int> > v;
	v.push_back( std::make_pair(cp.cpid, -2) );
	//v.push_back( std::make_pair(cp_goal.cpid, 0) );
	//fc.force_fixed_control_parameters(v);

	last_starttime = fctrl::get_curr_timeval();
	for(int i=0; i<100; i++) {
		std::cout << "---------- ENDING frame#" << fc.get_frame_number() << " ------------" << std::endl;
		//if(i%5 == 0)
		//	fc.override_control_parameter_in_next_frame(cp.cpid, 2);
		fc.frame_transition();
		std::cout << "STATS: frame#" << fc.get_frame_number() << " " << fc.get_allstatistics().print_all_string() << std::endl;
		ww();
	}

	return 0;
}
