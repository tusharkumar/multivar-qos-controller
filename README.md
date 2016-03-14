# Multi-Variate QoS Controller

## Purpose
A C++ library and API intended to be applied to "frame-oriented" immersive
applications that satisfy the "immersive assumptions". The controller is
applied to an application via API calls. The controller will adjust multiple
programmer-identified tunable parameters "X" at the beginning of a frame,
in an attempt to keep programmer-identified QoS outcomes "Y" is desired ranges.

## Dependences
Armadillo is used to execute math faster.

http://arma.sourceforge.net/

## Building
1. Modify the Makefile to provide the path to the Armadillo install directory.

2. Build libfctrl.a. This library will need to be linked with a user application.

    make

3. Basic test -- produces test program test.exe

    make test

## How to add controller to your application
1. Include "fctrl.h"

2. Link against libfctrl.a and libarmadillo.a.

3. See "multivar_test.cpp" as a simple example

4. Create controller, register application and execute
```C++ 
    // Create controller object
    fctrl::FeatureController fc;

    // Set frame budget -- e.g., 2% of the frame-time
    fc.set_frame_budget(frame_time_objective * 0.02);

    // Create QoS Objective -- repeat for each objective
    //  - Create objective object and set target: (desired, range) -- e.g., range = 30% of desired
    //  - Assign a callback function (returns double) that controller
    //    can call at the end of a frame to sample a QoS objective value
    //  - Add the objective object to the controller
    fctrl::Objective ob1(frame_time_objective, frame_time_objective*0.30);
    ob1.set_measurement_callback(measure_time);
    fc.add_objective( ob1 );

    // Create a tunable control parameter -- repeat for each control parameter
    //  - Create parameter object, specify "N" for integral range -N to N
    //  - Set callback function -- accepts an int x  (-N <= x <= N) from controller
    //    to apply to the application at the start of a frame
    //  - Add control parameter to the controller
    fctrl::ControlParameter cp1(2); // N = 2
    cp1.set_allcallback(Xset_callback);
    fc1.add_control_parameter(cp1);

    // Execute a frame transition
    // - QoS objectives will be sampled and
    //   Control parameters will be set for the next frame
    fc.frame_transition();

    // Print statistics about the controller's performance at any time
    std::cout << "STATS: frame#" << fc.get_frame_number() << " "
              << fc.get_allstatistics().print_all_string() << std::endl;
```
