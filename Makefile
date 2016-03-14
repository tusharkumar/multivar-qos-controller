#ARMADILLO_INCLUDE=/home/tushark/builds/armadillo/armadillo_install/usr/include
#ARMADILLO_LIB=/home/tushark/builds/armadillo/armadillo_install/usr/lib

ARMADILLO_INCLUDE=/c/Tushar/experiments/build/armadillo_install/usr/include
ARMADILLO_LIB=/c/Tushar/experiments/build/armadillo_install/usr/lib

#ARMADILLO_INCLUDE=/experiments/build/armadillo/armadillo_install/usr/include
#ARMADILLO_LIB=/experiments/build/armadillo/armadillo_install/usr/lib
##

TARGET_LIB=libfctrl.a
TEST_PROGRAM=test.exe

##

all: $(TARGET_LIB)

test: $(TEST_PROGRAM)

clean:
	rm -f *.o $(TEST_PROGRAM) $(TEST_PROGRAM).stackdump $(TARGET_LIB)

##

ACTUAL_CPP_FILES= \
	utils.cpp \
	controlparameter.cpp \
	objective.cpp \
	featurecontroller.cpp \
	frameobservation.cpp \
	modelstructure.cpp \
	llse.cpp \
	model.cpp \
	linearpredictionmodel.cpp \
	lqr.cpp \
	failuremodes.cpp \
	timing.cpp \
	samplestats.cpp

H_FILES= \
	fctrl.h \
	fctrl_version.h \
	debug_control.h \
	utils.h \
	controlparameter.h \
	objective.h \
	featurecontroller.h \
	frameobservation.h \
	modelstructure.h \
	llse.h \
	model.h \
	linearpredictionmodel.h \
	lqr.h \
	failuremodes.h \
	inputexplorer.h \
	timing.h \
	samplestats.h \
	samplestats_internal.h \
	stabilitycoveragemetrics.h

CFLAGS=-g -Wall -I $(ARMADILLO_INCLUDE)
LFLAGS=-L $(ARMADILLO_LIB) -larmadillo -llapack -lblas

##

CPP_FILES=fctrl_version.cpp $(ACTUAL_CPP_FILES)

O_FILES=$(CPP_FILES:.cpp=.o)

#fctrl_version.cpp: $(H_FILES) $(ACTUAL_CPP_FILES) Makefile
#	echo namespace fctrl { const char* fctrl_version = \"$$(basename $(CURDIR))\"\; } > fctrl_version.cpp

$(TARGET_LIB): $(O_FILES)
	ar -r $(TARGET_LIB) $(O_FILES)

%.o: %.cpp $(H_FILES)
	g++ -c $(CFLAGS) $<

$(TEST_PROGRAM): $(TARGET_LIB) multivar_test.cpp
	g++ $(CFLAGS) multivar_test.cpp -L. -lfctrl $(LFLAGS) -o $(TEST_PROGRAM)

