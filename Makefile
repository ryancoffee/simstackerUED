CXX 		= g++
CC 		= $(CXX)
EXEC      	= non_obj_propagator
SRC_FILES 	= data2d.cpp non_obj_propagator.cpp 
DEBUG_LEVEL     = -O3 -g
#EXTRA_CCFLAGS   = -Wall -c -DHAVE_INLINE -D_USE_MATH_DEFINES -std=gnu++11
EXTRA_CCFLAGS   = -W -c -DHAVE_INLINE -D_USE_MATH_DEFINES -std=gnu++14
INCLUDE_PATHS   = -I $(HOME)/computing/boost # $(HOME)/computing/wavelib/linuxshared 
CXXFLAGS        = $(DEBUG_LEVEL) $(EXTRA_CCFLAGS) $(INCLUDE_PATHS) $(LDLIBS) $(LDFLAGS)
CCFLAGS         = $(CXXFLAGS)
LDLIBS		= -lgsl -lgslcblas -lfftw3 -lm #-lopencv_core -lm
LDFLAGS         = -L. -L/usr/local/lib 
O_FILES         = $(SRC_FILES:%.cpp=%.o)
H_FILES		= $(SRC_FILES:%.cpp=%.hpp)

all: $(EXEC)
$(EXEC): $(O_FILES)

depend:
	makedepend -- $(CXXFLAGS) -- -Y $(SRC_FILES) $(H_FILES)

clean:
	$(RM) $(O_FILES) core *.rpo

very-clean: clean
	$(RM) $(EXEC)

# DO NOT DELETE

data2d.o: data2d.hpp
non_obj_propagator.o: ./non_obj_propagator.hpp Constants.hpp ./data2d.hpp
non_obj_propagator.o: Constants.hpp
