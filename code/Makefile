##################################################################################################
# makefile for a simple project
# available targets:
#      make
#      make clean
##################################################################################################
#  Edit the following lines according to your needs
##################################################################################################
PROGRAM	= main_prc
SOURCES = main_prc.cpp io_data.cpp initialize.cpp observables.cpp parameters.cpp
OBJECTS = $(SOURCES:.cpp=.o)

PROGRAM2 = SU2_start_config_gen_prc
SOURCES2 = SU2_start_config_gen_prc.cpp io_data.cpp initialize.cpp observables.cpp parameters.cpp
OBJECTS2 = $(SOURCES2:.cpp=.o)

PROGRAM3 = main_sim
SOURCES3 = main_sim.cpp io_data.cpp initialize.cpp observables.cpp parameters.cpp
OBJECTS3 = $(SOURCES3:.cpp=.o)

PROGRAM4 = SU2_start_config_gen_sim
SOURCES4 = SU2_start_config_gen_sim.cpp io_data.cpp initialize.cpp observables.cpp parameters.cpp
OBJECTS4 = $(SOURCES4:.cpp=.o)

PROGRAM5 = stats_dos_jn
SOURCES5 = stats_dos_jn.cpp parameters.cpp io_data.cpp
OBJECTS5 = $(SOURCES5:.cpp=.o)

CC	= gcc
CXX = g++
F77	= gfortran

WARNINGS = -Wall -Weffc++ -Woverloaded-virtual -W -Wfloat-equal -Wshadow \
           -Wredundant-decls -Winline
#  -Wunreachable-code
# CXXFLAGS = -ffast-math -O3 -funroll-all-loops -DNDEBUG -msse3 ${WARNING} -ftree-vectorizer-verbose=2
#CXXFLAGS = -ffast-math -O3 -funroll-all-loops -DNDEBUG ${WARNING} -ftree-vectorizer-verbose=2 \
           # -fopenmp
CXXFLAGS = -std=c++11 -march=native -pedantic -Wall -Wextra -Ofast -fopenmp
LDFLAGS = -fopenmp
#CXXFLAGS = -std=c++11 `gsl-config --cflags`
# CXXFLAGS = -g
#LDFLAGS  = -lcblas -fopenmp
#LDFLAGS = `gsl-config --libs`
#LDFLAGS =  -lgslcblas -lm -lgsl

# CFLAGS	= -ffast-math -O3 -DNDEBUG -msse3 -fopenmp -fdump-tree-vect-details
# CFLAGS	= -ffast-math -O3 -funroll-loops -DNDEBUG -msse3 -fopenmp -ftree-vectorizer-verbose=2
# #CFLAGS	= -ffast-math -O3 -DNDEBUG -msse3 -fopenmp
# FFLAGS	= -ffast-math -O3 -DNDEBUG -msse3 -fopenmp
# LFLAGS  = -ffast-math -O3 -DNDEBUG -msse3 -fopenmp

##################################################################################################
#    Don't change anything below this line
##################################################################################################
all_prc:	${PROGRAM}

${PROGRAM}:	${OBJECTS}
	$(CXX) $^ -o $@ ${LDFLAGS}

clean:
	rm -f ${PROGRAM} ${OBJECTS}

.cpp.o:
	$(CXX) -c $(CXXFLAGS) $<

.c.o:
	$(CC) -c $(CFLAGS) $<

.f.o:
	$(F77) -c $(FFLAGS) $<

all_cgprc:	${PROGRAM2}

${PROGRAM2}:	${OBJECTS2}
	$(CXX) $^ -o $@ ${LDFLAGS}

clean_cgprc:
	rm -f ${PROGRAM2} ${OBJECTS2}

all_sim:	${PROGRAM3}

${PROGRAM3}:	${OBJECTS3}
	$(CXX) $^ -o $@ ${LDFLAGS}

clean_sim:
	rm -f ${PROGRAM3} ${OBJECTS3}

all_cgsim:	${PROGRAM4}

${PROGRAM4}:	${OBJECTS4}
	$(CXX) $^ -o $@ ${LDFLAGS}

clean_cgsim:
	rm -f ${PROGRAM4} ${OBJECTS4}

all_stats:	${PROGRAM5}

${PROGRAM5}:	${OBJECTS5}
	$(CXX) $^ -o $@ ${LDFLAGS}

clean_stats:
	rm -f ${PROGRAM5} ${OBJECTS5}



###################################################################################################
#    some tools
# Cache behaviour (CXXFLAGS += -g  tracks down to source lines)
cache: ${PROGRAM}
	valgrind --tool=callgrind --simulate-cache=yes ./$^
#	kcachegrind callgrind.out.<pid> &

# Check for wrong memory accesses, memory leaks, ...
# use smaller data sets
mem: ${PROGRAM}
	valgrind -v --leak-check=yes --tool=memcheck --undef-value-errors=yes --track-origins=yes --log-file=$^.addr.out --show-reachable=yes ./$^

#  Simple run time profiling of your code
#  CXXFLAGS += -g -pg
#  LDFLAGS += -pg
prof: ${PROGRAM}
	./$^
	gprof -b ./$^ > gp.out
#	kprof -f gp.out -p gprof &
