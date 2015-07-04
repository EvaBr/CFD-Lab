#### The petsc environment can also be fixed in the makefile
# PETSC_DIR = 
# PETSC_ARCH = 
include ${PETSC_DIR}/conf/variables

# default gnu compiler (currently not used)
# CC = g++
# compiler wrapper for mac-cluster
#CC = mpiCC
# compiler on Ubuntu
#
#
CC = mpic++
CFLAGS = -Wall -g -Werror -O3
SRCDIR = ./
TESTDIR = ./tests/
BUILDDIR = build
INCLUDE = -I. -Istencils -Itests -Iturbulence_algebraic ${PETSC_CC_INCLUDES}


NSMAIN = main.o
TESTMAIN = tests/flow_field_test.o

OBJ = DataStructures.o Configuration.o 3rdparty/tinyxml2/tinyxml2.o SimpleTimer.o

NSOBJ = FlowField.o LinearSolver.o stencils/MaxUStencil.o stencils/MovingWallStencils.o stencils/PeriodicBoundaryStencils.o\
stencils/FGHStencil.o solvers/SORSolver.o solvers/PetscSolver.o \
stencils/RHSStencil.o stencils/VTKStencil.o stencils/VelocityStencil.o \
parallelManagers/PetscParallelManager.o parallelManagers/BufferStencils.o \
parallelManagers/PetscParallelConfiguration.o\
GlobalBoundaryFactory.o\
stencils/BFStepInitStencil.o stencils/NeumannBoundaryStencils.o stencils/BFInputStencils.o\


all: ns

ns: $(OBJ) $(NSOBJ) $(NSMAIN)
	$(CC) -o ns $(OBJ) $(NSOBJ) $(NSMAIN) $(PETSC_KSP_LIB) -lstdc++ $(CFLAGS)

test: $(OBJ) $(NSOBJ) $(TESTMAIN)
	$(CC) -o testbin $(OBJ) $(NSOBJ) $(TESTMAIN) $(PETSC_KSP_LIB) -stdc++ $(CFLAGS)

%.o: %.cpp
	$(CC) -c $(CFLAGS) $(INCLUDE) -o $*.o $*.cpp $(PETSC_KSP_LIB) -lstdc++

cleanall clean:
	for name in  ns main.o $(NSOBJ) testbin $(OBJ) tests/*.o; do \
	if [ -f $$name ]; then rm $$name; fi; \
	done;
