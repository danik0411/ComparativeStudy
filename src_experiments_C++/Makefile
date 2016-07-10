# MAKEFILE for mst
# 
# Usage:
#        make     - For standalone executable
#        make mpi - For mpi compatible executable
#

###a := $(findstring $MAKECMDGOALS, roc)

###ifneq (,$(findstring $(MAKECMDGOALS), roc))
    ###$(info I am in found)
###else
    ###$(info I am not in found)
###endif

$(info $$a is [${a}])

$(info $$MAKECMDGOALS is [${MAKECMDGOALS}])


ifeq ($(MAKECMDGOALS),mpi)
   OPT = $(INCL)  -O2 -g3 -DUSE_MPI
   CC  = mpiCC
else ifeq ($(MAKECMDGOALS),mpiroc)
   $(info i am in time)
   OPT = $(INCL)  -O2 -g3 -DUSE_MPI -DCALC_ROC
   CC  = mpiCC
else ifeq ($(MAKECMDGOALS),time)
   $(info i am in time)
   OPT = $(INCL)  -O2 -g3 -DREC_TIME
   CC  = g++
else ifeq ($(MAKECMDGOALS),timeroc)
$(info i am in time and roc)
   OPT = $(INCL)  -O2 -g3 -DREC_TIME -DCALC_ROC
   CC  = g++
else ifeq ($(MAKECMDGOALS),roctime)
$(info i am in roc and time)
   OPT = $(INCL)  -O2 -g3 -DREC_TIME -DCALC_ROC
   CC  = g++
else ifeq ($(MAKECMDGOALS),roc)
$(info i am in roc)
   OPT = $(INCL)  -O2 -g3 -DCALC_ROC
   CC  = g++
else
   OPT = $(INCL)  -O2 -g3
   CC  = g++
endif


EXECUTABLE=../optgamma

mpi: $(EXECUTABLE)
time: $(EXECUTABLE)
all: $(EXECUTABLE)
timeroc: $(EXECUTABLE)
roctime: $(EXECUTABLE)
mpiroc: $(EXECUTABLE)
roc: $(EXECUTABLE)

$(EXECUTABLE) : \
main.o \
parseArguments.o \
print.o \
utilities.o \
matrixOperations.o \
random.o \
dataGeneration.o \
calculateErrors.o \
featureSelection.o \
classifiers.o \
errorEstimation.o \
ROC.o
	$(CC) -o $(EXECUTABLE) *.o -lm $(OPT)

main.o : main.cpp
	$(CC) -c main.cpp -o main.o $(OPT)

parseArguments.o : parseArguments.cpp
	$(CC) -c parseArguments.cpp -o parseArguments.o $(OPT)

print.o : print.cpp
	$(CC) -c print.cpp -o print.o $(OPT)

utilities.o : utilities.cpp
	$(CC) -c utilities.cpp -o utilities.o $(OPT)

matrixOperations.o : matrixOperations.cpp
	$(CC) -c matrixOperations.cpp -o matrixOperations.o $(OPT)

random.o : random.cpp
	$(CC) -c random.cpp -o random.o $(OPT)

dataGeneration.o : dataGeneration.cpp
	$(CC) -c dataGeneration.cpp -o dataGeneration.o $(OPT)

calculateErrors.o : calculateErrors.cpp
	$(CC) -c calculateErrors.cpp -o calculateErrors.o $(OPT)

featureSelection.o : featureSelection.cpp
	$(CC) -c featureSelection.cpp -o featureSelection.o $(OPT)

classifiers.o : classifiers.cpp
	$(CC) -c classifiers.cpp -o classifiers.o $(OPT)

errorEstimation.o : errorEstimation.cpp
	$(CC) -c errorEstimation.cpp -o errorEstimation.o $(OPT)

ROC.o : ROC.cpp
	$(CC) -c ROC.cpp -o ROC.o $(OPT)

clean: 
	rm -f *.o
