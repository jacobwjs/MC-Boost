# CFLAGS for running
#CFLAGS = -Wall -v -O3

# CFLAGS for debugging
CFLAGS = -Wall -v -O0 -g

CC = g++
RM =/bin/rm -rf
LIBS =-lboost_thread


OBJS = main.o photon.o layer.o boundary.o medium.o absorber.o sphereAbsorber.o cylinderAbsorber.o vector3D.o logger.o detector.o circularDetector.o pressureMap.o
#OBJS := $(wildcard *.o)

.cpp.o:
	 $(CC) -c $(CFLAGS) $*.cpp


all : mc-boost 


mc-boost: $(OBJS)
	 $(CC) -o  $@ $(OBJS) $(CFLAGS) $(LIBS)


clean::
	 $(RM) mc-boost
	 $(RM) *.o

