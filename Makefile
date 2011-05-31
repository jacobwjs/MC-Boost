CFLAGS = -Wall
CC = g++
DEBUG = -g
RM =/bin/rm -rf
LIBS =-lboost_thread


OBJS = main.o photon.o layer.o boundary.o medium.o
#OBJS := $(wildcard *.o)

.cpp.o:
	 $(CC) -c $(CFLAGS) $*.cpp
#####


all : mc-boost 


mc-boost: $(OBJS)
	 $(CC) -o  $@ $(OBJS) $(CFLAGS) $(LIBS)


clean::
	 $(RM) mc-boost
	 $(RM) *.o

