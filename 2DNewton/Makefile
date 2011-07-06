# File: Makefile for Linux
# Author: J.K. Johnstone, Patrick Butler
# Last Modified: 02/09/10

ARCH	   = LINUX
SHELL      = /bin/bash
CC         = g++
HOME	   = /media/Data/Code/Grail/2DNewton/2DNewton/
LIBRARIES  = -lglut -lGLU -lGL -lm
# LIBRARIES = -lm
LDFLAGS    = -funroll-all-loops -O3

all: 2ddriver
.cpp: 
	${CC} -Wno-deprecated -g -o $@ $< ${LIBRARIES} ${LDFLAGS}
clean :
	rm -f *.o *~ a.out core
	ls -l

touch :
	touch 2ddriver.cpp

