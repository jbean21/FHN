########################################
# Makefile for PX425 2018 Assignment 5 #
########################################

# Detect whether on my macbook or not
# and change compiler accordingly
OS := $(shell uname)
ifeq ($(OS), Darwin)
	CC=g++-8 	# MacOS (on my macbook, am using g++-8, gcc-8, cc-8 etc)
else
	CC=g++		# Windows/Ubuntu
endif

CFLAGS=-O2 -fopenmp  

# Command to use for linking and executable
LD=$(CC)
EXE=a.out
LDFLAGS=-fopenmp

OBJECTS=SOURCEi.cpp ode.cpp
#OBJECTS=SourceMartin.cpp ode.cpp

# Default build target
perc : $(OBJECTS)
	#-----------------------------------------------------
	# NO!!!!
	# 1. CHANGE OUTPUT FILE FIRST IN 5 PLACES!!!!!
	# 2. CHANGE COUPLING TYPE AT BOTTOM OF SOURCEi.cpp
	# 3. CHANGE NUMBER OF DIMENSIONS
	$(LD) $(OBJECTS) $(LDFLAGS)
	# 5. RE-ARM THE MAKEFILE
	#-----------------------------------------------------

# Purge executable
clean :
	rm $(EXE)
	rm -f *.csv

