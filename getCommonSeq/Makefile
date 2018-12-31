objects = ComSeq.o  ClusterToFeature.o
CC = g++  -Wformat-extra-args -std=c++11  -g   #-Wall
ARMA_INCLUDE_FLAG = -I ../LinearAlgebra/include  -I ../ 
LIB_FLAGS = -lblas -llapack 
CXXFLAGS = $(ARMA_INCLUDE_FLAG)
install : $(objects)

ComSeq.o  : ComSeq.cpp ComSeq.h
	$(CC)$(CXXFLAGS)  -c  $<   $(LIB_FLAGS) 
ClusterToFeature.o : ClusterToFeature.cpp ComSeq.h 
	$(CC)$(CXXFLAGS)  -c  $<   $(LIB_FLAGS)  

clean :
	rm  $(objects) 

