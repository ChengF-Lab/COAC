objects =getCommonSeq/ComSeq.o getCommonSeq/ClusterToFeature.o
ARMA_INCLUDE_FLAG = -I  LinearAlgebra/include
LIB_FLAGS = -lblas -llapack 
CXXFLAGS = $(ARMA_INCLUDE_FLAG)
CC = g++  -std=c++11  -g  -I getCommonSeq/    #-Wall
all:    prepareall
	$(CC)$(CXXFLAGS)   -o   Cluster    main.cpp  $(objects) $(LIB_FLAGS)
	$(CC)$(CXXFLAGS)   -o   FeatureGen main_FeatureGen.cpp  $(objects) $(LIB_FLAGS) 
prepareall:    subsystem
subsystem:
	$(MAKE) -C getCommonSeq/
clean :  cleansub

	rm     Cluster  FeatureGen
cleansub :
	rm  $(objects)

