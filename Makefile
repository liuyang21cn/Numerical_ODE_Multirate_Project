# Makefile for Multirate project
# Yang Liu

# compiler & flags
CXX      = clang++ -std=c++11
CXXFLAGS = -O

# executable targets
all : ProtheroRobinson.exe AdvectionDiffusion.exe ReactionDiffusion.exe ODEsystem.exe Brusselator.exe

ProtheroRobinson.exe : ProtheroRobinson.cpp MBM.cpp MPRK2.cpp MPRK3.cpp mat.cpp 
	$(CXX) $(CXXFLAGS) -o $@ $^

AdvectionDiffusion.exe : AdvectionDiffusion.cpp MBM.cpp MPRK2.cpp MPRK3.cpp mat.cpp 
	$(CXX) $(CXXFLAGS) -o $@ $^

ReactionDiffusion.exe : ReactionDiffusion.cpp MBM.cpp MPRK2.cpp MPRK3.cpp mat.cpp 
	$(CXX) $(CXXFLAGS) -o $@ $^

ODEsystem.exe : ODEsystem.cpp MBM.cpp MPRK2.cpp MPRK3.cpp mat.cpp 
	$(CXX) $(CXXFLAGS) -o $@ $^

Brusselator.exe : Brusselator.cpp MBM.cpp MPRK2.cpp MPRK3.cpp mat.cpp erk4.cpp	
	$(CXX) $(CXXFLAGS) -o $@ $^

# utilities
clean :
	\rm -rf *.txt *.exe *~
