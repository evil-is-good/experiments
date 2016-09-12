#===============================================================================
#
#  TypeProject: obj
#  Description:
#
#===============================================================================

ProjectName := cgal_grid_generator

#==============================CONFIGURATION: debug/release/miscellaneous...===#

ConfigurationName := release

#================================================================DIRECTORIES===#

ProjectDir := /home/primat/projects/experiments/bloks/src/cgal_grid_generator

#==================================================INCDIRS (Format: IncDir*)===#

IncDirProjects = $(SupDirHome)
IncDirDealII   = $(SupDirDealII)/include
IncDirBoost    = $(SupDirDealII)/contrib/boost-1.49.0/include
IncDirTBB      = $(SupDirDealII)/contrib/tbb/tbb30_104oss/include 

#==================================================LIBDIRS (Format: LibDir*)===#

LibDirDealII = $(SupDirDealII)/lib
#LibDirZLib   = /usr/lib/ruby/1.9.1/x86_64-linux

#===================================================USER DIRS (Format UDir*)===#

Udir = 

#===========================================SUPPORTING DIRS (Format SupDir*)===#

SupDirDealII = /home/primat/deal.II
SupDirHome   = /home/primat

#======================================================================FILES===#

Executable := $(ProjectName).o

Sources    := $(ProjectName).cpp

DebugLibraryes    := deal_II.g tbb_debug

ReleaseLibraryes  := deal_II tbb gomp pthread CGAL_Core CGAL gmpxx mpfr gmp boost_thread-mt z stdc++

#=================================================================BASE FLAGS===#

DebugFlagsC   := -O0 -gstabs+ -pipe -std=c++11

ReleaseFlagsC := -std=c++11 -O2 -pipe -fopenmp -frounding-math

DebugFlagL    := -gstabs+

ReleaseFlagL  := 

AnalysisFlags := #-Weffc++ \
	             -Woverloaded-virtual -Wctor-dtor-privacy \
                 -Wnon-virtual-dtor -Wold-style-cast -Wconversion \
				 -Wsign-conversion -Winit-self -Wunreachable-code 

#=================================================================USER FLAGS===#

UserFlags :=
CGAL_libs :=

#==============================================================================#

ObjDir     := $(ProjectDir)
SourcesDir := $(ProjectDir)

CFlags :=
LFlags :=
ifeq ($(ConfigurationName), debug)
	CFlags := $(DebugFlagsC)
	LFlags := $(DebugFlagsL)
	Libs   := $(addprefix -l, $(DebugLibraryes))
endif
ifeq ($(ConfigurationName), release)
	CFlags := $(ReleaseFlagsC)
	LFlags := $(ReleaseFlagsL)
	Libs   := $(addprefix -l, $(ReleaseLibraryes))
endif

buff1   := $(filter IncDir%, $(.VARIABLES))
buff2   := $(foreach Dir, $(buff1), $($(Dir)))
IncDirs := $(addprefix -I, $(buff2)) -I.

buff1   := $(filter LibDir%, $(.VARIABLES))
buff2   := $(foreach Dir, $(buff1), $($(Dir)))
LibDirs := $(addprefix -L, $(buff2)) -L.


Objects    := $(patsubst %.cpp, %.o, $(Sources))
Objects    := $(addprefix $(ProjectDir), $(Objects))

Sources    := $(addprefix $(ProjectDir), $(Sources))

Executable := $(addprefix $(ProjectDir)/$(ConfigurationName)/, $(Executable)) 

#main:
#	@echo $(ObjDir)/Test1.o

$(ProjectName).o: $(ProjectName).cpp
	g++ -c $(CFlags) $< -o ../../obj/$@ $(IncDirs) $(LibDirs) $(Libs) 2> ~/projects/code.log 

analysis: $(Sources)
	$(CXX) $(AnalysisFlags) $^ $(IncDirs) $(LibDirs) $(Libs)

clean:
	$(RM) $(Objects) $(Executable)
	
# cgal: 
# 	$(CXX) $(CFlags) cgal.cpp -o cgal.exe $(IncDirs) $(LibDirs) $(Libs) 2> ~/projects/code.log
# 
# tria: 
# 	$(CXX) $(CFlags) tria.cpp -o tria.exe $(IncDirs) $(LibDirs) $(Libs) 2> ~/projects/code.log
	
.PHONI: run clean cgal tria
