#===============================================================================
#
#  TypeProject: exe
#  Description:
#
#===============================================================================

ProjectName := extension_speciment_test

#==============================CONFIGURATION: debug/release/miscellaneous...===#

ConfigurationName := release

#================================================================DIRECTORIES===#

ProjectDir := /home/primat/projects/tests/extension_speciment_test

#==================================================INCDIRS (Format: IncDir*)===#

IncDirProjects = $(SupDirHome)
IncDirDealII   = $(SupDirDealII)/include
IncDirBoost    = $(SupDirDealII)/contrib/boost-1.49.0/include
IncDirTBB      = $(SupDirDealII)/contrib/tbb/tbb30_104oss/include 

#==================================================LIBDIRS (Format: LibDir*)===#

LibDirDealII = $(SupDirDealII)/lib

#===================================================USER DIRS (Format UDir*)===#

Udir = 

#===========================================SUPPORTING DIRS (Format SupDir*)===#

SupDirHome   = /home/primat
SupDirDealII = /home/primat/deal.II

#======================================================================FILES===#

Executable := $(ProjectName).exe

Sources    := $(ProjectName).cpp

DebugLibraryes    := deal_II.g tbb_debug

# ReleaseLibraryes  := tbb gomp pthread gmpxx mpfr boost_thread-mt z

ReleaseLibraryes  := deal_II tbb gomp pthread CGAL_Core CGAL gmpxx mpfr gmp boost_thread z

# StaticLibs := /usr/lib/libCGAL_Core.a /usr/lib/libCGAL.a /usr/lib/x86_64-linux-gnu/libgmp.a /home/primat/deal.II/deal.II_bo_share/lib/libdeal_II.a 
StaticLibs := 

#=================================================================BASE FLAGS===#

DebugFlagsC   := -O2 -gstabs+ -pipe -std=c++11

ReleaseFlagsC := -std=c++11 -O2 -pipe -fopenmp -frounding-math

DebugFlagL    := -gstabs+

ReleaseFlagL  := 

AnalysisFlags := #-Weffc++ \
	             -Woverloaded-virtual -Wctor-dtor-privacy \
                 -Wnon-virtual-dtor -Wold-style-cast -Wconversion \
				 -Wsign-conversion -Winit-self -Wunreachable-code 

#=================================================================USER FLAGS===#

UserFlags :=

#==============================================================================#

ObjDir     := $(ProjectDir)/obj
SourcesDir := $(ProjectDir)/sources

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
Objects    := $(addprefix $(ProjectDir)/obj/, $(Objects))

Sources    := $(addprefix $(ProjectDir)/sources/, $(Sources))

Executable := $(addprefix $(ProjectDir)/$(ConfigurationName)/, $(Executable)) 

#main:
#	@echo $(ObjDir)/Test1.o
	
$(Executable): $(Objects)
	$(CXX) $(LFlags) $^ /home/primat/projects/cae/test/grid/grid.o $(StaticLibs) -o $@ $(LibDirs) $(Libs) 2> /home/primat/projects/code.log 

$(ObjDir)/%.o: $(SourcesDir)/%.cpp
	mkdir -p $(dir $@)
	$(CXX) -c $(CFlags) $< -o $@ $(IncDirs) $(LibDirs) $(Libs) 2> /home/primat/projects/code.log 

analysis: $(Sources)
	$(CXX) $(AnalysisFlags) $^ $(IncDirs) $(LibDirs) $(Libs)  

run:
	$(Executable)

clean:
	$(RM)  $(Objects) $(Executable)

chead:
	$(CXX) $(CFlags) -c pch.h $(IncDirs) $(LibDirs) $(Libs)

.PHONI: run clean
