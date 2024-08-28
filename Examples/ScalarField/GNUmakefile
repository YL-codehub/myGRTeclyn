GRTECLYN_HOME = $(realpath ../..)

include $(GRTECLYN_HOME)/Tools/GNUMake/Make.defaults

AMREX_HOME ?= $(realpath ../../../amrex)

include $(AMREX_HOME)/Tools/GNUMake/Make.defs

include ./Make.package
include $(AMREX_HOME)/Src/Base/Make.package
include $(AMREX_HOME)/Src/Boundary/Make.package
include $(AMREX_HOME)/Src/AmrCore/Make.package
include $(AMREX_HOME)/Src/Amr/Make.package
include $(AMREX_HOME)/Src/Extern/SWFFT/Make.package

GRTECLYN_HOME = $(realpath ../..)
GRTECLYN_SOURCE = $(GRTECLYN_HOME)/Source

src_dirs := $(GRTECLYN_SOURCE)/utils \
            $(GRTECLYN_SOURCE)/simd  \
            $(GRTECLYN_SOURCE)/CCZ4  \
            $(GRTECLYN_SOURCE)/Matter \
            $(GRTECLYN_SOURCE)/BoxUtils  \
            $(GRTECLYN_SOURCE)/GRTeclynCore  \
            $(GRTECLYN_SOURCE)/AMRInterpolator  \
            $(GRTECLYN_SOURCE)/TaggingCriteria  \
            $(GRTECLYN_SOURCE)/InitialConditions/ScalarFields  \

src_pack := $(foreach dir, $(src_dirs), $(dir)/Make.package)
include $(src_pack)
CEXE_sources += $(GRTECLYN_CEXE_sources)
CEXE_headers += $(GRTECLYN_CEXE_headers)
VPATH_LOCATIONS += $(src_dirs) $(AMREX_HOME)/Src/Extern/SWFFT
INCLUDE_LOCATIONS += $(realpath .) $(src_dirs)

include $(GRTECLYN_HOME)/Tools/GNUMake/Make.rules


INCLUDE_LOCATIONS += $(AMREX_HOME)/Src/Extern/SWFFT
LIBRARIES += -L$(FFTW_DIR) -lfftw3_mpi -lfftw3_omp -lfftw3
