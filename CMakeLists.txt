cmake_minimum_required(VERSION 3.30)
project(SIDM_order_nodrift_demao)

set(CMAKE_CXX_STANDARD 14)

# 1) Find MPI for C
find_package(MPI REQUIRED C)

# 2) Find GSL
find_package(GSL REQUIRED)

include_directories(simtemplate/codeSIDM1)
include_directories(simtemplate/input)

add_executable(SIDM_order_nodrift_demao
        simtemplate/codeSIDM1/accel.c
        simtemplate/codeSIDM1/allocate.c
        simtemplate/codeSIDM1/allvars.c
        simtemplate/codeSIDM1/allvars.h
        simtemplate/codeSIDM1/begrun.c
        simtemplate/codeSIDM1/check.C
        simtemplate/codeSIDM1/density.c
        simtemplate/codeSIDM1/domain.c
        simtemplate/codeSIDM1/doscatt.c
        simtemplate/codeSIDM1/driftfac.c
        simtemplate/codeSIDM1/endrun.c
        simtemplate/codeSIDM1/forcetree.c
        simtemplate/codeSIDM1/global.c
        simtemplate/codeSIDM1/gravtree.c
        simtemplate/codeSIDM1/gravtree_forcetest.c
        simtemplate/codeSIDM1/hydra.c
        simtemplate/codeSIDM1/init.c
        simtemplate/codeSIDM1/io.c
        simtemplate/codeSIDM1/longrange.c
        simtemplate/codeSIDM1/main.c
        simtemplate/codeSIDM1/ngb.c
        simtemplate/codeSIDM1/ngb_sidm.c
        simtemplate/codeSIDM1/peano.c
        simtemplate/codeSIDM1/pm_nonperiodic.c
        simtemplate/codeSIDM1/pm_periodic.c
        simtemplate/codeSIDM1/potential.c
        simtemplate/codeSIDM1/predict.c
        simtemplate/codeSIDM1/proto.h
        simtemplate/codeSIDM1/read_ic.c
        simtemplate/codeSIDM1/restart.c
        simtemplate/codeSIDM1/run.c
        simtemplate/codeSIDM1/sidm.c
        simtemplate/codeSIDM1/system.c
        simtemplate/codeSIDM1/tags.h
        simtemplate/codeSIDM1/timestep.c
        simtemplate/input/definitions.h
        simtemplate/input/functions.c
        simtemplate/input/functions.h
        simtemplate/input/io.c
        simtemplate/input/io.h
        simtemplate/input/routines.c
        simtemplate/input/routines.h
        simtemplate/input/spheric.c)

# Link your target with MPI (C interface) and GSL libraries
target_link_libraries(
        SIDM_order_nodrift_demao
        PRIVATE
        MPI::MPI_C         # From find_package(MPI ...)
        GSL::gsl           # From find_package(GSL ...)
        GSL::gslcblas      # GSL's CBLAS component
)