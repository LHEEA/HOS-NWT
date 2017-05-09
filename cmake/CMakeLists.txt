CMAKE_MINIMUM_REQUIRED(VERSION 2.8.8)
PROJECT(HOS-NWT C Fortran)

SET(CMAKE_VERBOSE_MAKEFILE OFF)
SET(BUILD_SHARED_LIBS ON)

# Define manually these variables to access fftw and lapack libraries
# FFTW_ROOT
# LAPACK_ROOT

################################################################################
SET(CMAKE_MODULE_PATH   ${CMAKE_CURRENT_SOURCE_DIR}/cmake)
SET(BENCHMARK_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/../Benchmark)
SET(SRC_DIRECTORY       ${CMAKE_CURRENT_SOURCE_DIR}/../sources)
STRING(REPLACE "cmake/../" "" SRC_DIRECTORY ${SRC_DIRECTORY})
MESSAGE(STATUS "Sources to be used are located in ${SRC_DIRECTORY}")

IF(NOT DEFINED CMAKE_COLOR_MAKEFILE)
    SET(CMAKE_COLOR_MAKEFILE OFF)
ENDIF()

IF(CMAKE_SIZEOF_VOID_P EQUAL 8)
    SET(ARCH64_BITS TRUE)
ELSE()
    SET(ARCH64_BITS FALSE)
ENDIF()

IF(NOT DEFINED CMAKE_BUILD_TYPE)
    SET(CMAKE_BUILD_TYPE Debug)
ENDIF()

MESSAGE(STATUS "CMAKE_BUILD_TYPE = ${CMAKE_BUILD_TYPE}")
STRING(TOUPPER "${CMAKE_BUILD_TYPE}" CMAKE_BUILD_TYPE_UPPER)
IF(CMAKE_BUILD_TYPE_UPPER MATCHES DEBUG)
    SET(BUILD_SHARED_LIBS OFF)
ENDIF()

IF(${CMAKE_BUILD_TYPE_UPPER} STREQUAL "COVERAGE")
    MESSAGE(STATUS "Adding coverage")
    FIND_PACKAGE(codecov)
ENDIF()

IF(DEFINED INSTALL_PREFIX)
    SET(CMAKE_INSTALL_PREFIX ${INSTALL_PREFIX})
ELSE()
    SET(CMAKE_INSTALL_PREFIX ../install${CMAKE_BUILD_TYPE})
ENDIF()

################################################################################
IF(NOT DEFINED FFTW_USE_STATIC_LIBS)
    SET(FFTW_USE_STATIC_LIBS TRUE)
ENDIF()
FIND_PACKAGE(FFTW REQUIRED)
MESSAGE(STATUS "FFTW_INCLUDES  : ${FFTW_INCLUDES}")
MESSAGE(STATUS "FFTW_LIBRARIES : ${FFTW_LIBRARIES}")
MESSAGE(STATUS "FFTW_LIB       : ${FFTW_LIB}")
MESSAGE(STATUS "FFTWF_LIB      : ${FFTWF_LIB}")
MESSAGE(STATUS "FFTWL_LIB      : ${FFTWL_LIB}")

IF(DEFINED LAPACK_ROOT)
    FIND_PACKAGE(LAPACK REQUIRED
                 HINTS ${LAPACK_ROOT})
ELSE()
    FIND_PACKAGE(LAPACK REQUIRED)
ENDIF()
################################################################################

MESSAGE(STATUS "CMAKE_GENERATOR = ${CMAKE_GENERATOR}")
MESSAGE(STATUS "CMAKE_INSTALL_PREFIX = ${CMAKE_INSTALL_PREFIX}")

INCLUDE(CMakeOutputDirectories)

ENABLE_LANGUAGE(Fortran)
INCLUDE(CMakeCompilationOptionsFortran)

################################################################################
INCLUDE(CMakeSources.cmake)
################################################################################

MESSAGE(STATUS "CMAKE_CURRENT_BINARY_DIR ${CMAKE_CURRENT_BINARY_DIR}")
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_BINARY_DIR})

ADD_EXECUTABLE(HOS_NWT ${src})
TARGET_LINK_LIBRARIES(HOS_NWT lapack blas ${FFTW_LIBRARIES})

IF(${CMAKE_BUILD_TYPE_UPPER} STREQUAL "COVERAGE")
    add_coverage(HOS_NWT)
ENDIF()

INSTALL(TARGETS HOS_NWT
        RUNTIME DESTINATION ${RUNTIME_OUTPUT_DIRECTORY}
        LIBRARY DESTINATION ${LIBRARY_OUTPUT_DIRECTORY})
SET(ADDITIONAL_LIB)
LIST(APPEND ADDITIONAL_LIB blas lapack)
FOREACH(f ${ADDITIONAL_LIB})
    IF(EXISTS ${f})
        INSTALL(FILES ${f} DESTINATION ${LIBRARY_OUTPUT_DIRECTORY})
    ENDIF()
ENDFOREACH()

ENABLE_TESTING()
MACRO(ADD_A_NEW_BENCHMARK test_name)
    EXECUTE_PROCESS(COMMAND ${CMAKE_COMMAND} -E make_directory ${BENCHMARK_DIRECTORY}/${test_name}/Results)
    # SET_DIRECTORY_PROPERTIES(PROPERTIES ADDITIONAL_MAKE_CLEAN_FILES TRUE)
    ADD_TEST(NAME ${test_name}
            COMMAND HOS_NWT input_HOS-NWT.dat
            WORKING_DIRECTORY ${BENCHMARK_DIRECTORY}/${test_name})
ENDMACRO()

# ADD_A_NEW_BENCHMARK(Irreg_2D)
# ADD_A_NEW_BENCHMARK(Irreg_3D)
# ADD_A_NEW_BENCHMARK(Regular)
# ADD_A_NEW_BENCHMARK(Regular_3D)
# ADD_A_NEW_BENCHMARK(Sloshing)

FILE(WRITE "${CMAKE_CURRENT_BINARY_DIR}/input_HOS-NWT.dat"
"----------------------- Tank dimensions -----------------------
Length(m) of the wavetank    :: xlen             :: 120.0d0
Beam(m) of the wavetank      :: ylen             :: 100.d0
Depth(m) of the wavetank     :: h                :: 5.d0
------------------------ Select case --------------------------
Choice of computed case      :: icase            :: 2
------------------ Sloshing case (icase = 1) ------------------
Number of the natural mode   :: islosh           :: 3
Amplitude (m) of the mode    :: aslosh           :: 0.2d0
--------------- Monochromatic case (icase = 2) ----------------
Amplitude (m) of wave train  :: amp_mono         :: 0.1d0
Frequency (Hz) of wave train :: nu_mono          :: 0.5588d0
Angle (deg) from x-axis      :: theta_mono       :: 0.d0
Phasis (rad) of wave train   :: ph_mono          :: 0.d0
Directional wmk type         :: ibat             :: 2
Wave target distance (ibat=3):: xd_mono          :: 18.d0
-------------------- File case (icase = 3) --------------------
Name of frequency file       :: file_name        :: Hs4_Tp10
Frequency cut-off            :: i_cut            :: 0
Low cut_off frequency (Hz)   :: nuc_low          :: 0.3d0
High cut_off frequency (Hz)  :: nuc_high         :: 3.d0
------------------ Irregular wave (icase=4) -------------------
Significant wave height (m)  :: Hs               :: 0.05d0
Peak period (s)              :: Tp               :: 1.d0
Shape factor (Jonswap)       :: gamma            :: 3.3d0
Seed number for random numb. :: iseed            :: 513
-------------------- Wavemaker definition ---------------------
Nonlinear order of wavemaker :: i_wmk            :: 2
Type (1: piston, 2: hinged)  :: igeom            :: 2
Rotation axis distance       :: d_hinge          :: 1.d0
Time ramp                    :: iramp            :: 0
Time ramp duration (s)       :: Tramp            :: 5.0d0
----------------------- Numerical beach -----------------------
Absorption numerical beach   :: iabsnb           :: 1
Beginning front num. beach   :: xabsf            :: 0.8d0
Absorption strength front    :: coeffabsf        :: 1.d0
------------- Elevation/Velocity-Pressure probes --------------
Use of probes                :: iprobes          :: 1
Filename of probe positions  :: pro_file         :: prob.inp
---------------------- Time-integration -----------------------
Duration of the simulation   :: T_stop           :: 1.0d0
Time tolerance: RK 4(5)      :: toler            :: 1.e-3
Output frequency             :: f_out            :: 1.d0
--------------------------- Output ----------------------------
Output: 1-dim. ; 0-nondim.   :: idim             :: 1
free surface plot            :: i3d              :: 0
modes plot                   :: imodes           :: 0
wavemaker motion plot        :: iwmk             :: 0
Swense output 1='yes',0='no' :: i_sw             :: 0
")

FILE(MAKE_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/Results")
ADD_TEST(NAME test_Irreg_2D
         COMMAND HOS_NWT input_HOS-NWT.dat
         WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})

INCLUDE(CMakePack)
