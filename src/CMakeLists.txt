cmake_minimum_required (VERSION 2.6)

set (CMAKE_CXX_FLAGS_RELEASE "-O3")
set (CMAKE_CXX_FLAGS_DEBUG  "-g -O2")
set (CMAKE_CXX_COMPILER ${MPI_COMPILER})

project(LW)
set (LW_VERSION_MAJOR 0)
set (LW_VERSION_MINOR 1)
set (PACKAGE \"lw\")
set (PACKAGE_NAME \"LW\")
set (PACKAGE_STRING \"LW\ 0.1\")
set (PACKAGE_TARNAME \"LW\")
set (PACKAGE_VERSION \"0.1\")
set (VERSION \"0.1\")

set(CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/CMakeModules")

find_package (IPPL REQUIRED)
find_package (GSL REQUIRED)
find_package (TCLAP REQUIRED)
find_package (Boost REQUIRED)

set (CMAKE_CXX_FLAGS
     "-O3 -std=c++0x -Wall -DUSEH5FEDV2 -DIPPL_MPI -DMPICH_SKIP_MPICXX -DIPPL_DEBUG -DIPPL_DONT_POOL -DIPPL_USE_XDIV_RNG -DIPPL_LINUX -DIPPL_NO_STRINGSTREAM -DPETE_BITWISE_COPY -DIPPL_HAS_TEMPLATED_COMPLEX -DIPPL_USE_STANDARD_HEADERS -DIPPL_USE_PARTIAL_SPECIALIZATION -DIPPL_STDSTL -DIPPL_LONGLONG -DPARALLEL_IO -DNOCTAssert -Drestrict=__restrict__ -DH5_USE_16_API ${CMAKE_CXX_FLAGS}")

set (LW_SRCS
    PlotBinaryVTK.cpp
    BinaryVTKFile.cpp
    utils.cpp
    Physics.cpp
    Communicator.cpp
)

include_directories(
    ${CMAKE_CURRENT_SOURCE_PATH}
    ${IPPL_INCLUDE_DIR}
    ${TCLAP_INCLUDE_DIR}
    ${BOOST_INCLUDE_DIR}
)

set (LW_LIBS
    m
    ${IPPL_LIBRARY}
    ${BOOST_LIBRARY}
    ${GSL_LIBRARY}
    ${GSL_LIBRARY_CBLAS}
#    gfortran
    z
)

add_executable( lw ${LW_SRCS} main.cpp )
target_link_libraries( lw ${LW_LIBS})
