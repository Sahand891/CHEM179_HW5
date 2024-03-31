cmake_minimum_required(VERSION 3.26)
project(HW5)

set(CMAKE_CXX_STANDARD 17)

include_directories(.)

add_executable(HW5
        Fock_Matrix.h
        Fock_Matrix.cpp
        Diag_Conv_CNDO2.cpp
        Diag_Conv_CNDO2.h
        Data_Abstraction.cpp
        Data_Abstraction.h
        Utils.cpp
        Utils.h
        main.cpp
        Position_Gradient.cpp
        Position_Gradient.h)

include_directories(/Users/sahandadibnia/armadillo-12.6.7/include)

# Find LAPACK
find_package(LAPACK REQUIRED)

# Add Armadillo and LAPACK libraries
find_library(ARMA_LIBRARY armadillo)
target_link_libraries(HW5 ${ARMA_LIBRARY} ${LAPACK_LIBRARIES})
