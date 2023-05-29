cmake_minimum_required(VERSION 3.19)
project(Project1)

set(CMAKE_CXX_STANDARD 14)

add_executable(Project1 main.cpp p1solvers.cpp p1solvers.h mat.cpp mat.h vect.cpp vect.h matrix_csr.cpp matrix_csr.h)