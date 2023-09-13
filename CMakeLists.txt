cmake_minimum_required(VERSION 3.0.0)
project(delta VERSION 0.1.0 LANGUAGES C CXX)
file(GLOB SOURCES "src/*.cpp")
add_executable(NJL2 main.cpp ${SOURCES})
target_include_directories(NJL2 PRIVATE include)
set(CPACK_PROJECT_NAME ${PROJECT_NAME})
set(CPACK_PROJECT_VERSION ${PROJECT_VERSION})
include(CPack)




