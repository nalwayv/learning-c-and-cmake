cmake_minimum_required(VERSION 3.19)
project(app VERSION 1.0.0)

# c version
set(CMAKE_C_STANDARD 11)
set(CMAKE_C_STANDARD_REQUIRED 11)

# main exe
add_executable(c_cmake main.c)

# include directories
add_subdirectory(include)
target_link_libraries(c_cmake include)
target_link_directories(c_cmake PRIVATE include)

