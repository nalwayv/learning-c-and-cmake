cmake_minimum_required(VERSION 3.19)
project(app VERSION 1.0.0)

# c version
set(CMAKE_C_STANDARD 11)
set(CMAKE_C_STANDARD_REQUIRED 11)

# main exe
add_executable(c_cmake main.c)

# include directories
add_subdirectory(glad)
add_subdirectory(src)
add_subdirectory(external/glfw )

target_include_directories(c_cmake
    PUBLIC src
	PUBLIC glad
    PUBLIC external/glfw/include
)

# link directories
target_link_directories(c_cmake
    PRIVATE src
	PRIVATE glad
    PRIVATE external/glfw/src
)

target_link_libraries(c_cmake glad src glfw)
