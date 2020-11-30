/*
 * SIMPLE cmake learning project
 *
 * CMAKE COMMAND
 * -------------
 * * Build:
 *      "MSVS"
 *      cmake -S . -B ./build
 *
 *      "MINGW"
 *      cmake -S . -B ./build -G "MinGW Makefiles"
 *
 *      Recompile:
 *          cmake -C ./build
 *
 *  * Compile:
 *      cmake --build ./build
 *
 *  * Run:
 *      ./build/Debug/c_cmake.exe
 * */

#include <stdio.h>
#include <stdlib.h>

#include "include/foo.h"

int main(void)
{
    printf("hello world\n");
    printf("NUMBER: %d\n", foo_func());

    return EXIT_SUCCESS;
}
