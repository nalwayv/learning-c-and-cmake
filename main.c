/*
 * SIMPLE cmake/glfw oject
 *
 * CMAKE COMMAND
 * -------------
 * * Build:
 *      "MSVS": cmake -S . -B ./build
 *      "MINGW": cmake -S . -B ./build -G "MinGW Makefiles"
 *
 *      *Recompile:
 *          cmake -C ./build
 *
 *      *CMAKE_COMMANDS:
 *          "-DCMAKE_EXPORT_COMPILE_COMMANDS=1": export compile_commands.json
 *
 *  * Compile:
 *      cmake --build ./build
 *
 *  Example:
 *      cmake -DGLFW_BUILD_DOCS=OFF -DCMAKE_EXPORT_COMPILE_COMMANDS=1 -S . -B ./build/ -G "MinGW Makefiles"
 * */

#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <assert.h>
#include <string.h>
#include <time.h>
#include <stdarg.h>
#include <math.h>

#include "glad/glad.h"
#include <GLFW/glfw3.h>

#include "src/vec.h"
#include "src/matrix.h"

// #define DEBUG_MODE

#define GL_LOG_FILE "./gl.log"
#define FRAG_FILE "./shader.frag"
#define VERT_FILE "./shader.vert"
#define not !
#define and &&
#define or ||
#define internal static
#define global_var static
#define local_persist static

const GLint WIDTH = 640;
const GLint HEIGHT = 480;
const char* NAME = "GLFW CMAKE";

/* SHADERS */

/* *
 * layout(location=): where in vbo to get data
 * in: get from previous stage of render pipeline
 * out: send info to next stage of render pipeline
 * uniform: global with data passed in from cpu
 * */

/* --- */

/* GLOBALS */

global_var GLfloat points[] = {
//   x     y    z
    -0.5, -0.5, 0.0, // px
     0.5, -0.5, 0.0, // py
     0.0,  0.5, 0.0  // pz
};

global_var GLfloat colour[] = {
//  r     g     b
    1.0f, 0.0f, 0.0f, // px
    0.0f, 1.0f, 0.0f, // py 
    0.0f, 0.0f, 1.0f  // pz
};

/*
global_var const GLfloat points_colour[] = {
//   x      y     z      r     g     b
    -0.5f, -0.5f, 0.0f,  1.0f, 0.0f, 0.0f, // px
     0.5f, -0.5f, 0.0f,  0.0f, 1.0f, 0.0f, // py
     0.0f,  0.5f, 0.0f,  0.0f, 0.0f, 1.0f, // pz
};
*/

// transform on x axis
global_var GLfloat transform_4X4[] = {
//  x     y     z     t
    1.0f, 0.0f, 0.0f, 0.0f, // vx
    0.0f, 1.0f, 0.0f, 0.0f, // vy
    0.0f, 0.0f, 1.0f, 0.0f, // vz
    0.5f, 0.0f, 0.0f, 1.0f  // transform 1
};

/* --- */

/* *
 * clamp value between min and max values
 *
 * @param value target.
 * @param min value target.
 * @param max value target.
 * @return value between min and max.
 * */
internal int
int_clamp(int value, int min_value, int max_value)
{
    return fmax(min_value, fmin(max_value, value));
}

/* *
 * clear current log file info and start again with current time
 * */
internal bool
restart_gl_log(void)
{
    FILE* fp = NULL;
    errno_t err = fopen_s(&fp, GL_LOG_FILE, "w");
    if(err)
    {
        fprintf_s(stderr, "ERROR: could not open file %s\n", GL_LOG_FILE);
        return false;
    }
   
    time_t now = time(NULL);
    char date_buffer[128];
    ctime_s(date_buffer, sizeof(*date_buffer) * 128, &now); 
    
    fprintf_s(fp, "GL_LOG: time: %s\n", date_buffer);
    
    fclose(fp);
    
    return true;
}

/* *
 * Log to gl.log file
 *
 * @param message a string.
 * @return true if message was appanded to file.
 * */
internal bool
gl_log(const char* message, ...)
{
    FILE* fp = NULL;
    errno_t err = fopen_s(&fp, GL_LOG_FILE, "a");
    if(err)
    {
        fprintf_s(stderr, "ERROR: could not open file %s\n", GL_LOG_FILE);
        return false;
    }

    va_list list;
    va_start(list, message);
    fprintf_s(fp, message, va_arg(list, const char*));
    va_end(list);
    
    fclose(fp);

    return true;
}

/* *
 * log error to gl.log
 *
 * @param message string.
 * @return true is message was writen to file.
 * */
internal bool
gl_log_error(const char* message, ...)
{
    FILE *fp = NULL;
    errno_t err = fopen_s(&fp, GL_LOG_FILE, "a");
    if(err)
    {
        fprintf_s(stderr, "ERROR: could not open file %s\n", GL_LOG_FILE);
        return false;
    }

    va_list list;
    va_start(list, message);
    fprintf_s(fp, message, va_arg(list, const char*));
    va_end(list);

    fclose(fp);
    
    return true;
}

/* *
 * print shader info log to conole
 *
 * @param idx of shader program.
 * */
internal void
print_shader_infolog(GLuint idx)
{
    int len = 0;
    char log[1024];
    glad_glGetShaderInfoLog(idx, 1024, &len, log);
    fprintf_s(stderr, "Shader Info Log:\nidx: %u\nlog:%s\n", idx, log);
}

/* *
 * print program info log to console
 *
 * @param program id generated by glCreateProgram.
 * */
internal void
print_program_infolog(GLuint program)
{
    int len = 0;
    char log[1024];
    glad_glGetProgramInfoLog(program, 1024, &len, log);
    fprintf_s(stderr, "Program Info Log:\nidx: %u \nlog: %s\n", program, log);
}

/* *
 * compile shader 
 *
 * @param shader id number from glCreateShader.
 * @param shader string source.
 * @param if it should check compile status for errors.
 * @return true if no errors found.
 * */
internal bool 
compile_shader(GLuint shader, const char* source, bool check_compile_status)
{
    glad_glShaderSource(shader, 1, &source, NULL);
    glad_glCompileShader(shader);

    if(check_compile_status)
    {
        GLint check = -1;
        glad_glGetShaderiv(shader, GL_COMPILE_STATUS, &check);   
        if(check != GL_TRUE)
        {
            fprintf(stderr, "ERROR: failed to create shader\n");
            print_shader_infolog(shader);
            return false;
        }
    }

    return true;
}

/* *
 * link shaders to shader program
 *
 * @param program id generated from glCreateProgram.
 * @param vert shader id from glCreateShader.
 * @param frag shader id from glCreateShader.
 * @param if it should check link status for errors.
 * @return true if no errors found.
 * */
internal bool
link_shaders(GLuint program, GLuint vert, GLuint frag, bool check_link_status)
{
    glad_glAttachShader(program, frag);
    glad_glAttachShader(program, vert);
    glad_glLinkProgram(program);
    
    if(check_link_status)
    {
        int check = -1;
        glad_glGetProgramiv(program, GL_LINK_STATUS, &check);
        if(check != GL_TRUE)
        {
            fprintf(stderr, "ERROR: failed to link shader program\n");
            print_program_infolog(program);
            return false;
        }
    }

    return true;
}

/* *
 * get text from file and copy into buffer
 *
 * @param file_path of text file.
 * @param buffer to store text.
 * @return true if file was successfully read.
 * */
internal bool 
get_text_from_file(const char* file_path, char *buffer)
{
    FILE *fp = NULL;
    errno_t err = fopen_s(&fp, file_path, "r");
    if(err)
    {
        fprintf_s(stderr, "ERROR: failed to open file %s\n", file_path);
        return false;
    }

    char c;
    int idx = 0;
    while((c = fgetc(fp)) != EOF)
    {
        buffer[idx] = c;
        ++idx;
    }
    buffer[idx] = '\0';
 
    fclose(fp);

    return true;
}

/*
 * create shaders vert and frag from file and link to shader program
 *
 * @param program id generated from glCreateProgram.
 * @param v_file path for vertex shader data.
 * @param f_file path for fragment shader data.
 * @return true if shaders were created and linked to program.
 * */
internal bool
create_shaders_and_link_to_program(GLuint program, const char* v_file, const char* f_file)
{
    char vert_buffer[512];
    if(not get_text_from_file(v_file, vert_buffer))
    {
        return false;
    }

    GLuint v_shader = glad_glCreateShader(GL_VERTEX_SHADER);
    if(not compile_shader(v_shader, vert_buffer, true))
    {
        return false;
    }

    char frag_buffer[512];
    if(not get_text_from_file(f_file, frag_buffer))
    {
        return false;
    }

    GLuint f_shader = glad_glCreateShader(GL_FRAGMENT_SHADER);
    if(not compile_shader(f_shader, frag_buffer, true))
    {
        return false;
    }

    if(not link_shaders(program, v_shader, f_shader, true))
    {
        return false;
    } 

    // CLEAN UP
    glad_glDeleteShader(v_shader);
    glad_glDeleteShader(f_shader);
    
    return true;
}

/* *
 * print title to console
 * */
internal void
title(void)
{
    printf("===========================\n");
    printf(" CMAKE OPEN_GL\n");
    printf("===========================\n");
}


internal void
check_matrix(void)
{
    struct mat4 m4 = mat4_new(1,  2, 3, 4,
                              5,  6, 7, 8,
                              9, 10,11,12,
                              13,14,15,16);
    struct mat4 m4_inv = mat4_inverse(&m4);
    mat4_print(&m4_inv);
}

/* --- */

/* CALLBACKS */

internal void
error_callback(int error, const char* description)
{
    gl_log_error("ERROR: id: %i, msg: %s\n", error, description);
}

internal void
key_callback(GLFWwindow *window, int key, int scancode, int action, int modes)
{
    if(key == GLFW_KEY_Q and action == GLFW_PRESS)
    {
        glfwSetWindowShouldClose(window, GLFW_TRUE);
    }
}

/* --- */

int main(void)
{
    title();
    
    check_matrix();

    assert(restart_gl_log());
    gl_log("GLFW START: version %s\n", glfwGetVersionString());
    glfwSetErrorCallback(error_callback);


    if(not glfwInit())
    {
        gl_log_error("glfw init failed to start\n");
        return EXIT_FAILURE;
    }

    glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 3);
    glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 3);
    
    GLFWwindow *window = glfwCreateWindow(WIDTH, HEIGHT, NAME, NULL, NULL);
    if(not window)
    {
        gl_log_error("glfw window failed to init\n");
        return EXIT_FAILURE;
    }
   
    glfwSetKeyCallback(window, key_callback);
    glfwMakeContextCurrent(window);
  
    /* GLAD */
    if(not gladLoadGL())
    {
        gl_log_error("glad failed to start\n");
        return EXIT_FAILURE;
    }

    const GLubyte *renderer = glGetString(GL_RENDERER);
    const GLubyte *version = glGetString(GL_VERSION);

    printf("RENDERER:: %s\nVERSION:: %s\n", renderer, version);

    glad_glEnable(GL_DEPTH_TEST);
    glad_glDepthFunc(GL_LESS);

    /* --- */
    
    // VBO (vertex buffer object) array of data
    //  position
    GLuint points_vbo = 0;
    glad_glGenBuffers(1, &points_vbo);
    glad_glBindBuffer(GL_ARRAY_BUFFER, points_vbo);
    glad_glBufferData(GL_ARRAY_BUFFER, sizeof(points), points, GL_STATIC_DRAW);

    //  colour
    GLuint colour_vbo = 0;
    glad_glGenBuffers(1, &colour_vbo);
    glad_glBindBuffer(GL_ARRAY_BUFFER, colour_vbo);
    glad_glBufferData(GL_ARRAY_BUFFER, sizeof(colour), colour, GL_STATIC_DRAW);
    
    // VAO
    GLuint vao = 0;
    glad_glGenVertexArrays(1, &vao);
    glad_glBindVertexArray(vao);
    
    //  vao position
    glad_glEnableVertexAttribArray(0);
    glad_glBindBuffer(GL_ARRAY_BUFFER, points_vbo);
    glad_glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 0, NULL);

    //  vao colour
    glad_glEnableVertexAttribArray(1);
    glad_glBindBuffer(GL_ARRAY_BUFFER, colour_vbo);
    glad_glVertexAttribPointer(1, 3, GL_FLOAT, GL_FALSE, 0, NULL);
   
    
    // SHADER PROGRAM
    GLuint program = glad_glCreateProgram();
    if(not create_shaders_and_link_to_program(program, VERT_FILE, FRAG_FILE))
    {
        return EXIT_FAILURE;
    }

    /* --- */

    // speed
    // double speed = 0.2;
    // double last_pos = 0.0;

    // fps
    double previous = glfwGetTime();
    int counter = 0;
 
    while(!glfwWindowShouldClose(window))
    {
        double current = glfwGetTime();
        double elapsed = current - previous;
        if(elapsed > 0.25)
        {
            previous = current;
            
            double fps = (double)counter / elapsed;
            char buffer[128];
            sprintf_s(buffer, 128, "FPS: %.2f", fps);
            glfwSetWindowTitle(window, buffer);

            counter = 0;
        }
        ++counter;

        // move by updating transform matrix
        // speed = fabs(last_pos) > 1.0 ? -speed : speed;
        // transform_4X4[12] = elapsed * speed + last_pos;
        // last_pos = transform_4X4[12];


        glad_glClearColor(0.1, 0.1, 0.1, 1.0);
        glad_glUseProgram(program);
        glad_glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
        glad_glViewport(0, 0, WIDTH, HEIGHT);

        /* DRAW */
        
        // int transform_local = glad_glGetUniformLocation(program, "transform");
        // glad_glUniformMatrix4fv(transform_local, 1, GL_FALSE, transform_4X4);
        
        glad_glBindVertexArray(vao);
        glad_glDrawArrays(GL_TRIANGLES, 0, 3);
       
        // wireframe mode 
        // glad_glPolygonMode(GL_FRONT, GL_LINE);

        /* --- */

        glfwSwapBuffers(window);
        glfwPollEvents();
    }

    // CLEAN UP WINDOW
    glfwDestroyWindow(window);
    glfwTerminate();

    return EXIT_SUCCESS;
}
