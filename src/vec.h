#ifndef _VEC_
#define _VEC_

struct vec2
{
    struct {
        float x;
        float y;
    };

    float arr[2];
};
typedef struct vec2 vec2;


struct vec3
{
    struct {
        float x;
        float y;
        float z;
    };

    float arr[3];
};
typedef struct vec3 vec3;


struct vec4
{
    struct {
        float x;
        float y;
        float z;
        float w;
    };

    float arr[4];
};
typedef struct vec4 vec4;

/* --- VEC 2 --- */

/* *
 * create a new vec2
 * */
vec2 vec2_new(float x, float y);


/* *
 * create a new zero vec2
 * */
vec2 vec2_zero();

/* *
 * add vec2 b to vec2 a
 * */
vec2 vec2_add(const vec2 *a, const vec2 *b);

/* *
 * sub vec2 b from vec2 a
 * */
vec2 vec2_sub(const vec2 *a, const vec2 *b);

/* *
 * scale vec2 by value
 * */
vec2 vec2_mul(const vec2 *v2, float by);

/* *
 * divide vec2 by value
 * */
vec2 vec2_div(const vec2 *v2, float by);

/* *
 * get dot product between vec2 a and vec2 b
 * */
float vec2_dot(const vec2 *a, const vec2 *b);

/* *
 * get length of vec2
 * */
float vec2_len(const vec2 *v2);

/* *
 * get vec2 normal
 * */
vec2 vec2_normal(const vec2 *v2);

/* *
 * multiply matrix 2x2 by vec2
 * */
vec2 vec2_mat2(const vec2 *v2, const float *mat_arr);


/* --- VEC 3 --- */


/* *
 * create a new vec3
 * */
vec3 vec3_new(float x, float y, float z);

/* *
 * create a new zero vec3
 * */
vec3 vec3_zero();

/* *
 * add vec3 b to vec3 a
 * */
vec3 vec3_add(const vec3 *a, const vec3 *b);

/* *
 * sub vec3 b from vec3 a
 * */
vec3 vec3_sub(const vec3 *a, const vec3 *b);

/* *
 *  scale vec3 a by value
 * */
vec3 vec3_mul(const vec3 *v3, float by);

/* *
 *  divide vec3 a by value
 * */
vec3 vec3_div(const vec3 *v3, float by);

/* *
 * get dot product between vec3 a and vec3 b
 * */
float vec3_dot(const vec3 *a, const vec3 *b);

/* *
 * get cross product between vec3 a and vec3 b
 * */
vec3 vec3_cross(const vec3 *a, const vec3 *b);

/* *
 * get length of vec3 a
 * */
float vec3_len(const vec3 *v3);

/* *
 * get normal vec3
 * */
vec3 vec3_normal(const vec3 *v3);

/* *
 * multiply matrix 3x3 by vec3;
 * */
vec3 vec3_mat3(const vec3 *v3, const float* mat_arr);

/* --- VEC 4 --- */


/* *
 * create new vec4
 * */
vec4 vec4_new(float x, float y, float z, float w);

/* *
 * create new zero vec4
 * */
vec4 vec4_zero();

/* *
 * add vec4 b to vec4 a
 * */
vec4 vec4_add(const vec4 *a, const vec4 *b);

/* *
 * sub vec4 b from vec4 a
 * */
vec4 vec4_sub(const vec4 *a, const vec4 *b);

/* *
 *  scale vec4 a by value
 * */
vec4 vec4_mul(const vec4 *v4, float by);

/* *
 *  divide vec4 a by value
 * */
vec4 vec4_div(const vec4 *v4, float by);

/* *
 * get dot product between vec4 a and vec4 b
 * */
float vec4_dot(const vec4 *a, const vec4 *b);

/* *
 * get length of vec4 a
 * */
float vec4_len(const vec4 *v4);

/* *
 * get normal vec4
 * */
vec4 vec4_normal(const vec4 *v4);


/* *
 * multiply matrix 4x4 by vec4
 * */
vec4 vec4_mat4(const vec4 *v4, const float *mat_arr);

#endif /* _VEC_ */
