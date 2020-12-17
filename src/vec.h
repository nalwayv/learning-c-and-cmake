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
        float a;
    };

    float arr[4];
};
typedef struct vec4 vec4;

/* --- VEC 2 --- */

/* *
 * create new vec2
 * */
vec2 vec2_new(float x, float y);

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


/* --- VEC 3 --- */


/* *
 * create new vec3
 * */
vec3 vec3_new(float x, float y, float z);

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


/* --- VEC 4 --- */


/* *
 * create new vec4
 * */
vec4 vec4_new(float x, float y, float z, float a);

#endif /* _VEC_ */
