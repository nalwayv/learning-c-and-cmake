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

struct vec3
{
    struct {
        float x;
        float y;
        float z;
    };

    float arr[3];
};

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


/* --- VEC 2 --- */

/* *
 * create new vec2
 * */
struct vec2 vec2_new(float x, float y);

/* *
 * add vec2 b to vec2 a
 * */
struct vec2 vec2_add(const struct vec2 *a, const struct vec2 *b);

/* *
 * sub vec2 b from vec2 a
 * */
struct vec2 vec2_sub(const struct vec2 *a, const struct vec2 *b);

/* *
 * scale vec2 by value
 * */
struct vec2 vec2_mul(const struct vec2 *v2, float by);

/* *
 * divide vec2 by value
 * */
struct vec2 vec2_div(const struct vec2 *v2, float by);

/* *
 * get dot product between vec2 a and vec2 b
 * */
float vec2_dot(const struct vec2 *a, const struct vec2 *b);

/* *
 * get length of vec2
 * */
float vec2_len(const struct vec2 *v2);

/* *
 * get vec2 normal
 * */
struct vec2 vec2_normal(const struct vec2 *v2);


/* --- VEC 3 --- */


/* *
 * create new vec3
 * */
struct vec3 vec3_new(float x, float y, float z);

/* *
 * add vec3 b to vec3 a
 * */
struct vec3 vec3_add(const struct vec3 *a, const struct vec3 *b);

/* *
 * sub vec3 b from vec3 a
 * */
struct vec3 vec3_sub(const struct vec3 *a, const struct vec3 *b);

/* *
 *  scale vec3 a by value
 * */
struct vec3 vec3_mul(const struct vec3 *v3, float by);

/* *
 *  divide vec3 a by value
 * */
struct vec3 vec3_div(const struct vec3 *v3, float by);

/* *
 * get dot product between vec3 a and vec3 b
 * */
float vec3_dot(const struct vec3 *a, const struct vec3 *b);

/* *
 * get cross product between vec3 a and vec3 b
 * */
struct vec3 vec3_cross(const struct vec3 *a, const struct vec3 *b);

/* *
 * get length of vec3 a
 * */
float vec3_len(const struct vec3 *v3);

/* *
 * get normal vec3 
 * */
struct vec3 vec3_normal(const struct vec3 *v3);


/* --- VEC 4 --- */


/* *
 * create new vec4
 * */
struct vec4 vec4_new(float x, float y, float z, float a);

#endif /* _VEC_ */
