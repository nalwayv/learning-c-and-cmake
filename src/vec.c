#include "vec.h"

#include "assert.h"
#include "math.h"


/* --- VEC 2 --- */


struct vec2 vec2_new(float x, float y)
{
    struct vec2 v2 = (struct vec2){.x=x, .y=y};
    return v2;
}

struct vec2 vec2_add(const struct vec2 *a, const struct vec2 *b)
{
    return vec2_new(a->x + b->x, a->y + b->y);
}

struct vec2 vec2_sub(const struct vec2 *a, const struct vec2 *b)
{
    return vec2_new(a->x - b->x, a->y - b->y);
}

struct vec2 vec2_mul(const struct vec2 *v2, float by)
{
    return vec2_new(v2->x * by, v2->y * by);
}

struct vec2 vec2_div(const struct vec2 *v2, float by)
{
    assert(by != 0.0);
    return vec2_new(v2->x / by, v2->y / by);
}

float vec2_dot(const struct vec2 *a, const struct vec2 *b)
{
    return (a->x * b->x) + (a->y * b->y);
}

float vec2_len(const struct vec2 *v2)
{
    float dot = vec2_dot(v2, v2);
    assert(dot != 0.0f);
    return sqrtf(dot);
}

struct vec2 vec2_normal(const struct vec2 *v2)
{
    return vec2_div(v2, vec2_len(v2));
}

/* --- VEC 3 --- */

struct vec3 vec3_new(float x, float y, float z)
{
    struct vec3 v3 = (struct vec3){.x=x, .y=y, .z=z};
    return v3;
}

struct vec3 vec3_add(const struct vec3 *a, const struct vec3 *b)
{
    return vec3_new(a->x + b->x, a->y + b->y, a->z + b->z);
}

struct vec3 vec3_sub(const struct vec3 *a, const struct vec3 *b)
{
    return vec3_new(a->x - b->x, a->y - b->y, a->z - b->z);
}

struct vec3 vec3_mul(const struct vec3 *v3, float by)
{
    return vec3_new(v3->x * by, v3->y * by, v3->z * by);
}

struct vec3 vec3_div(const struct vec3 *v3, float by)
{
    assert(by != 0.0);
    return vec3_new(v3->x / by, v3->y / by, v3->z / by);
}

float vec3_dot(const struct vec3 *a, const struct vec3 *b)
{
    return (a->x * b->x) + (a->y * b->y) + (a->z  * b->z);
}

struct vec3 vec3_cross(const struct vec3 *a, const struct vec3 *b)
{
    float x = (a->y * b->z) - (a->z * b->y);
    float y = (a->z * b->x) - (a->x * b->z);
    float z = (a->x * b->y) - (a->y * b->x);
    return vec3_new(x, y, z);
}

float vec3_len(const struct vec3 *v3)
{
    float dot = vec3_dot(v3, v3);
    assert(dot != 0.0f);
    return sqrtf(dot);
}

struct vec3 vec3_normal(const struct vec3 *v3)
{
    return vec3_div(v3, vec3_len(v3));
}


/* --- VEC 4 --- */


struct vec4 vec4_new(float x, float y, float z, float a)
{
    struct vec4 v4 = (struct vec4){.x=x, .y=y, .z=z, .a=a};
    return v4;
}
