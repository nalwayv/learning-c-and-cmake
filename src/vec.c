#include "vec.h"

#include "assert.h"
#include "math.h"
#include <math.h>

/*
 * HELPERS
 *
static void
matrix_multiply_by_vector(float *dest, const float *vec, const float *mat, int n)
{
    for(int x = 0; x < n; ++x)
    {
        float v = 0.0f;
        for(int y = 0; y < n; ++y)
        {
            v +=  mat[n * x + y] * vec[x];
        }
        dest[x] = v; 
    }
}
*/

/* --- VEC 2 --- */

vec2 vec2_new(float x, float y)
{
    vec2 v2 = (vec2){.x=x, .y=y};
    return v2;
}

vec2 vec2_zero()
{
    return vec2_new(0.0f, 0.0f);
}

vec2 vec2_add(const vec2 *a, const vec2 *b)
{
    return vec2_new(a->x + b->x, a->y + b->y);
}

vec2 vec2_sub(const vec2 *a, const vec2 *b)
{
    return vec2_new(a->x - b->x, a->y - b->y);
}

vec2 vec2_mul(const vec2 *v2, float by)
{
    return vec2_new(v2->x * by, v2->y * by);
}

vec2 vec2_div(const vec2 *v2, float by)
{
    assert(by != 0.0);
    return vec2_new(v2->x / by, v2->y / by);
}

float vec2_dot(const vec2 *a, const vec2 *b)
{
    return (a->x * b->x) + (a->y * b->y);
}

float vec2_len(const vec2 *v2)
{
    float dot = vec2_dot(v2, v2);
    assert(dot != 0.0f);
    return sqrtf(dot);
}

vec2 vec2_normal(const vec2 *v2)
{
    float dot = vec2_dot(v2,v2);
    assert(dot != 0.0f);
    float by = 1.0f / sqrtf(dot);
    vec2 result;
    result.x = v2->x * by;
    result.y = v2->y * by;
    return result;
}

vec2 vec2_mat2(const vec2 *v2, const float *mat_arr)
{
    float x = mat_arr[0]*v2->x + mat_arr[1]*v2->y;
    float y = mat_arr[2]*v2->x + mat_arr[3]*v2->y;
    return vec2_new(x, y);
}

/* --- VEC 3 --- */

vec3 vec3_new(float x, float y, float z)
{
    vec3 v3 = (vec3){.x=x, .y=y, .z=z};
    return v3;
}

vec3 vec3_zero()
{
    return vec3_new(0.0f, 0.0f, 0.0f);
}

vec3 vec3_add(const vec3 *a, const vec3 *b)
{
    return vec3_new(a->x + b->x, a->y + b->y, a->z + b->z);
}

vec3 vec3_sub(const vec3 *a, const vec3 *b)
{
    return vec3_new(a->x - b->x, a->y - b->y, a->z - b->z);
}

vec3 vec3_mul(const vec3 *v3, float by)
{
    return vec3_new(v3->x * by, v3->y * by, v3->z * by);
}

vec3 vec3_div(const vec3 *v3, float by)
{
    assert(by != 0.0);
    return vec3_new(v3->x / by, v3->y / by, v3->z / by);
}

float vec3_dot(const vec3 *a, const vec3 *b)
{
    return (a->x * b->x) + (a->y * b->y) + (a->z  * b->z);
}

vec3 vec3_cross(const vec3 *a, const vec3 *b)
{
    float x = (a->y * b->z) - (a->z * b->y);
    float y = (a->z * b->x) - (a->x * b->z);
    float z = (a->x * b->y) - (a->y * b->x);
    return vec3_new(x, y, z);
}

float vec3_len(const vec3 *v3)
{
    float dot = vec3_dot(v3, v3);
    assert(dot != 0.0f);
    return sqrtf(dot);
}

vec3 vec3_normal(const vec3 *v3)
{
    float dot = vec3_dot(v3, v3);
    assert(dot != 0.0f);
    float by = 1.0f / sqrtf(dot);

    vec3 result;
    result.x = v3->x * by;
    result.y = v3->y * by;
    result.z = v3->z * by;
    return result;
}

vec3 vec3_mat3(const vec3 *v3, const float* mat_arr)
{
    float x = mat_arr[0]*v3->x + mat_arr[1]*v3->y + mat_arr[2]*v3->z;
    float y = mat_arr[3]*v3->x + mat_arr[4]*v3->y + mat_arr[5]*v3->z;
    float z = mat_arr[6]*v3->x + mat_arr[7]*v3->y + mat_arr[8]*v3->z;
    return vec3_new(x,y,z);
}

/* --- VEC 4 --- */


vec4 vec4_new(float x, float y, float z, float w)
{
    vec4 v4 = (vec4){.x=x, .y=y, .z=z, .w=w};
    return v4;
}

vec4 vec4_zero()
{
    return vec4_new(0.0f, 0.0f, 0.0f, 0.0f);
}

vec4 vec4_add(const vec4 *a, const vec4 *b)
{
    return vec4_new(a->x + b->x, a->y + b->y, a->z + b->z , a->w + b->w);
}

vec4 vec4_sub(const vec4 *a, const vec4 *b)
{
    return vec4_new(a->x - b->x, a->y - b->y, a->z - b->z , a->w - b->w);
}

vec4 vec4_mul(const vec4 *v4, float by)
{
    return vec4_new(v4->x * by, v4->y * by, v4->z * by, v4->w * by);
}

vec4 vec4_div(const vec4 *v4, float by)
{
    assert(by != 0.0f);
    return vec4_new(v4->x / by, v4->y / by, v4->z / by, v4->w / by);
}

float vec4_dot(const vec4 *a, const vec4 *b)
{
    return (a->x * b->x) + (a->y * b->y) + (a->z * b->z) + (a->w * b->w);
}

float vec4_len(const vec4 *v4)
{
    float dot = vec4_dot(v4,v4);
    assert(dot != 0.0f);
    return sqrtf(dot);
}

vec4 vec4_normal(const vec4 *v4)
{
    float dot = vec4_dot(v4, v4);
    assert(dot != 0.0f);
    float by = 1.0f / sqrtf(dot);

    vec4 result;
    result.x = v4->x * by;
    result.y = v4->y * by;
    result.z = v4->z * by;
    result.w = v4->w * by;
    return result;
}

vec4 vec4_mat4(const vec4 *v4, const float *mat_arr)
{
    float x = mat_arr[0]*v4->x + mat_arr[1]*v4->y + mat_arr[2]*v4->z + mat_arr[3]*v4->w;
    float y = mat_arr[4]*v4->x + mat_arr[5]*v4->y + mat_arr[6]*v4->z + mat_arr[7]*v4->w;
    float z = mat_arr[8]*v4->x + mat_arr[9]*v4->y + mat_arr[10]*v4->z + mat_arr[11]*v4->w;
    float w = mat_arr[12]*v4->x + mat_arr[13]*v4->y + mat_arr[14]*v4->z + mat_arr[15]*v4->w;
    return vec4_new(x, y, z, w); 
}
