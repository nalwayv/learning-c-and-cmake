/* *
 *
 * */
#include <stdio.h>
#include <math.h>

#include "matrix.h"
#include "utils.h"


/* MAT 2 */

mat2 mat2_zero(){
    mat2 m2 = (mat2){
        ._11 = 0.0f, ._12 = 0.0f,
        ._21 = 0.0f, ._22 = 0.0f
    };
    return m2;
}

mat2 mat2_identity(){
    mat2 m2 = (mat2){
        ._11 = 1.0f, ._12 = 0.0f,
        ._21 = 0.0f, ._22 = 1.0f
    };
    return m2;
}

mat2 mat2_new(float _11, float _12,
              float _21, float _22)
{
    mat2 m2 = (mat2){
        ._11 = _11, ._12 = _12,
        ._21 = _21, ._22 = _22
    };
    return m2;
}

void
mat2_add(mat2 *dest, const mat2 *a, const mat2 *b)
{
    for(int i = 0; i < 4; ++i)
    {
        dest->arr[i] = a->arr[i] + b->arr[i];
    }
}

void
mat2_sub(mat2 *dest, const mat2 *a, const mat2 *b)
{
    for(int i = 0; i < 4; ++i)
    {
        dest->arr[i] = a->arr[i] - b->arr[i];
    }
}

void
mat2_scale(mat2 *dest, const mat2 *m2, float by)
{
    for(int i = 0; i < 4; ++i)
    {
        dest->arr[i] = m2->arr[i] * by;
    }
}

void
mat2_mul(mat2 *dest, const mat2 *a, const mat2 *b)
{

    dest->_11 = a->_11 * b->_11 + a->_12 * b->_21;
    dest->_12 = a->_11 * b->_12 + a->_12 * b->_22;
    dest->_21 = a->_21 * b->_11 + a->_22 * b->_21;
    dest->_22 = a->_21 * b->_12 + a->_22 * b->_22;
}

void
mat2_transpose(mat2 *dest, const mat2 *m2)
{
    /* *
     * 11 12 -> 11 21
     * 21 22    12 22
     * */

    dest->_11 = m2->_11;
    dest->_12 = m2->_21;
    dest->_21 = m2->_12;
    dest->_22 = m2->_22;
}

void
mat2_cofactor(mat2 *dest, const mat2 *m2)
{
    dest->_11 =  m2->_11;
    dest->_12 = -m2->_12;
    dest->_21 = -m2->_21;
    dest->_22 =  m2->_22;
}

float
mat2_determinant(const mat2 *m2)
{
    return (m2->_11 * m2->_22) - (m2->_12 * m2->_21);
}

bool
mat2_inverse(mat2 *dest, const mat2 *m2)
{
    float m2_det = mat2_determinant(m2);
    if(m2_det == 0.0f)
    {
        return false;
    }
    float det = 1.0f / m2_det;

    dest->_11 =  (m2->_22) * det;
    dest->_12 = -(m2->_12) * det;
    dest->_21 = -(m2->_21) * det;
    dest->_22 =  (m2->_11) * det;

    return true;
}

void
mat2_print(const mat2 *m2)
{
    printf("%.3f, %.3f \n", m2->_11, m2->_12);
    printf("%.3f, %.3f \n", m2->_21, m2->_22);
}

/* MAT 3 */

mat3
mat3_zero(){
    mat3 m3 = (mat3){
        ._11 = 0.0f, ._12 = 0.0f, ._13 = 0.0f,
        ._21 = 0.0f, ._22 = 0.0f, ._23 = 0.0f,
        ._31 = 0.0f, ._32 = 0.0f, ._33 = 0.0f
    };
    return m3;
}

mat3
mat3_identity(){
    mat3 m3 = (mat3){
        ._11 = 1.0f, ._12 = 0.0f, ._13 = 0.0f,
        ._21 = 0.0f, ._22 = 1.0f, ._23 = 0.0f,
        ._31 = 0.0f, ._32 = 0.0f, ._33 = 1.0f
    };
    return m3;
}

mat3
mat3_new(float _11, float _12, float _13,
         float _21, float _22, float _23,
         float _31, float _32, float _33)
{
    mat3 m3 = (mat3){
        ._11 = _11, ._12 = _12, ._13 = _13,
        ._21 = _21, ._22 = _22, ._23 = _23,
        ._31 = _31, ._32 = _32, ._33 = _33
    };
    return m3;
}

void
mat3_add(mat3 *dest, const mat3 *a, const mat3 *b)
{
    for(int i = 0; i < 9; ++i)
    {
        dest->arr[i] = a->arr[i] + b->arr[i];
    }
}

void
mat3_sub(mat3 *dest, const mat3 *a, const mat3 *b)
{
    for(int i = 0; i < 9; ++i)
    {
        dest->arr[i] = a->arr[i] - b->arr[i];
    }
}

void
mat3_scale(mat3 *dest, const mat3 *m3, float by)
{
    for(int i = 0; i < 9; ++i)
    {
        dest->arr[i] = m3->arr[i] * by;
    }
}

void
mat3_mul(mat3 *dest, const mat3 *a, const mat3 *b)
{

    dest->_11 = a->_11 * b->_11 + a->_12 * b->_21 + a->_13 * b->_31;
    dest->_12 = a->_11 * b->_12 + a->_12 * b->_22 + a->_13 * b->_32;
    dest->_13 = a->_11 * b->_13 + a->_12 * b->_23 + a->_13 * b->_33;
    dest->_21 = a->_21 * b->_11 + a->_22 * b->_21 + a->_23 * b->_31;
    dest->_22 = a->_21 * b->_12 + a->_22 * b->_22 + a->_23 * b->_32;
    dest->_23 = a->_21 * b->_13 + a->_22 * b->_23 + a->_23 * b->_33;
    dest->_31 = a->_31 * b->_11 + a->_32 * b->_21 + a->_33 * b->_31;
    dest->_32 = a->_31 * b->_12 + a->_32 * b->_22 + a->_33 * b->_32;
    dest->_33 = a->_31 * b->_13 + a->_32 * b->_23 + a->_33 * b->_33;

}

void
mat3_transpose(mat3 *dest, const struct mat3 *m3)
{
    /* *
     * 11 12 13    11 21 31
     * 21 22 23 -> 12 22 32
     * 31 32 33    13 23 33
     * */
    dest->_11 = m3->_11;
    dest->_12 = m3->_21;
    dest->_13 = m3->_31;
    dest->_21 = m3->_12;
    dest->_22 = m3->_22;
    dest->_23 = m3->_32;
    dest->_31 = m3->_13;
    dest->_32 = m3->_23;
    dest->_33 = m3->_33;
}

void
mat3_cofactor(mat3 *dest, const struct mat3 *m3)
{
    dest->_11 =  m3->_11;
    dest->_12 = -m3->_12;
    dest->_13 =  m3->_13;
    dest->_21 = -m3->_21;
    dest->_22 =  m3->_22;
    dest->_23 = -m3->_23;
    dest->_31 =  m3->_31;
    dest->_32 = -m3->_32;
    dest->_33 =  m3->_33;
}

float mat3_determinant(const struct mat3 *m3)
{
    float a = m3->_11 * (m3->_22 * m3->_33 - m3->_23 * m3->_32);
    float b = -m3->_12 * (m3->_21 * m3->_33 - m3->_23 * m3->_31);
    float c = m3->_13 * (m3->_21 * m3->_32 - m3->_22 * m3->_31);
    return a + b + c;
}

bool
mat3_inverse(mat3 *dest, const mat3 *m3)
{
    float m3_det =  mat3_determinant(m3);
    if(m3_det == 0.0f)
    {
        return false;
    }
    float det = 1.0f / m3_det;

    dest->_11 =  (m3->_22 * m3->_33 - m3->_23 * m3->_32) * det;
    dest->_12 = -(m3->_12 * m3->_33 - m3->_32 * m3->_13) * det;
    dest->_13 =  (m3->_12 * m3->_23 - m3->_22 * m3->_13) * det;
    dest->_21 = -(m3->_21 * m3->_33 - m3->_31 * m3->_23) * det;
    dest->_22 =  (m3->_11 * m3->_33 - m3->_13 * m3->_31) * det;
    dest->_23 = -(m3->_11 * m3->_23 - m3->_21 * m3->_13) * det;
    dest->_31 =  (m3->_21 * m3->_32 - m3->_31 * m3->_22) * det;
    dest->_32 = -(m3->_11 * m3->_32 - m3->_31 * m3->_12) * det;
    dest->_33 =  (m3->_11 * m3->_22 - m3->_12 * m3->_21) * det;

    return true;
}

void
mat3_rotate_x(mat3 *dest, float by)
{
    by = DEG2RAD(by);

    dest->_11 = 1.0f;
    dest->_12 = 0.0f;
    dest->_13 = 0.0f;
    dest->_21 = 0.0f;
    dest->_22 = cosf(by);
    dest->_23 = sinf(by);
    dest->_31 = 0.0f;
    dest->_32 = -sinf(by);
    dest->_33 = cosf(by);
}

void
mat3_rotate_y(mat3 *dest, float by)
{
    by = DEG2RAD(by);

    dest->_11 = cosf(by);
    dest->_12 = 0.0f;
    dest->_13 = -sinf(by);
    dest->_21 = 0.0f;
    dest->_22 = 1.0f;
    dest->_23 = 0.0f;
    dest->_31 = sinf(by);
    dest->_32 = 0.0f;
    dest->_33 = cosf(by);
}

void
mat3_rotate_z(mat3 *dest, float by)
{
    by = DEG2RAD(by);

    dest->_11 = cosf(by);
    dest->_12 = sinf(by);
    dest->_13 = 0.0f;
    dest->_21 = -sinf(by);
    dest->_22 = cosf(by);
    dest->_23 = 0.0f;
    dest->_31 = 0.0f;
    dest->_32 = 0.0f;
    dest->_33 = 1.0f;
}

void
mat3_print(const mat3 *m3)
{
    printf("%.3f, %.3f, %.3f \n", m3->_11, m3->_12, m3->_13);
    printf("%.3f, %.3f, %.3f \n", m3->_21, m3->_22, m3->_23);
    printf("%.3f, %.3f, %.3f \n", m3->_31, m3->_32, m3->_33);
}

/* MAT 4 */

mat4 mat4_zero(){
    mat4 m4 = (mat4){
        ._11 = 0.0f, ._12 = 0.0f, ._13 = 0.0f, ._14 = 0.0f,
        ._21 = 0.0f, ._22 = 0.0f, ._23 = 0.0f, ._24 = 0.0f,
        ._31 = 0.0f, ._32 = 0.0f, ._33 = 0.0f, ._34 = 0.0f,
        ._41 = 0.0f, ._42 = 0.0f, ._43 = 0.0f, ._44 = 0.0f
    };
    return m4;
}

mat4 mat4_identity(){
    mat4 m4 = (mat4){
        ._11 = 1.0f, ._12 = 0.0f, ._13 = 0.0f, ._14 = 0.0f,
        ._21 = 0.0f, ._22 = 1.0f, ._23 = 0.0f, ._24 = 0.0f,
        ._31 = 0.0f, ._32 = 0.0f, ._33 = 1.0f, ._34 = 0.0f,
        ._41 = 0.0f, ._42 = 0.0f, ._43 = 0.0f, ._44 = 1.0f
    };
    return m4;
}

mat4 mat4_new(float _11, float _12, float _13, float _14,
              float _21, float _22, float _23, float _24,
              float _31, float _32, float _33, float _34,
              float _41, float _42, float _43, float _44)
{
    mat4 m4 = (mat4){
        ._11 = _11, ._12 = _12, ._13 = _13, ._14 = _14,
        ._21 = _21, ._22 = _22, ._23 = _23, ._24 = _24,
        ._31 = _31, ._32 = _32, ._33 = _33, ._34 = _34,
        ._41 = _41, ._42 = _42, ._43 = _43, ._44 = _44
    };
    return m4;
}

void
mat4_add(mat4 *dest, const mat4 *a, const mat4 *b)
{
    for(int i = 0; i < 16; ++i)
    {
        dest->arr[i] = a->arr[i] + b->arr[i];
    }
}

void
mat4_sub(mat4 *dest, const mat4 *a, const mat4 *b)
{
    for(int i = 0; i < 16; ++i)
    {
        dest->arr[i] = a->arr[i] - b->arr[i];
    }
}

void
mat4_scale(mat4 *dest, const mat4 *a, float by)
{
    for(int i = 0; i < 16; ++i)
    {
        dest->arr[i] = a->arr[i] * by;
    }
}

void
mat4_mul(mat4 *dest, const mat4 *a, const mat4 *b)
{
    dest->_11 = a->_11 * b->_11 + a->_12 * b->_21 + a->_13 * b->_31 + a->_14 * b->_41;
    dest->_12 = a->_11 * b->_12 + a->_12 * b->_22 + a->_13 * b->_32 + a->_14 * b->_42;
    dest->_13 = a->_11 * b->_13 + a->_12 * b->_23 + a->_13 * b->_33 + a->_14 * b->_43;
    dest->_14 = a->_11 * b->_14 + a->_12 * b->_24 + a->_13 * b->_34 + a->_14 * b->_44;
    dest->_21 = a->_21 * b->_11 + a->_22 * b->_21 + a->_23 * b->_31 + a->_24 * b->_41;
    dest->_22 = a->_21 * b->_12 + a->_22 * b->_22 + a->_23 * b->_32 + a->_24 * b->_42;
    dest->_23 = a->_21 * b->_13 + a->_22 * b->_23 + a->_23 * b->_33 + a->_24 * b->_43;
    dest->_24 = a->_21 * b->_14 + a->_22 * b->_24 + a->_23 * b->_34 + a->_24 * b->_44;
    dest->_31 = a->_31 * b->_11 + a->_32 * b->_21 + a->_33 * b->_31 + a->_34 * b->_41;
    dest->_32 = a->_31 * b->_12 + a->_32 * b->_22 + a->_33 * b->_32 + a->_34 * b->_42;
    dest->_33 = a->_31 * b->_13 + a->_32 * b->_23 + a->_33 * b->_33 + a->_34 * b->_43;
    dest->_34 = a->_31 * b->_14 + a->_32 * b->_24 + a->_33 * b->_34 + a->_34 * b->_44;
    dest->_41 = a->_41 * b->_11 + a->_42 * b->_21 + a->_43 * b->_31 + a->_44 * b->_41;
    dest->_42 = a->_41 * b->_12 + a->_42 * b->_22 + a->_43 * b->_32 + a->_44 * b->_42;
    dest->_43 = a->_41 * b->_13 + a->_42 * b->_23 + a->_43 * b->_33 + a->_44 * b->_43;
    dest->_44 = a->_41 * b->_14 + a->_42 * b->_24 + a->_43 * b->_34 + a->_44 * b->_44;
}

void
mat4_transpose(mat4 *dest, const mat4 *m4)
{
    /* *
     * 11 12 13 14    11 21 31 41
     * 21 22 23 24 -> 12 22 32 42
     * 31 32 33 34    13 23 33 43
     * 41 42 43 44    14 24 34 44
     * */
    dest->_11 = m4->_11;
    dest->_12 = m4->_21;
    dest->_13 = m4->_31;
    dest->_14 = m4->_41;
    dest->_21 = m4->_12;
    dest->_22 = m4->_22;
    dest->_23 = m4->_32;
    dest->_24 = m4->_42;
    dest->_31 = m4->_13;
    dest->_32 = m4->_23;
    dest->_33 = m4->_33;
    dest->_34 = m4->_34;
    dest->_41 = m4->_14;
    dest->_42 = m4->_24;
    dest->_43 = m4->_34;
    dest->_44 = m4->_44;
}

void
mat4_cofactor(mat4 *dest, const mat4 *m4)
{
    dest->_11 =  m4->_11;
    dest->_12 = -m4->_12;
    dest->_13 =  m4->_13;
    dest->_14 = -m4->_14;
    dest->_21 = -m4->_21;
    dest->_22 =  m4->_22;
    dest->_23 = -m4->_23;
    dest->_24 =  m4->_24;
    dest->_31 =  m4->_31;
    dest->_32 = -m4->_32;
    dest->_33 =  m4->_33;
    dest->_34 = -m4->_34;
    dest->_41 = -m4->_41;
    dest->_42 =  m4->_42;
    dest->_43 = -m4->_43;
    dest->_44 =  m4->_44;
}

float
mat4_determinant(const mat4 *m4)
{
    float a =  m4->_11 * (m4->_22 * (m4->_33 * m4->_44 - m4->_34 * m4->_43) + -m4->_23 * (m4->_32 * m4->_44 - m4->_34 * m4->_42) + m4->_24 * (m4->_32 * m4->_43 - m4->_33 * m4->_42));
    float b = -m4->_12 * (m4->_21 * (m4->_33 * m4->_44 - m4->_34 * m4->_43) + -m4->_23 * (m4->_31 * m4->_44 - m4->_34 * m4->_41) + m4->_24 * (m4->_31 * m4->_43 - m4->_33 * m4->_41));
    float c =  m4->_13 * (m4->_21 * (m4->_32 * m4->_44 - m4->_34 * m4->_42) + -m4->_22 * (m4->_31 * m4->_44 - m4->_34 * m4->_41) + m4->_24 * (m4->_31 * m4->_42 - m4->_32 * m4->_41));
    float d = -m4->_14 * (m4->_21 * (m4->_32 * m4->_43 - m4->_33 * m4->_42) + -m4->_22 * (m4->_31 * m4->_43 - m4->_33 * m4->_41) + m4->_23 * (m4->_31 * m4->_42 - m4->_32 * m4->_41));
    return a + b + c + d;
}

bool
mat4_inverse(mat4 *dest, const mat4 *m4)
{
    float m4_det = mat4_determinant(m4);
    if(m4_det == 0.0f)
    {
        return false;
    }
    float det = 1.0f / m4_det;

	dest->_11 = (m4->_22 * m4->_33 * m4->_44 + m4->_23 * m4->_34 * m4->_42 + m4->_24 * m4->_32 * m4->_43 - m4->_22 * m4->_34 * m4->_43 - m4->_23 * m4->_32 * m4->_44 - m4->_24 * m4->_33 * m4->_42) * det;
	dest->_12 = (m4->_12 * m4->_34 * m4->_43 + m4->_13 * m4->_32 * m4->_44 + m4->_14 * m4->_33 * m4->_42 - m4->_12 * m4->_33 * m4->_44 - m4->_13 * m4->_34 * m4->_42 - m4->_14 * m4->_32 * m4->_43) * det;
	dest->_13 = (m4->_12 * m4->_23 * m4->_44 + m4->_13 * m4->_24 * m4->_42 + m4->_14 * m4->_22 * m4->_43 - m4->_12 * m4->_24 * m4->_43 - m4->_13 * m4->_22 * m4->_44 - m4->_14 * m4->_23 * m4->_42) * det;
	dest->_14 = (m4->_12 * m4->_24 * m4->_33 + m4->_13 * m4->_22 * m4->_34 + m4->_14 * m4->_23 * m4->_32 - m4->_12 * m4->_23 * m4->_34 - m4->_13 * m4->_24 * m4->_32 - m4->_14 * m4->_22 * m4->_33) * det;
    dest->_21 = (m4->_21 * m4->_34 * m4->_43 + m4->_23 * m4->_31 * m4->_44 + m4->_24 * m4->_33 * m4->_41 - m4->_21 * m4->_33 * m4->_44 - m4->_23 * m4->_34 * m4->_41 - m4->_24 * m4->_31 * m4->_43) * det;
	dest->_22 = (m4->_11 * m4->_33 * m4->_44 + m4->_13 * m4->_34 * m4->_41 + m4->_14 * m4->_31 * m4->_43 - m4->_11 * m4->_34 * m4->_43 - m4->_13 * m4->_31 * m4->_44 - m4->_14 * m4->_33 * m4->_41) * det;
	dest->_23 = (m4->_11 * m4->_24 * m4->_43 + m4->_13 * m4->_21 * m4->_44 + m4->_14 * m4->_23 * m4->_41 - m4->_11 * m4->_23 * m4->_44 - m4->_13 * m4->_24 * m4->_41 - m4->_14 * m4->_21 * m4->_43) * det;
	dest->_24 = (m4->_11 * m4->_23 * m4->_34 + m4->_13 * m4->_24 * m4->_31 + m4->_14 * m4->_21 * m4->_33 - m4->_11 * m4->_24 * m4->_33 - m4->_13 * m4->_21 * m4->_34 - m4->_14 * m4->_23 * m4->_31) * det;
    dest->_31 = (m4->_21 * m4->_32 * m4->_44 + m4->_22 * m4->_34 * m4->_41 + m4->_24 * m4->_31 * m4->_42 - m4->_21 * m4->_34 * m4->_42 - m4->_22 * m4->_31 * m4->_44 - m4->_24 * m4->_32 * m4->_41) * det;
	dest->_32 = (m4->_11 * m4->_34 * m4->_42 + m4->_12 * m4->_31 * m4->_44 + m4->_14 * m4->_32 * m4->_41 - m4->_11 * m4->_32 * m4->_44 - m4->_12 * m4->_34 * m4->_41 - m4->_14 * m4->_31 * m4->_42) * det;
	dest->_33 = (m4->_11 * m4->_22 * m4->_44 + m4->_12 * m4->_24 * m4->_41 + m4->_14 * m4->_21 * m4->_42 - m4->_11 * m4->_24 * m4->_42 - m4->_12 * m4->_21 * m4->_44 - m4->_14 * m4->_22 * m4->_41) * det;
	dest->_34 = (m4->_11 * m4->_24 * m4->_32 + m4->_12 * m4->_21 * m4->_34 + m4->_14 * m4->_22 * m4->_31 - m4->_11 * m4->_22 * m4->_34 - m4->_12 * m4->_24 * m4->_31 - m4->_14 * m4->_21 * m4->_32) * det;
    dest->_41 = (m4->_21 * m4->_33 * m4->_42 + m4->_22 * m4->_31 * m4->_43 + m4->_23 * m4->_32 * m4->_41 - m4->_21 * m4->_32 * m4->_43 - m4->_22 * m4->_33 * m4->_41 - m4->_23 * m4->_31 * m4->_42) * det;
    dest->_42 = (m4->_11 * m4->_32 * m4->_43 + m4->_12 * m4->_33 * m4->_41 + m4->_13 * m4->_31 * m4->_42 - m4->_11 * m4->_33 * m4->_42 - m4->_12 * m4->_31 * m4->_43 - m4->_13 * m4->_32 * m4->_41) * det;
	dest->_43 = (m4->_11 * m4->_23 * m4->_42 + m4->_12 * m4->_21 * m4->_43 + m4->_13 * m4->_22 * m4->_41 - m4->_11 * m4->_22 * m4->_43 - m4->_12 * m4->_23 * m4->_41 - m4->_13 * m4->_21 * m4->_42) * det;
	dest->_44 = (m4->_11 * m4->_22 * m4->_33 + m4->_12 * m4->_23 * m4->_31 + m4->_13 * m4->_21 * m4->_32 - m4->_11 * m4->_23 * m4->_32 - m4->_12 * m4->_21 * m4->_33 - m4->_13 * m4->_22 * m4->_31) * det;

    return true;
}

void
mat4_translation(mat4 *dest, float x, float y, float z)
{
    dest->_11 = 1.0f;
    dest->_12 = 0.0f;
    dest->_13 = 0.0f;
    dest->_14 = 0.0f;
    dest->_21 = 0.0f;
    dest->_22 = 1.0f;
    dest->_23 = 0.0f;
    dest->_24 = 0.0f;
    dest->_31 = 0.0f;
    dest->_32 = 0.0f;
    dest->_33 = 1.0f;
    dest->_34 = 0.0f;
    dest->_41 = x;
    dest->_42 = y;
    dest->_43 = z;
    dest->_44 = 1.0f;
}

void
mat4_rotate_x(mat4 *dest, float by)
{
    by = DEG2RAD(by);

    dest->_11 = 1.0f;
    dest->_12 = 0.0f;
    dest->_13 = 0.0f;
    dest->_14 = 0.0f;
    dest->_21 = 0.0f;
    dest->_22 = cosf(by);
    dest->_23 = sinf(by);
    dest->_24 = 0.0f;
    dest->_31 = 0.0f;
    dest->_32 = -sinf(by);
    dest->_33 = cosf(by);
    dest->_34 = 0.0f;
    dest->_41 = 0.0f;
    dest->_42 = 0.0f;
    dest->_43 = 0.0f;
    dest->_44 = 1.0f;
}

void
mat4_rotate_y(mat4 *dest, float by)
{
    by = DEG2RAD(by);

    dest->_11 = cosf(by);
    dest->_12 = 0.0f;
    dest->_13 = -sinf(by);
    dest->_14 = 0.0f;
    dest->_21 = 0.0f;
    dest->_22 = 1.0f;
    dest->_23 = 0.0f;
    dest->_24 = 0.0f;
    dest->_31 = sinf(by);
    dest->_32 = 0.0f;
    dest->_33 = cosf(by);
    dest->_34 = 0.0f;
    dest->_41 = 0.0f;
    dest->_42 = 0.0f;
    dest->_43 = 0.0f;
    dest->_44 = 1.0f;
}

void
mat4_rotate_z(mat4 *dest, float by)
{
    by = DEG2RAD(by);

    dest->_11 = cosf(by);
    dest->_12 = sinf(by);
    dest->_13 = 0.0f;
    dest->_14 = 0.0f;
    dest->_21 = -sinf(by);
    dest->_22 = cosf(by);
    dest->_23 = 0.0f;
    dest->_24 = 0.0f;
    dest->_31 = 0.0f;
    dest->_32 = 0.0f;
    dest->_33 = 1.0f;
    dest->_34 = 0.0f;
    dest->_41 = 0.0f;
    dest->_42 = 0.0f;
    dest->_43 = 0.0f;
    dest->_44 = 1.0f;
}

void
mat4_print(const mat4 *m4)
{
    printf("%.3f, %.3f, %.3f, %.3f \n", m4->_11, m4->_12, m4->_13, m4->_14);
    printf("%.3f, %.3f, %.3f, %.3f \n", m4->_21, m4->_22, m4->_23, m4->_24);
    printf("%.3f, %.3f, %.3f, %.3f \n", m4->_31, m4->_32, m4->_33, m4->_34);
    printf("%.3f, %.3f, %.3f, %.3f \n", m4->_41, m4->_42, m4->_43, m4->_44);
}

