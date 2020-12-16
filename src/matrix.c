/* *
 * TODO(...)
 * */
#include <stdbool.h>
#include <stdio.h>
#include <math.h>
#include <assert.h>

#include "matrix.h"

static bool
matrix_multiply(float *out, const float *arrA, int colA, int rowA, const float *arrB, int rowB, int colB)
{
    if(colA != rowB)
    {
        return false;
    }

    for(int i = 0; i < rowA; ++i)
    {
        for (int j = 0; j < colB; ++j)
        {
            int id = (colB * i) + j;
            out[id] = 0.0f;

            for(int k = 0; k < rowB; ++k)
            {
                int idA = (colA * i) + k;
                int idB = (colB * k) + j;
                out[id] += arrA[idA] * arrB[idB];
            }
        }
    }

    return true;
}

static void
matrix_transpose(float *out, const float* arr, int cols, int rows)
{
    int n = cols*rows;
    for(int i = 0; i < n; ++i)
    {
        int r = i / rows;
        int c = i % rows;
        int id = (cols * c) + r;
        out[i] = arr[id];
    }
}

static void
matrix_print(const float *arr, int row, int col)
{
    for(int i = 0; i < row; ++i)
    {
        for(int j = 0; j < col; ++j)
        {
            int id = (col * i) + j;
            printf("%.4f ", arr[id]);
        }
        printf("\n");
    }
}


/* MAT 2 */

struct mat2
mat2_zero(){
    struct mat2 m2 = (struct mat2){
        ._11 = 0.0f, ._12 = 0.0f,
        ._21 = 0.0f, ._22 = 0.0f
    };
    return m2;
}

struct mat2
mat2_identity(){
    struct mat2 m2 = (struct mat2){
        ._11 = 1.0f, ._12 = 0.0f,
        ._21 = 0.0f, ._22 = 1.0f
    };
    return m2;
}

struct mat2
mat2_new(float _11, float _12,
         float _21, float _22)
{
    struct mat2 m2 = (struct mat2){
        ._11 = _11, ._12 = _12,
        ._21 = _21, ._22 = _22
    };
    return m2;
}

struct mat2
mat2_add(const struct mat2 *a, const struct mat2 *b)
{
    struct mat2 m2;
    for(int i = 0; i < 4; ++i)
    {
        m2.arr[i] = a->arr[i] + b->arr[i];
    }
    return m2;
}

struct mat2
mat2_sub(const struct mat2 *a, const struct mat2 *b)
{
    struct mat2 m2;
    for(int i = 0; i < 4; ++i)
    {
        m2.arr[i] = a->arr[i] - b->arr[i];
    }
    return m2;
}

struct mat2
mat2_scale(const struct mat2 *a, float by)
{
    struct mat2 m2;
    for(int i = 0; i < 4; ++i)
    {
        m2.arr[i] = a->arr[i] * by;
    }
    return m2;
}

struct mat2
mat2_mul(const struct mat2 *a, const struct mat2 *b)
{
    struct mat2 dest;
    bool check = matrix_multiply(dest.arr, a->arr, 2, 2, b->arr, 2, 2);
    assert(check != false);
    return dest;
}

struct mat2
mat2_transpose(const struct mat2 *m2)
{
    /* *
     * a b -> a c
     * c d    b d
     * */
    struct mat2 dest;
    matrix_transpose(dest.arr, m2->arr, 2, 2);
    return dest;
}

struct mat2
mat2_cofactor(const struct mat2 *m2)
{
    return mat2_new(m2->_11, -m2->_12, 
                    -m2->_21, m2->_22);
}

float
mat2_determinant(const struct mat2 *m2)
{
    return (m2->_11 * m2->_22) - (m2->_12 * m2->_21);
}

struct mat2
mat2_inverse(const struct mat2 *m2)
{
    float det = 1.0f / mat2_determinant(m2);
    assert(det != 0.0f);

    float a =  (m2->_22) * det; 
    float b = -(m2->_12) * det;
    float c = -(m2->_21) * det;
    float d =  (m2->_11) * det;
    
    return mat2_new(a, b, 
                    c, d);
}

void
mat2_print(const struct mat2 *m2)
{
    matrix_print(m2->arr, 2, 2);
}

/* MAT 3 */

struct mat3
mat3_zero(){
    struct mat3 m3 = (struct mat3){
        ._11 = 0.0f, ._12 = 0.0f, ._13 = 0.0f,
        ._21 = 0.0f, ._22 = 0.0f, ._23 = 0.0f,
        ._31 = 0.0f, ._32 = 0.0f, ._33 = 0.0f
    };
    return m3;
}

struct mat3
mat3_identity(){
    struct mat3 m3 = (struct mat3){
        ._11 = 1.0f, ._12 = 0.0f, ._13 = 0.0f,
        ._21 = 0.0f, ._22 = 1.0f, ._23 = 0.0f,
        ._31 = 0.0f, ._32 = 0.0f, ._33 = 1.0f
    };
    return m3;
}

struct mat3
mat3_new(float _11, float _12, float _13,
         float _21, float _22, float _23,
         float _31, float _32, float _33)
{
    struct mat3 m3 = (struct mat3){
        ._11 = _11, ._12 = _12, ._13 = _13,
        ._21 = _21, ._22 = _22, ._23 = _23,
        ._31 = _31, ._32 = _32, ._33 = _33
    };
    return m3;
}

struct mat3
mat3_add(const struct mat3 *a, const struct mat3 *b)
{
    struct mat3 m3;
    for(int i = 0; i < 9; ++i)
    {
        m3.arr[i] = a->arr[i] + b->arr[i];
    }
    return m3;
}

struct mat3
mat3_sub(const struct mat3 *a, const struct mat3 *b)
{
    struct mat3 m3;
    for(int i = 0; i < 9; ++i)
    {
        m3.arr[i] = a->arr[i] - b->arr[i];
    }
    return m3;
}

struct mat3
mat3_scale(const struct mat3 *a, float by)
{
    struct mat3 m3;
    for(int i = 0; i < 9; ++i)
    {
        m3.arr[i] = a->arr[i] * by;
    }
    return m3;
}

struct mat3
mat3_mul(const struct mat3 *a, const struct mat3 *b)
{
    struct mat3 dest;
    bool check = matrix_multiply(dest.arr, a->arr, 3, 3, b->arr, 3, 3);
    assert(check != false);
    return dest;
}

struct mat3
mat3_transpose(const struct mat3 *m3)
{
    /* *
     * a b c    a d g
     * d e f -> b e h
     * g h i    c f i
     * */
    struct mat3 dest;
    matrix_transpose(dest.arr, m3->arr, 3, 3);
    return dest;
}

struct mat3
mat3_cofactor(const struct mat3 *m3)
{
    return mat3_new( m3->_11, -m3->_12,  m3->_13,
                    -m3->_21,  m3->_22, -m3->_23,
                     m3->_31, -m3->_32,  m3->_33);
}

float mat3_determinant(const struct mat3 *m3)
{
    float a = m3->_11 * (m3->_22 * m3->_33 - m3->_23 * m3->_32);
    float b = -m3->_12 * (m3->_21 * m3->_33 - m3->_23 * m3->_31);
    float c = m3->_13 * (m3->_21 * m3->_32 - m3->_22 * m3->_31);
    return a + b + c;
}

struct mat3
mat3_inverse(const struct mat3 *m3)
{   
    float det = 1.0f / mat3_determinant(m3);
    assert(det != 0.0f);

    float a =  (m3->_22 * m3->_33 - m3->_23 * m3->_32) * det;
    float b = -(m3->_12 * m3->_33 - m3->_32 * m3->_13) * det;
    float c =  (m3->_12 * m3->_23 - m3->_22 * m3->_13) * det;

    float d = -(m3->_21 * m3->_33 - m3->_31 * m3->_23) * det;
    float e =  (m3->_11 * m3->_33 - m3->_13 * m3->_31) * det;
    float f = -(m3->_11 * m3->_23 - m3->_21 * m3->_13) * det;

    float g =  (m3->_21 * m3->_32 - m3->_31 * m3->_22) * det;
    float h = -(m3->_11 * m3->_32 - m3->_31 * m3->_12) * det;
    float i =  (m3->_11 * m3->_22 - m3->_12 * m3->_21) * det;

    return mat3_new(a, b, c, 
                    d, e, f, 
                    g, h, i);
}

void
mat3_print(const struct mat3 *m3)
{
    matrix_print(m3->arr, 3, 3);
}

/* MAT 4 */

struct mat4
mat4_zero(){
    struct mat4 m4 = (struct mat4){
        ._11 = 0.0f, ._12 = 0.0f, ._13 = 0.0f, ._14 = 0.0f,
        ._21 = 0.0f, ._22 = 0.0f, ._23 = 0.0f, ._24 = 0.0f,
        ._31 = 0.0f, ._32 = 0.0f, ._33 = 0.0f, ._34 = 0.0f,
        ._41 = 0.0f, ._42 = 0.0f, ._43 = 0.0f, ._44 = 0.0f
    };
    return m4;
}

struct mat4
mat4_identity(){
    struct mat4 m4 = (struct mat4){
        ._11 = 1.0f, ._12 = 0.0f, ._13 = 0.0f, ._14 = 0.0f,
        ._21 = 0.0f, ._22 = 1.0f, ._23 = 0.0f, ._24 = 0.0f,
        ._31 = 0.0f, ._32 = 0.0f, ._33 = 1.0f, ._34 = 0.0f,
        ._41 = 0.0f, ._42 = 0.0f, ._43 = 0.0f, ._44 = 1.0f
    };
    return m4;
}

struct mat4
mat4_new(float _11, float _12, float _13, float _14,
         float _21, float _22, float _23, float _24,
         float _31, float _32, float _33, float _34,
         float _41, float _42, float _43, float _44)
{
    struct mat4 m4 = (struct mat4){
        ._11 = _11, ._12 = _12, ._13 = _13, ._14 = _14,
        ._21 = _21, ._22 = _22, ._23 = _23, ._24 = _24,
        ._31 = _31, ._32 = _32, ._33 = _33, ._34 = _34,
        ._41 = _41, ._42 = _42, ._43 = _43, ._44 = _44
    };
    return m4;
}

struct mat4
mat4_add(const struct mat4 *a, const struct mat4 *b)
{
    struct mat4 m4 = mat4_zero();
    for(int i = 0; i < 16; ++i)
    {
        m4.arr[i] = a->arr[i] + b->arr[i];
    }
    return m4;
}

struct mat4
mat4_sub(const struct mat4 *a, const struct mat4 *b)
{
    struct mat4 m4 = mat4_zero();
    for(int i = 0; i < 16; ++i)
    {
        m4.arr[i] = a->arr[i] - b->arr[i];
    }
    return m4;
}

struct mat4
mat4_scale(const struct mat4 *a, float by)
{
    struct mat4 m4 = mat4_zero();
    for(int i = 0; i < 16; ++i)
    {
        m4.arr[i] = a->arr[i] * by;
    }
    return m4;
}

struct mat4
mat4_mul(const struct mat4 *a, const struct mat4 *b)
{
    struct mat4 dest = mat4_zero();
    bool check = matrix_multiply(dest.arr, a->arr, 4, 4, b->arr, 4, 4);
    assert(check != false);
    return dest;
}

struct mat4
mat4_transpose(const struct mat4 *m4)
{
    /* *
     * a b c d    a e i m
     * e f g h -> b f j n
     * i j k l    c g k o
     * m n o p    d h l p
     * */
    struct mat4 dest = mat4_zero();
    matrix_transpose(dest.arr, m4->arr, 4, 4);
    return dest;
}

struct mat4
mat4_cofactor(const struct mat4 *m4)
{
    return mat4_new( m4->_11, -m4->_12,  m4->_13, -m4->_14,
                    -m4->_21,  m4->_22, -m4->_23,  m4->_24,
                     m4->_31, -m4->_32,  m4->_33, -m4->_34,
                    -m4->_41,  m4->_42, -m4->_43,  m4->_44);
}

float
mat4_determinant(const struct mat4 *m4)
{
    float a = m4->_11 * (m4->_22 * (m4->_33 * m4->_44 - m4->_34 * m4->_43) + -m4->_23 * (m4->_32 * m4->_44 - m4->_34 * m4->_42) + m4->_24 * (m4->_32 * m4->_43 - m4->_33 * m4->_42));
    float b = -m4->_12 * (m4->_21 * (m4->_33 * m4->_44 - m4->_34 * m4->_43) + -m4->_23 * (m4->_31 * m4->_44 - m4->_34 * m4->_41) + m4->_24 * (m4->_31 * m4->_43 - m4->_33 * m4->_41));
    float c = m4->_13 * (m4->_21 * (m4->_32 * m4->_44 - m4->_34 * m4->_42) + -m4->_22 * (m4->_31 * m4->_44 - m4->_34 * m4->_41) + m4->_24 * (m4->_31 * m4->_42 - m4->_32 * m4->_41));
    float d = -m4->_14 * (m4->_21 * (m4->_32 * m4->_43 - m4->_33 * m4->_42) + -m4->_22 * (m4->_31  * m4->_43 - m4->_33 * m4->_41) + m4->_23 * (m4->_31 * m4->_42 - m4->_32 * m4->_41));
    return a + b + c + d;
}

struct mat4
mat4_inverse(const struct mat4 *m4)
{
    float det = 1.0f / mat4_determinant(m4);
    assert(det != 0.0f);

	float a = (m4->_22 * m4->_33 * m4->_44 + m4->_23 * m4->_34 * m4->_42 + m4->_24 * m4->_32 * m4->_43 - m4->_22 * m4->_34 * m4->_43 - m4->_23 * m4->_32 * m4->_44 - m4->_24 * m4->_33 * m4->_42) * det;
	float b = (m4->_12 * m4->_34 * m4->_43 + m4->_13 * m4->_32 * m4->_44 + m4->_14 * m4->_33 * m4->_42 - m4->_12 * m4->_33 * m4->_44 - m4->_13 * m4->_34 * m4->_42 - m4->_14 * m4->_32 * m4->_43) * det;
	float c = (m4->_12 * m4->_23 * m4->_44 + m4->_13 * m4->_24 * m4->_42 + m4->_14 * m4->_22 * m4->_43 - m4->_12 * m4->_24 * m4->_43 - m4->_13 * m4->_22 * m4->_44 - m4->_14 * m4->_23 * m4->_42) * det;
	float d = (m4->_12 * m4->_24 * m4->_33 + m4->_13 * m4->_22 * m4->_34 + m4->_14 * m4->_23 * m4->_32 - m4->_12 * m4->_23 * m4->_34 - m4->_13 * m4->_24 * m4->_32 - m4->_14 * m4->_22 * m4->_33) * det;
	float e = (m4->_21 * m4->_34 * m4->_43 + m4->_23 * m4->_31 * m4->_44 + m4->_24 * m4->_33 * m4->_41 - m4->_21 * m4->_33 * m4->_44 - m4->_23 * m4->_34 * m4->_41 - m4->_24 * m4->_31 * m4->_43) * det;
	float f = (m4->_11 * m4->_33 * m4->_44 + m4->_13 * m4->_34 * m4->_41 + m4->_14 * m4->_31 * m4->_43 - m4->_11 * m4->_34 * m4->_43 - m4->_13 * m4->_31 * m4->_44 - m4->_14 * m4->_33 * m4->_41) * det;
	float g = (m4->_11 * m4->_24 * m4->_43 + m4->_13 * m4->_21 * m4->_44 + m4->_14 * m4->_23 * m4->_41 - m4->_11 * m4->_23 * m4->_44 - m4->_13 * m4->_24 * m4->_41 - m4->_14 * m4->_21 * m4->_43) * det;
	float h = (m4->_11 * m4->_23 * m4->_34 + m4->_13 * m4->_24 * m4->_31 + m4->_14 * m4->_21 * m4->_33 - m4->_11 * m4->_24 * m4->_33 - m4->_13 * m4->_21 * m4->_34 - m4->_14 * m4->_23 * m4->_31) * det;
	float i = (m4->_21 * m4->_32 * m4->_44 + m4->_22 * m4->_34 * m4->_41 + m4->_24 * m4->_31 * m4->_42 - m4->_21 * m4->_34 * m4->_42 - m4->_22 * m4->_31 * m4->_44 - m4->_24 * m4->_32 * m4->_41) * det;
	float j = (m4->_11 * m4->_34 * m4->_42 + m4->_12 * m4->_31 * m4->_44 + m4->_14 * m4->_32 * m4->_41 - m4->_11 * m4->_32 * m4->_44 - m4->_12 * m4->_34 * m4->_41 - m4->_14 * m4->_31 * m4->_42) * det;
	float k = (m4->_11 * m4->_22 * m4->_44 + m4->_12 * m4->_24 * m4->_41 + m4->_14 * m4->_21 * m4->_42 - m4->_11 * m4->_24 * m4->_42 - m4->_12 * m4->_21 * m4->_44 - m4->_14 * m4->_22 * m4->_41) * det;
	float l = (m4->_11 * m4->_24 * m4->_32 + m4->_12 * m4->_21 * m4->_34 + m4->_14 * m4->_22 * m4->_31 - m4->_11 * m4->_22 * m4->_34 - m4->_12 * m4->_24 * m4->_31 - m4->_14 * m4->_21 * m4->_32) * det;
	float m = (m4->_21 * m4->_33 * m4->_42 + m4->_22 * m4->_31 * m4->_43 + m4->_23 * m4->_32 * m4->_41 - m4->_21 * m4->_32 * m4->_43 - m4->_22 * m4->_33 * m4->_41 - m4->_23 * m4->_31 * m4->_42) * det;
    float n = (m4->_11 * m4->_32 * m4->_43 + m4->_12 * m4->_33 * m4->_41 + m4->_13 * m4->_31 * m4->_42 - m4->_11 * m4->_33 * m4->_42 - m4->_12 * m4->_31 * m4->_43 - m4->_13 * m4->_32 * m4->_41) * det;
	float o = (m4->_11 * m4->_23 * m4->_42 + m4->_12 * m4->_21 * m4->_43 + m4->_13 * m4->_22 * m4->_41 - m4->_11 * m4->_22 * m4->_43 - m4->_12 * m4->_23 * m4->_41 - m4->_13 * m4->_21 * m4->_42) * det;
	float p = (m4->_11 * m4->_22 * m4->_33 + m4->_12 * m4->_23 * m4->_31 + m4->_13 * m4->_21 * m4->_32 - m4->_11 * m4->_23 * m4->_32 - m4->_12 * m4->_21 * m4->_33 - m4->_13 * m4->_22 * m4->_31) * det;

    return mat4_new(a, b, c, d,
                    e, f, g, h,
                    i, j, k, l,
                    m, n, o, p);
}

void
mat4_print(const struct mat4 *m4)
{
    matrix_print(m4->arr, 4, 4);
}

