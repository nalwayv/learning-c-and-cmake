#ifndef _MATRIX_H_
#define _MATRIX_H_

struct mat2 {
    union {
        struct {
            float _11, _12,
                  _21, _22;
        };

        float arr[4];
    };
};
typedef struct mat2 mat2;

struct mat3 {
    union {
        struct{
            float _11, _12, _13,
                  _21, _22, _23,
                  _31, _32, _33;
        };

        float arr[9];
    };
};
typedef struct mat3 mat3;

struct mat4 {
    union {
        struct {
            float _11, _12, _13, _14,
                  _21, _22, _23, _24,
                  _31, _32, _33, _34,
                  _41, _42, _43, _44;
        };

        float arr[16];
    };
};
typedef struct mat4 mat4;

/* MAT 2 */

mat2 mat2_zero();
mat2 mat2_identity();
mat2 mat2_new(float _11, float _12,
              float _21, float _22);
void mat2_add(mat2 *dest, const mat2 *a, const mat2 *b);
void mat2_sub(mat2 *dest, const mat2 *a, const mat2 *b);
void mat2_scale(mat2 *dest, const mat2 *m2, float by);
void mat2_mul(mat2 *dest, const mat2 *a, const mat2 *b);
void mat2_transpose(mat2 *dest, const mat2 *m2);
void mat2_cofactor(mat2 *dest, const mat2 *m2);
float mat2_determinant(const mat2 *m2);
int mat2_inverse(mat2 *dest, const mat2 *m2);
void mat2_print(const mat2 *m2);

/* MAT 3 */

mat3 mat3_zero();
mat3 mat3_identity();
mat3 mat3_new(float _11, float _12, float _13,
              float _21, float _22, float _23,
              float _31, float _32, float _33);
void mat3_add(mat3 *dest, const mat3 *a, const mat3 *b);
void mat3_sub(mat3 *dest, const mat3 *a, const mat3 *b);
void mat3_scale(mat3 *dest, const mat3 *m3, float by);
void mat3_mul(mat3 *dest, const mat3 *m3, const mat3 *b);
void mat3_transpose(mat3 *dest, const mat3 *m3);
void mat3_cofactor(mat3 *dest, const mat3 *m3);
float mat3_determinant(const mat3 *m3);
int  mat3_inverse(mat3 *dest, const mat3 *m3);
void mat3_print(const mat3 *m3);

/* MAT 4*/

mat4 mat4_zero();
mat4 mat4_identity();
mat4 mat4_new(float _11, float _12, float _13, float _14,
              float _21, float _22, float _23, float _24,
              float _31, float _32, float _33, float _34,
              float _41, float _42, float _43, float _44);
void mat4_add(mat4 *dest, const mat4 *a, const mat4 *b);
void mat4_sub(mat4 *dest, const mat4 *a, const mat4 *b);
void mat4_scale(mat4 *dest, const mat4 *a, float by);
void mat4_mul(mat4 *dest, const mat4 *a, const mat4 *b);
void mat4_transpose(mat4 *dest, const mat4 *m4);
void mat4_cofactor(mat4 *dest, const mat4 *m4);
float mat4_determinant(const mat4 *m4);
int mat4_inverse(mat4 *dest, const mat4 *m4);
void mat4_print(const mat4 *m4);

#endif /**/
