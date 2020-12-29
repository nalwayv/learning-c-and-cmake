#ifndef _MATRIX_H_
#define _MATRIX_H_

#include <stdbool.h>

/**
 * matrix 2x2
 * */
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

/* *
 * matrix 3x3
 * */
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

/* *
 * matrix 4x4
 * */
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

/* zero matrix 2x2 */
mat2 mat2_zero();
/* identity matrix 2x2 */
mat2 mat2_identity();
/* create new matrix 2x2 */
mat2 mat2_new(float _11, float _12, float _21, float _22);
/* copy matrix 2x2 into dest */
void mat2_cpy(mat2 *dest, const mat2 *m2);
/* add matrix 2x2 */
void mat2_add(mat2 *dest, const mat2 *a, const mat2 *b);
/* sub matrix 2x2 */
void mat2_sub(mat2 *dest, const mat2 *a, const mat2 *b);
/* scale matrix 2x2 by value */
void mat2_scale(mat2 *dest, const mat2 *m2, float by);
/* multiply matrix 2x2 by another matrix 2x2 */
void mat2_mul(mat2 *dest, const mat2 *a, const mat2 *b);
/* transpose matrix 2x2 */
void mat2_transpose(mat2 *dest, const mat2 *m2);
/* get cofactor of matrix 2x2 */
void mat2_cofactor(mat2 *dest, const mat2 *m2);
/* cut a matrix 2x2 into single float value */
float mat2_cut(const mat2 *m2, int row, int col);
/* get minor for matrix 2x2 */
void mat2_minor(mat2 *dest, const mat2 *m2);
/* get adjugate of matrix 2x2 */
void mat2_adjugate(mat2 *dest, const mat2 *m2);
/* get determinant of matrix 2x2 */
float mat2_determinant(const mat2 *m2);
/* inverse matrix 2x2 */
bool mat2_inverse(mat2 *dest, const mat2 *m2);
/* rotate matrix 2x2 on z axis */
void mat2_rotate_z(mat2 *dest, float by);
/* print matrix 2x2 to console */
void mat2_print(const mat2 *m2);

/* MAT 3 */

/* zero matrix 3x3 */
mat3 mat3_zero();
/* identity matrix 3x3 */
mat3 mat3_identity();
/* create new matrix 3x3 */
mat3 mat3_new(float _11, float _12, float _13, float _21, float _22, float _23,
              float _31, float _32, float _33);
/* copy matrix 3x3 into dest */
void mat3_cpy(mat3 *dest, const mat3 *m3);
/* add matrix 3x3 */
void mat3_add(mat3 *dest, const mat3 *a, const mat3 *b);
/* sub matrix 3x3 */
void mat3_sub(mat3 *dest, const mat3 *a, const mat3 *b);
/* scale matrix 3x3 by value */
void mat3_scale(mat3 *dest, const mat3 *m3, float by);
/* multiply matrix 3x3 by another matrix 3x3 */
void mat3_mul(mat3 *dest, const mat3 *m3, const mat3 *b);
/* transpose matrix 3x3 */
void mat3_transpose(mat3 *dest, const mat3 *m3);
/* get cofactor of matrix 3x3 */
void mat3_cofactor(mat3 *dest, const mat3 *m3);
/* cut a matrix 3x3 into a matrix 2x2 */
void mat3_cut(mat2 *dest, const mat3 *m3, int row, int col);
/* get minor matrix 3x3 */
void mat3_minor(mat3 *dest, const mat3 *m3);
/* get adjugate of matrix 3x3 */
void mat3_adjugate(mat3 *dest, const mat3 *m3);
/* get detrminant of matrix 3x3 */
float mat3_determinant(const mat3 *m3);
/* inverse matrix 3x3 */
bool mat3_inverse(mat3 *dest, const mat3 *m3);
/* scale matrix 3x3 by vector 3 */
void mat3_scaling(mat3 *dest, float x, float y, float z);
/* rotate matrix 3x3 on x axis */
void mat3_rotate_x(mat3 *dest, float by);
/* rotate matrix 3x3 on y axis */
void mat3_rotate_y(mat3 *dest, float by);
/* rotate matrix 3x3 on z axis */
void mat3_rotate_z(mat3 *dest, float by);
/* rotate matrix 3x3 by pitch yaw and roll */
void mat3_rotation(mat3 *dest, float pitch, float yaw, float roll);
/* print matrix 3x3 to console */
void mat3_print(const mat3 *m3);

/* MAT 4 */

/* zero matrix 4x4 */
mat4 mat4_zero();
/* identity matrix 5x4 */
mat4 mat4_identity();
/* create new matrix 4x4 */
mat4 mat4_new(float _11, float _12, float _13, float _14, float _21, float _22,
              float _23, float _24, float _31, float _32, float _33, float _34,
              float _41, float _42, float _43, float _44);
/* copy matrix 4x4 into dest */
void mat4_cpy(mat4 *dest, const mat4 *v4);
/* add matrix 4x4 */
void mat4_add(mat4 *dest, const mat4 *a, const mat4 *b);
/* sub matrix 4x4 */
void mat4_sub(mat4 *dest, const mat4 *a, const mat4 *b);
/* scale matrix 4x4 by value */
void mat4_scale(mat4 *dest, const mat4 *a, float by);
/* multiply matrix 4x4 by another matrix 4x4 */
void mat4_mul(mat4 *dest, const mat4 *a, const mat4 *b);
/* transpose matrix 4x4 */
void mat4_transpose(mat4 *dest, const mat4 *m4);
/* get cofactor of matrix 4x4 */
void mat4_cofactor(mat4 *dest, const mat4 *m4);
/* cut a matrix 4x4 into a matrix 3x3 */
void mat4_cut(mat3 *dest, const mat4 *m4, int row, int col);
/* get minor matrix 4x4 */
void mat4_minor(mat4 *dest, const mat4 *m4);
/* get adjugate of matrix 4x4 */
void mat4_adjugate(mat4 *dest, const mat4 *m4);
/* get determinant of matrix 4x4 */
float mat4_determinant(const mat4 *m4);
/* inverse matrix 4x4 */
bool mat4_inverse(mat4 *dest, const mat4 *m4);
/* translate matrix 4x4 by vector 4 */
void mat4_translation(mat4 *dest, float x, float y, float z);
/* scale matrix 4x4 by vector 4 */
void mat4_scaling(mat4 *dest, float x, float y, float z);
/* rotate matrix 4x4 on x axis */
void mat4_rotate_x(mat4 *dest, float by);
/* rotate matrix 4x4 on y axis */
void mat4_rotate_y(mat4 *dest, float by);
/* rotate matrix 4x4 on z axis */
void mat4_rotate_z(mat4 *dest, float by);
/* rotate matrix 4x4 by pitch yaw and roll */
void mat4_rotation(mat4 *dest, float pitch, float yaw, float roll);
/* print matrix 4x4 to console */
void mat4_print(const mat4 *m4);

#endif /**/
