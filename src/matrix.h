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

/* MAT 2 */

struct mat2 mat2_zero();
struct mat2 mat2_identity();
struct mat2 mat2_new(float _11, float _12,
                     float _21, float _22);
struct mat2 mat2_add(const struct mat2 *a, const struct mat2 *b);
struct mat2 mat2_sub(const struct mat2 *a, const struct mat2 *b);
struct mat2 mat2_scale(const struct mat2 *a, float by);
struct mat2 mat2_mul(const struct mat2 *a, const struct mat2 *b);
struct mat2 mat2_transpose(const struct mat2 *m2);
struct mat2 mat2_cofactor(const struct mat2 *m2);
float mat2_determinant(const struct mat2 *m2);
struct mat2 mat2_inverse(const struct mat2 *m2);
void mat2_print(const struct mat2 *m2);

/* MAT 3 */

struct mat3 mat3_zero();
struct mat3 mat3_identity();
struct mat3 mat3_new(float _11, float _12, float _13,
                     float _21, float _22, float _23,
                     float _31, float _32, float _33);
struct mat3 mat3_add(const struct mat3 *a, const struct mat3 *b);
struct mat3 mat3_sub(const struct mat3 *a, const struct mat3 *b);
struct mat3 mat3_scale(const struct mat3 *a, float by);
struct mat3 mat3_mul(const struct mat3 *a, const struct mat3 *b);
struct mat3 mat3_transpose(const struct mat3 *m3);
struct mat3 mat3_cofactor(const struct mat3 *m3);
float mat3_determinant(const struct mat3 *m3);
struct mat3 mat3_inverse(const struct mat3 *m3);
void mat3_print(const struct mat3 *m3);

/* MAT 4*/

struct mat4 mat4_zero();
struct mat4 mat4_identity();
struct mat4 mat4_new(float _11, float _12, float _13, float _14,
                     float _21, float _22, float _23, float _24,
                     float _31, float _32, float _33, float _34,
                     float _41, float _42, float _43, float _44);
struct mat4 mat4_add(const struct mat4 *a, const struct mat4 *b);
struct mat4 mat4_sub(const struct mat4 *a, const struct mat4 *b);
struct mat4 mat4_scale(const struct mat4 *a, float by);
struct mat4 mat4_mul(const struct mat4 *a, const struct mat4 *b);
struct mat4 mat4_transpose(const struct mat4 *m4);
struct mat4 mat4_cofactor(const struct mat4 *m4);
float mat4_determinant(const struct mat4 *m4);
struct mat4 mat4_inverse(const struct mat4 *m4);
void mat4_print(const struct mat4 *m4);

#endif /**/
