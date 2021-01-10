/* *
 * matrix
 * */
#include <math.h>
#include <stdio.h>

#include "matrix.h"
#include "utils.h"

/* HELPERS */

internal bool matrix_compar(const float *a, const float *b, int row, int col) {
  int n = row*col;
  for (int i = 0; i < n; ++i) {
    if (a[i] != b[i]) {
      return false;
    }
  }
  return true;
}

internal void matrix_cpy(float *dest, const float *src, int row, int col) {
  int n = row*col;
  for (int i = 0; i < n; ++i) {
    dest[i] = src[i];
  }
}

/* generic transpose matrix */
internal void matrix_transpose(float *dest, const float *src, int row,
                               int col) {
  int n = row * col;
  for (int i = 0; i < n; ++i) {
    int r = i / row;
    int c = i % row;
    int id = col * c + r;
    dest[i] = src[id];
  }
}

/* generic cofactor matrix */
internal void matrix_cofactor(float *dest, const float *src, int row, int col) {
  for (int i = 0; i < row; ++i) {
    for (int j = 0; j < col; ++j) {
      int id = col * i + j;
      // + or -
      dest[id] = src[id] * powf(-1.0, i + j);
    }
  }
}

/* generic multiply matrix */
internal void matrix_multiply(float *dest, const float *mat_a,
                              const float *mat_b, int row, int col) {
  for (int i = 0; i < row; ++i) {
    for (int j = 0; j < col; ++j) {

      int id = col * i + j;
      dest[id] = 0.0f;

      for (int k = 0; k < row; ++k) {
        int id_a = col * i + k;
        int id_b = col * k + j;
        dest[id] += mat_a[id_a] * mat_b[id_b];
      }
    }
  }
}

/* MAT 2 */

mat2 mat2_zero() {
  mat2 m2 = (mat2){._11 = 0.0f, ._12 = 0.0f, ._21 = 0.0f, ._22 = 0.0f};
  return m2;
}

mat2 mat2_identity() {
  mat2 m2 = (mat2){._11 = 1.0f, ._12 = 0.0f, ._21 = 0.0f, ._22 = 1.0f};
  return m2;
}

mat2 mat2_new(float _11, float _12, float _21, float _22) {
  mat2 m2 = (mat2){._11 = _11, ._12 = _12, ._21 = _21, ._22 = _22};
  return m2;
}

void mat2_cpy(mat2 *dest, const mat2 *m2) {
  dest->_11 = m2->_11;
  dest->_12 = m2->_12;

  dest->_21 = m2->_21;
  dest->_22 = m2->_22;
}

mat2 mat2_add(const mat2 *a, const mat2 *b) {
  return mat2_new(a->_11 + b->_11, a->_12 + b->_12, a->_21 + b->_21,
                  a->_22 + b->_22);
}

mat2 mat2_sub(const mat2 *a, const mat2 *b) {
  return mat2_new(a->_11 - b->_11, a->_12 - b->_12, a->_21 - b->_21,
                  a->_22 - b->_22);
}

mat2 mat2_scale(const mat2 *m2, float by) {
  return mat2_new(m2->_11 * by, m2->_12 * by, m2->_21 * by, m2->_22 * by);
}

mat2 mat2_mul(const mat2 *a, const mat2 *b) {
  float _11 = a->_11 * b->_11 + a->_12 * b->_21;
  float _12 = a->_11 * b->_12 + a->_12 * b->_22;

  float _21 = a->_21 * b->_11 + a->_22 * b->_21;
  float _22 = a->_21 * b->_12 + a->_22 * b->_22;

  return mat2_new(_11, _12, _21, _22);
}

mat2 mat2_transpose(const mat2 *m2) {
  /* *
   * |a b| -> |a c|
   * |c d|    |b d|
   * */
  return mat2_new(m2->_11, m2->_21, m2->_12, m2->_22);
}

float mat2_cut(const mat2 *m2, int row, int col) {
  if (row == 1 and col == 1)
    return m2->_11;
  if (row == 1 and col == 0)
    return m2->_12;
  if (row == 0 and col == 1)
    return m2->_21;
  if (row == 0 and col == 0)
    return m2->_22;

  return -1.0;
}

mat2 mat2_cofactor(const mat2 *m2) {
  /* *
   * | + - |
   * | - + |
   * */
  return mat2_new(m2->_11, -m2->_12, -m2->_21, m2->_22);
}

mat2 mat2_minor(const mat2 *m2) {
  return mat2_new(m2->_22, m2->_12, m2->_21, m2->_11);
}

mat2 mat2_adjugate(const mat2 *m2) {
  mat2 cof = mat2_cofactor(m2);
  return mat2_transpose(&cof);
}

float mat2_determinant(const mat2 *m2) {
  float a = m2->_11;
  float b = m2->_12;

  float c = m2->_21;
  float d = m2->_22;

  return (a * d) - (b * c);
}

bool mat2_inverse(mat2 *dest, const mat2 *m2) {
  float m2_det = mat2_determinant(m2);
  if (m2_det == 0.0f) {
    return false;
  }
  float det = 1.0f / m2_det;

  dest->_11 = (m2->_22) * det;
  dest->_12 = -(m2->_12) * det;

  dest->_21 = -(m2->_21) * det;
  dest->_22 = (m2->_11) * det;

  return true;
}

mat2 mat2_rotate_z(const mat2 *m2, float by) {
  /* *
   * |  cos sin |
   * | -sin cos |
   * */
  by = DEG2RAD(by);
  float c = cosf(by);
  float s = sinf(by);
  mat2 m = mat2_new(c, s, -s, c);
  return mat2_mul(m2, &m);
}

void mat2_print(const mat2 *m2) {
  printf("%.3f, %.3f \n", m2->_11, m2->_12);
  printf("%.3f, %.3f \n", m2->_21, m2->_22);
}

/* MAT 3 */

mat3 mat3_zero() {
  mat3 m3 = (mat3){._11 = 0.0f,
                   ._12 = 0.0f,
                   ._13 = 0.0f,
                   ._21 = 0.0f,
                   ._22 = 0.0f,
                   ._23 = 0.0f,
                   ._31 = 0.0f,
                   ._32 = 0.0f,
                   ._33 = 0.0f};
  return m3;
}

mat3 mat3_identity() {
  mat3 m3 = (mat3){._11 = 1.0f,
                   ._12 = 0.0f,
                   ._13 = 0.0f,
                   ._21 = 0.0f,
                   ._22 = 1.0f,
                   ._23 = 0.0f,
                   ._31 = 0.0f,
                   ._32 = 0.0f,
                   ._33 = 1.0f};
  return m3;
}

mat3 mat3_new(float _11, float _12, float _13, float _21, float _22, float _23,
              float _31, float _32, float _33) {
  mat3 m3 = (mat3){._11 = _11,
                   ._12 = _12,
                   ._13 = _13,
                   ._21 = _21,
                   ._22 = _22,
                   ._23 = _23,
                   ._31 = _31,
                   ._32 = _32,
                   ._33 = _33};
  return m3;
}

void mat3_cpy(mat3 *dest, const mat3 *m3) {
  dest->_11 = m3->_11;
  dest->_12 = m3->_12;
  dest->_13 = m3->_13;

  dest->_21 = m3->_21;
  dest->_22 = m3->_22;
  dest->_23 = m3->_23;

  dest->_31 = m3->_31;
  dest->_32 = m3->_32;
  dest->_33 = m3->_33;
}

mat3 mat3_add(const mat3 *a, const mat3 *b) {
  return mat3_new(a->_11 + b->_11, a->_12 + b->_12, a->_13 + b->_13,
                  a->_21 + b->_21, a->_22 + b->_22, a->_23 + b->_23,
                  a->_31 + b->_31, a->_32 + b->_32, a->_33 + b->_33);
}

mat3 mat3_sub(const mat3 *a, const mat3 *b) {
  return mat3_new(a->_11 - b->_11, a->_12 - b->_12, a->_13 - b->_13,
                  a->_21 - b->_21, a->_22 - b->_22, a->_23 - b->_23,
                  a->_31 - b->_31, a->_32 - b->_32, a->_33 - b->_33);
}

mat3 mat3_scale(const mat3 *m3, float by) {
  return mat3_new(m3->_11 * by, m3->_12 * by, m3->_13 * by, m3->_21 * by,
                  m3->_22 * by, m3->_23 * by, m3->_31 * by, m3->_32 * by,
                  m3->_33 * by);
}

mat3 mat3_mul(const mat3 *a, const mat3 *b) {
  /* *
   * | a b c |   | a b c |   | (aa+bd+cg) (ab+be+ch) (ac+bf+ci) |
   * | d e f | x | d e f | = | (da+ed+fg) (db+ee+fh) (dc+ef+fi) |
   * | g h i |   | g h i |   | (ga+hd+ig) (gb+he+ih) (gc+hf+ii) |
   * */
  float _11 = a->_11 * b->_11 + a->_12 * b->_21 + a->_13 * b->_31;
  float _12 = a->_11 * b->_12 + a->_12 * b->_22 + a->_13 * b->_32;
  float _13 = a->_11 * b->_13 + a->_12 * b->_23 + a->_13 * b->_33;

  float _21 = a->_21 * b->_11 + a->_22 * b->_21 + a->_23 * b->_31;
  float _22 = a->_21 * b->_12 + a->_22 * b->_22 + a->_23 * b->_32;
  float _23 = a->_21 * b->_13 + a->_22 * b->_23 + a->_23 * b->_33;

  float _31 = a->_31 * b->_11 + a->_32 * b->_21 + a->_33 * b->_31;
  float _32 = a->_31 * b->_12 + a->_32 * b->_22 + a->_33 * b->_32;
  float _33 = a->_31 * b->_13 + a->_32 * b->_23 + a->_33 * b->_33;

  return mat3_new(_11, _12, _13, _21, _22, _23, _31, _32, _33);
}

mat3 mat3_transpose(const struct mat3 *m3) {
  /* *
   * | a b c |    | a d g |
   * | d e f | -> | b e h |
   * | g h i |    | c f i |
   * */
  return mat3_new(m3->_11, m3->_21, m3->_31, m3->_12, m3->_22, m3->_32, m3->_13,
                  m3->_23, m3->_33);
}

mat3 mat3_cofactor(const mat3 *m3) {
  /* *
   * | + - + |
   * | - + - |
   * | + - + |
   * */
  return mat3_new(m3->_11, -m3->_12, m3->_13, -m3->_21, m3->_22, -m3->_23,
                  m3->_31, -m3->_32, m3->_33);
}

mat2 mat3_cut(const mat3 *m3, int row, int col) {
  /* *
   * r = 0
   * c = 0
   *     0 1 2
   * 0 | a b c |    | -   -  - |
   * 1 | d e f | -> | - | e  f |
   * 2 | g h i |    | - | h  i |
   * */
  mat2 result;

  int n = 0;
  for (int i = 0; i < 3; ++i) {
    for (int j = 0; j < 3; ++j) {
      if (i == row or j == col) {
        continue;
      }
      int id = n;
      ++n;
      result.arr[id] = m3->arr[3 * i + j];
    }
  }

  return result;
}

mat3 mat3_minor(const mat3 *m3) {
  /* *
   * | a b c |    | (ei - fh) (di - fg) (dh - eg) |
   * | d e f | -> | (bi - ch) (ai - cg) (ah - bg) |
   * | g h i |    | (bf - ce) (af - cd) (ae - bd) |
   * */
  mat3 result;
  for (int i = 0; i < 3; ++i) {
    for (int j = 0; j < 3; ++j) {
      mat2 cut = mat3_cut(m3, i, j);
      int id = 3 * i + j;
      result.arr[id] = mat2_determinant(&cut);
    }
  }
  return result;
}

mat3 mat3_adjugate(const mat3 *m3) {
  mat3 cof = mat3_cofactor(m3);
  return mat3_transpose(&cof);
}

float mat3_determinant(const mat3 *m3) {
  float a = m3->_11;
  float b = m3->_12;
  float c = m3->_13;

  float d = m3->_21;
  float e = m3->_22;
  float f = m3->_23;

  float g = m3->_31;
  float h = m3->_32;
  float i = m3->_33;
  return a * (e * i - f * h) - b * (d * i - f * g) + c * (d * h - e * g);
}

bool mat3_inverse(mat3 *dest, const mat3 *m3) {
  float m3_det = mat3_determinant(m3);
  if (m3_det == 0.0f) {
    return false;
  }
  float det = 1.0f / m3_det;

  dest->_11 = (m3->_22 * m3->_33 - m3->_23 * m3->_32) * det;
  dest->_12 = (m3->_13 * m3->_32 - m3->_12 * m3->_33) * det;
  dest->_13 = (m3->_12 * m3->_23 - m3->_13 * m3->_22) * det;

  dest->_21 = (m3->_23 * m3->_31 - m3->_21 * m3->_33) * det;
  dest->_22 = (m3->_11 * m3->_33 - m3->_13 * m3->_31) * det;
  dest->_23 = (m3->_13 * m3->_21 - m3->_11 * m3->_23) * det;

  dest->_31 = (m3->_21 * m3->_32 - m3->_22 * m3->_31) * det;
  dest->_32 = (m3->_12 * m3->_31 - m3->_11 * m3->_32) * det;
  dest->_33 = (m3->_11 * m3->_22 - m3->_12 * m3->_21) * det;

  return true;
}

mat3 mat3_scaling(const mat3 *m3, const vec3 *v3) {
  /* *
   * | x  0  0 |
   * | 0  y  0 |
   * | 0  0  z |
   * */
  mat3 m = mat3_zero();
  m._11 = v3->x;
  m._22 = v3->y;
  m._33 = v3->z;

  return mat3_mul(m3, &m);
}

mat3 mat3_rotate_x(const mat3 *m3, float by) {
  /* *
   * | 1    0    0 |
   * | 0  cos  sin |
   * | 0 -sin  cos |
   * */
  by = DEG2RAD(by);
  float c = cosf(by);
  float s = sinf(by);

  mat3 m = mat3_identity();
  m._22 = c;
  m._23 = s;
  m._32 = -s;
  m._33 = c;

  return mat3_mul(m3, &m);
}

mat3 mat3_rotate_y(const mat3 *m3, float by) {
  /* *
   * | cos  0  -sin |
   * |   0  1     0 |
   * | sin  0   cos |
   * */
  by = DEG2RAD(by);
  float c = cosf(by);
  float s = sinf(by);

  mat3 m = mat3_identity();
  m._11 = c;
  m._13 = -s;
  m._31 = s;
  m._33 = c;

  return mat3_mul(m3, &m);
}

mat3 mat3_rotate_z(const mat3 *m3, float by) {
  /* *
   * |  cos  sin    0 |
   * | -sin  cos    0 |
   * |    0    0    1 |
   * */
  by = DEG2RAD(by);
  float c = cosf(by);
  float s = sinf(by);
  mat3 m = mat3_identity();

  m._11 = c;
  m._12 = s;
  m._21 = -s;
  m._22 = c;

  return mat3_mul(m3, &m);
}

void mat3_print(const mat3 *m3) {
  printf("%.3f, %.3f, %.3f \n", m3->_11, m3->_12, m3->_13);
  printf("%.3f, %.3f, %.3f \n", m3->_21, m3->_22, m3->_23);
  printf("%.3f, %.3f, %.3f \n", m3->_31, m3->_32, m3->_33);
}

/* MAT 4 */

mat4 mat4_zero() {
  mat4 m4 = (mat4){._11 = 0.0f,
                   ._12 = 0.0f,
                   ._13 = 0.0f,
                   ._14 = 0.0f,
                   ._21 = 0.0f,
                   ._22 = 0.0f,
                   ._23 = 0.0f,
                   ._24 = 0.0f,
                   ._31 = 0.0f,
                   ._32 = 0.0f,
                   ._33 = 0.0f,
                   ._34 = 0.0f,
                   ._41 = 0.0f,
                   ._42 = 0.0f,
                   ._43 = 0.0f,
                   ._44 = 0.0f};
  return m4;
}

mat4 mat4_identity() {
  mat4 m4 = (mat4){._11 = 1.0f,
                   ._12 = 0.0f,
                   ._13 = 0.0f,
                   ._14 = 0.0f,
                   ._21 = 0.0f,
                   ._22 = 1.0f,
                   ._23 = 0.0f,
                   ._24 = 0.0f,
                   ._31 = 0.0f,
                   ._32 = 0.0f,
                   ._33 = 1.0f,
                   ._34 = 0.0f,
                   ._41 = 0.0f,
                   ._42 = 0.0f,
                   ._43 = 0.0f,
                   ._44 = 1.0f};
  return m4;
}

mat4 mat4_new(float _11, float _12, float _13, float _14, float _21, float _22,
              float _23, float _24, float _31, float _32, float _33, float _34,
              float _41, float _42, float _43, float _44) {
  mat4 m4 = (mat4){._11 = _11,
                   ._12 = _12,
                   ._13 = _13,
                   ._14 = _14,
                   ._21 = _21,
                   ._22 = _22,
                   ._23 = _23,
                   ._24 = _24,
                   ._31 = _31,
                   ._32 = _32,
                   ._33 = _33,
                   ._34 = _34,
                   ._41 = _41,
                   ._42 = _42,
                   ._43 = _43,
                   ._44 = _44};
  return m4;
}

void mat4_cpy(mat4 *dest, const mat4 *m4) {
  dest->_11 = m4->_11;
  dest->_12 = m4->_12;
  dest->_13 = m4->_13;
  dest->_14 = m4->_14;

  dest->_21 = m4->_21;
  dest->_22 = m4->_22;
  dest->_23 = m4->_23;
  dest->_24 = m4->_24;

  dest->_31 = m4->_31;
  dest->_32 = m4->_32;
  dest->_33 = m4->_33;
  dest->_34 = m4->_34;

  dest->_41 = m4->_41;
  dest->_42 = m4->_42;
  dest->_43 = m4->_43;
  dest->_44 = m4->_44;
}

mat4 mat4_add(const mat4 *a, const mat4 *b) {
  return mat4_new(
      a->_11 + b->_11, a->_12 + b->_12, a->_13 + b->_13, a->_14 + b->_14,
      a->_21 + b->_21, a->_22 + b->_22, a->_23 + b->_23, a->_24 + b->_24,
      a->_31 + b->_31, a->_32 + b->_32, a->_33 + b->_33, a->_34 + b->_34,
      a->_41 + b->_41, a->_42 + b->_42, a->_43 + b->_43, a->_44 + b->_44);
}

mat4 mat4_sub(const mat4 *a, const mat4 *b) {
  return mat4_new(
      a->_11 - b->_11, a->_12 - b->_12, a->_13 - b->_13, a->_14 - b->_14,
      a->_21 - b->_21, a->_22 - b->_22, a->_23 - b->_23, a->_24 - b->_24,
      a->_31 - b->_31, a->_32 - b->_32, a->_33 - b->_33, a->_34 - b->_34,
      a->_41 - b->_41, a->_42 - b->_42, a->_43 - b->_43, a->_44 - b->_44);
}

mat4 mat4_scale(const mat4 *m4, float by) {
  return mat4_new(m4->_11 * by, m4->_12 * by, m4->_13 * by, m4->_14 * by,
                  m4->_21 * by, m4->_22 * by, m4->_23 * by, m4->_24 * by,
                  m4->_31 * by, m4->_32 * by, m4->_33 * by, m4->_34 * by,
                  m4->_41 * by, m4->_42 * by, m4->_43 * by, m4->_44 * by);
}

mat4 mat4_mul(const mat4 *a, const mat4 *b) {
  float _11 =
      a->_11 * b->_11 + a->_12 * b->_21 + a->_13 * b->_31 + a->_14 * b->_41;
  float _12 =
      a->_11 * b->_12 + a->_12 * b->_22 + a->_13 * b->_32 + a->_14 * b->_42;
  float _13 =
      a->_11 * b->_13 + a->_12 * b->_23 + a->_13 * b->_33 + a->_14 * b->_43;
  float _14 =
      a->_11 * b->_14 + a->_12 * b->_24 + a->_13 * b->_34 + a->_14 * b->_44;
  float _21 =
      a->_21 * b->_11 + a->_22 * b->_21 + a->_23 * b->_31 + a->_24 * b->_41;
  float _22 =
      a->_21 * b->_12 + a->_22 * b->_22 + a->_23 * b->_32 + a->_24 * b->_42;
  float _23 =
      a->_21 * b->_13 + a->_22 * b->_23 + a->_23 * b->_33 + a->_24 * b->_43;
  float _24 =
      a->_21 * b->_14 + a->_22 * b->_24 + a->_23 * b->_34 + a->_24 * b->_44;
  float _31 =
      a->_31 * b->_11 + a->_32 * b->_21 + a->_33 * b->_31 + a->_34 * b->_41;
  float _32 =
      a->_31 * b->_12 + a->_32 * b->_22 + a->_33 * b->_32 + a->_34 * b->_42;
  float _33 =
      a->_31 * b->_13 + a->_32 * b->_23 + a->_33 * b->_33 + a->_34 * b->_43;
  float _34 =
      a->_31 * b->_14 + a->_32 * b->_24 + a->_33 * b->_34 + a->_34 * b->_44;
  float _41 =
      a->_41 * b->_11 + a->_42 * b->_21 + a->_43 * b->_31 + a->_44 * b->_41;
  float _42 =
      a->_41 * b->_12 + a->_42 * b->_22 + a->_43 * b->_32 + a->_44 * b->_42;
  float _43 =
      a->_41 * b->_13 + a->_42 * b->_23 + a->_43 * b->_33 + a->_44 * b->_43;
  float _44 =
      a->_41 * b->_14 + a->_42 * b->_24 + a->_43 * b->_34 + a->_44 * b->_44;
  return mat4_new(_11, _12, _13, _14, _21, _22, _23, _24, _31, _32, _33, _34,
                  _41, _42, _43, _44);
}

mat4 mat4_transpose(const mat4 *m4) {
  /* *
   * | a b c d |    | a e i m |
   * | e f g h | -> | b f j n |
   * | i j k l |    | c g k o |
   * | m n o p |    | d h l p |
   * */
  return mat4_new(m4->_11, m4->_21, m4->_31, m4->_41, m4->_12, m4->_22, m4->_32,
                  m4->_42, m4->_13, m4->_23, m4->_33, m4->_43, m4->_14, m4->_24,
                  m4->_34, m4->_44);
}

mat4 mat4_cofactor(const mat4 *m4) {
  /* *
   * |  +  -  +  - |
   * |  -  +  -  + |
   * |  +  -  +  - |
   * |  -  +  -  + |
   * */
  return mat4_new(m4->_11, -m4->_12, m4->_13, -m4->_14, -m4->_21, m4->_22,
                  -m4->_23, m4->_24, m4->_31, -m4->_32, m4->_33, -m4->_34,
                  -m4->_41, m4->_42, -m4->_43, m4->_44);
}

mat3 mat4_cut(const mat4 *m4, int row, int col) {
  /* *
   * r = 0
   * c = 0
   *     0 1 2 3
   * 0 | a b c d |    | -   -  -  - |
   * 1 | e f g h | -> | - | f  g  h |
   * 2 | i j k l |    | - | j  k  l |
   * 3 | m n o p |    | - | n  o  p |
   * */
  mat3 result;
  int n = 0;
  for (int i = 0; i < 4; ++i) {
    for (int j = 0; j < 4; ++j) {
      if (i == row or j == col) {
        continue;
      }
      int t = n;
      ++n;
      result.arr[t] = m4->arr[4 * i + j];
    }
  }
  return result;
}

mat4 mat4_minor(const mat4 *m4) {
  mat4 result;
  for (int i = 0; i < 4; ++i) {
    for (int j = 0; j < 4; ++j) {
      mat3 cut = mat4_cut(m4, i, j);
      result.arr[4 * i + j] = mat3_determinant(&cut);
    }
  }
  return result;
}

mat4 mat4_adjugate(const mat4 *m4) {
  mat4 cof = mat4_cofactor(m4);
  return mat4_transpose(&cof);
}

float mat4_determinant(const mat4 *m4) {
  float a = m4->_11;
  float b = m4->_12;
  float c = m4->_13;
  float d = m4->_14;
  float e = m4->_21;
  float f = m4->_22;
  float g = m4->_23;
  float h = m4->_24;
  float i = m4->_31;
  float j = m4->_32;
  float k = m4->_33;
  float l = m4->_34;
  float m = m4->_41;
  float n = m4->_42;
  float o = m4->_43;
  float p = m4->_44;

  return d * g * j * m - c * h * j * m - d * f * k * m + b * h * k * m +
         c * f * l * m - b * g * l * m - d * g * i * n + c * h * i * n +
         d * e * k * n - a * h * k * n - c * e * l * n + a * g * l * n +
         d * f * i * o - b * h * i * o - d * e * j * o + a * h * j * o +
         b * e * l * o - a * f * l * o - c * f * i * o + b * g * i * p +
         c * e * j * p - a * g * j * p - b * e * k * p + a * f * k * p;
}

bool mat4_inverse(mat4 *dest, const mat4 *m4) {
  float m4_det = mat4_determinant(m4);
  if (m4_det == 0.0f) {
    return false;
  }
  float det = 1.0f / m4_det;

  dest->_11 = (m4->_22 * m4->_33 * m4->_44 + m4->_23 * m4->_34 * m4->_42 +
               m4->_24 * m4->_32 * m4->_43 - m4->_22 * m4->_34 * m4->_43 -
               m4->_23 * m4->_32 * m4->_44 - m4->_24 * m4->_33 * m4->_42) *
              det;
  dest->_12 = (m4->_12 * m4->_34 * m4->_43 + m4->_13 * m4->_32 * m4->_44 +
               m4->_14 * m4->_33 * m4->_42 - m4->_12 * m4->_33 * m4->_44 -
               m4->_13 * m4->_34 * m4->_42 - m4->_14 * m4->_32 * m4->_43) *
              det;
  dest->_13 = (m4->_12 * m4->_23 * m4->_44 + m4->_13 * m4->_24 * m4->_42 +
               m4->_14 * m4->_22 * m4->_43 - m4->_12 * m4->_24 * m4->_43 -
               m4->_13 * m4->_22 * m4->_44 - m4->_14 * m4->_23 * m4->_42) *
              det;
  dest->_14 = (m4->_12 * m4->_24 * m4->_33 + m4->_13 * m4->_22 * m4->_34 +
               m4->_14 * m4->_23 * m4->_32 - m4->_12 * m4->_23 * m4->_34 -
               m4->_13 * m4->_24 * m4->_32 - m4->_14 * m4->_22 * m4->_33) *
              det;
  dest->_21 = (m4->_21 * m4->_34 * m4->_43 + m4->_23 * m4->_31 * m4->_44 +
               m4->_24 * m4->_33 * m4->_41 - m4->_21 * m4->_33 * m4->_44 -
               m4->_23 * m4->_34 * m4->_41 - m4->_24 * m4->_31 * m4->_43) *
              det;
  dest->_22 = (m4->_11 * m4->_33 * m4->_44 + m4->_13 * m4->_34 * m4->_41 +
               m4->_14 * m4->_31 * m4->_43 - m4->_11 * m4->_34 * m4->_43 -
               m4->_13 * m4->_31 * m4->_44 - m4->_14 * m4->_33 * m4->_41) *
              det;
  dest->_23 = (m4->_11 * m4->_24 * m4->_43 + m4->_13 * m4->_21 * m4->_44 +
               m4->_14 * m4->_23 * m4->_41 - m4->_11 * m4->_23 * m4->_44 -
               m4->_13 * m4->_24 * m4->_41 - m4->_14 * m4->_21 * m4->_43) *
              det;
  dest->_24 = (m4->_11 * m4->_23 * m4->_34 + m4->_13 * m4->_24 * m4->_31 +
               m4->_14 * m4->_21 * m4->_33 - m4->_11 * m4->_24 * m4->_33 -
               m4->_13 * m4->_21 * m4->_34 - m4->_14 * m4->_23 * m4->_31) *
              det;
  dest->_31 = (m4->_21 * m4->_32 * m4->_44 + m4->_22 * m4->_34 * m4->_41 +
               m4->_24 * m4->_31 * m4->_42 - m4->_21 * m4->_34 * m4->_42 -
               m4->_22 * m4->_31 * m4->_44 - m4->_24 * m4->_32 * m4->_41) *
              det;
  dest->_32 = (m4->_11 * m4->_34 * m4->_42 + m4->_12 * m4->_31 * m4->_44 +
               m4->_14 * m4->_32 * m4->_41 - m4->_11 * m4->_32 * m4->_44 -
               m4->_12 * m4->_34 * m4->_41 - m4->_14 * m4->_31 * m4->_42) *
              det;
  dest->_33 = (m4->_11 * m4->_22 * m4->_44 + m4->_12 * m4->_24 * m4->_41 +
               m4->_14 * m4->_21 * m4->_42 - m4->_11 * m4->_24 * m4->_42 -
               m4->_12 * m4->_21 * m4->_44 - m4->_14 * m4->_22 * m4->_41) *
              det;
  dest->_34 = (m4->_11 * m4->_24 * m4->_32 + m4->_12 * m4->_21 * m4->_34 +
               m4->_14 * m4->_22 * m4->_31 - m4->_11 * m4->_22 * m4->_34 -
               m4->_12 * m4->_24 * m4->_31 - m4->_14 * m4->_21 * m4->_32) *
              det;
  dest->_41 = (m4->_21 * m4->_33 * m4->_42 + m4->_22 * m4->_31 * m4->_43 +
               m4->_23 * m4->_32 * m4->_41 - m4->_21 * m4->_32 * m4->_43 -
               m4->_22 * m4->_33 * m4->_41 - m4->_23 * m4->_31 * m4->_42) *
              det;
  dest->_42 = (m4->_11 * m4->_32 * m4->_43 + m4->_12 * m4->_33 * m4->_41 +
               m4->_13 * m4->_31 * m4->_42 - m4->_11 * m4->_33 * m4->_42 -
               m4->_12 * m4->_31 * m4->_43 - m4->_13 * m4->_32 * m4->_41) *
              det;
  dest->_43 = (m4->_11 * m4->_23 * m4->_42 + m4->_12 * m4->_21 * m4->_43 +
               m4->_13 * m4->_22 * m4->_41 - m4->_11 * m4->_22 * m4->_43 -
               m4->_12 * m4->_23 * m4->_41 - m4->_13 * m4->_21 * m4->_42) *
              det;
  dest->_44 = (m4->_11 * m4->_22 * m4->_33 + m4->_12 * m4->_23 * m4->_31 +
               m4->_13 * m4->_21 * m4->_32 - m4->_11 * m4->_23 * m4->_32 -
               m4->_12 * m4->_21 * m4->_33 - m4->_13 * m4->_22 * m4->_31) *
              det;

  return true;
}

mat4 mat4_translation(const mat4 *m4, const vec3 *v3) {
  /* *
   * row order
   * |  1  0  0  0 |
   * |  0  1  0  0 |
   * |  0  0  1  0 |
   * ---------------
   * |  x  y  z  1 |
   * */
  mat4 m = mat4_identity();

  m._41 = v3->x;
  m._42 = v3->y;
  m._43 = v3->z;

  return mat4_mul(m4, &m);
}

mat4 mat4_scaling(const mat4 *m4, const vec3 *v3) {
  /* *
   * | x  0  0  0 |
   * | 0  y  0  0 |
   * | 0  0  z  0 |
   * | 0  0  0  1 |
   * */
  mat4 m = mat4_identity();

  m._11 = v3->x;
  m._22 = v3->y;
  m._33 = v3->z;

  return mat4_mul(m4, &m);
}

mat4 mat4_rotate_x(const mat4 *m4, float by) {
  /* *
   * | 1    0     0  0 |
   * | 0  cos  -sin  0 |
   * | 0  sin   cos  0 |
   * | 0    0     0  1 |
   * */
  by = DEG2RAD(by);
  float c = cosf(by);
  float s = sinf(by);
  mat4 m = mat4_identity();

  m._22 = c;
  m._23 = -s;
  m._32 = s;
  m._33 = c;

  return mat4_mul(m4, &m);
}

mat4 mat4_rotate_y(const mat4 *m4, float by) {
  /* *
   * |  cos  0  sin  0 |
   * |    0  1    0  0 |
   * | -sin  0  cos  0 |
   * |    0  0    0  1 |
   * */
  by = DEG2RAD(by);
  float c = cosf(by);
  float s = sinf(by);
  mat4 m = mat4_identity();

  m._11 = c;
  m._13 = s;
  m._31 = -s;
  m._33 = c;

  return mat4_mul(m4, &m);
}

mat4 mat4_rotate_z(const mat4 *m4, float by) {
  /* *
   * |  cos  sin  0  0 |
   * | -sin  cos  0  0 |
   * |    0    0  1  0 |
   * |    0    0  0  1 |
   * */
  by = DEG2RAD(by);
  float c = cosf(by);
  float s = sinf(by);
  mat4 m = mat4_identity();

  m._11 = c;
  m._12 = s;
  m._21 = -s;
  m._22 = c;

  return mat4_mul(m4, &m);
}

// TODO
internal mat4 mat4_rotate(const mat4 *m4, float by, const vec3 *v3) {
  /* https://en.wikipedia.org/wiki/Rotation_matrix
   * c = cos
   * s = sin
   * t = 1-c
   * | (c + x*x*t),   (x*y*t - z*s), (x*z*t + y*s) |  |
   * | (y*x*t + z*s), (c + y*y*t),   (y*z*t - x*s) |  |
   * | (z*x*t - y*s), (z*y*t + x*s), (c + z*z*t)   |  |
   * --------------------------------------------------
   * |                                             |  |
   */
  by = DEG2RAD(by);
  float x = v3->x;
  float y = v3->y;
  float z = v3->z;
  float c = cosf(by);
  float s = sinf(by);
  float t = 1.0 - c;

  mat4 m = mat4_identity();
  m._11 = c + (x*x*t);
  m._12 = (x*y*t) - (z*s);
  m._13 = (x*z*t) + (y*s);

  m._21 = (y*x*t) + (z*s);
  m._22 = c + (y*y*t);
  m._23 = (y*z*t) - (x*s);

  m._31 = (z*x*t) - (y*s);
  m._32 = (z*y*t) + (x*s);
  m._33 = c + (z*z*t);


  return mat4_identity();
}

void mat4_print(const mat4 *m4) {
  printf("%.3f, %.3f, %.3f, %.3f \n", m4->_11, m4->_12, m4->_13, m4->_14);
  printf("%.3f, %.3f, %.3f, %.3f \n", m4->_21, m4->_22, m4->_23, m4->_24);
  printf("%.3f, %.3f, %.3f, %.3f \n", m4->_31, m4->_32, m4->_33, m4->_34);
  printf("%.3f, %.3f, %.3f, %.3f \n", m4->_41, m4->_42, m4->_43, m4->_44);
}
