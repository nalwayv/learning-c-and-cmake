/* *
 * matrix
 * */
#include <math.h>
#include <stdio.h>

#include "matrix.h"
#include "utils.h"

/* HELPERS */

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
  for (int i = 0; i < 4; ++i) {
    dest->arr[i] = m2->arr[i];
  }
}

void mat2_add(mat2 *dest, const mat2 *a, const mat2 *b) {
  for (int i = 0; i < 4; ++i) {
    dest->arr[i] = a->arr[i] + b->arr[i];
  }
}

void mat2_sub(mat2 *dest, const mat2 *a, const mat2 *b) {
  for (int i = 0; i < 4; ++i) {
    dest->arr[i] = a->arr[i] - b->arr[i];
  }
}

void mat2_scale(mat2 *dest, const mat2 *m2, float by) {
  for (int i = 0; i < 4; ++i) {
    dest->arr[i] = m2->arr[i] * by;
  }
}

void mat2_mul(mat2 *dest, const mat2 *a, const mat2 *b) {
  dest->_11 = a->_11 * b->_11 + a->_12 * b->_21;
  dest->_12 = a->_11 * b->_12 + a->_12 * b->_22;
  dest->_21 = a->_21 * b->_11 + a->_22 * b->_21;
  dest->_22 = a->_21 * b->_12 + a->_22 * b->_22;
}

void mat2_transpose(mat2 *dest, const mat2 *m2) {
  /* *
   * |a b| -> |a c|
   * |c d|    |b d|
   * */
  dest->_11 = m2->_11;
  dest->_12 = m2->_21;
  dest->_21 = m2->_12;
  dest->_22 = m2->_22;
}

float mat2_cut(const mat2 *m2, int row, int col){
  if(row == 1 and col == 1) return m2->_11;
  if(row == 1 and col == 0) return m2->_12;
  if(row == 0 and col == 1) return m2->_21;
  if(row == 0 and col == 0) return m2->_22;

  return -1.0;
}

void mat2_cofactor(mat2 *dest, const mat2 *m2) {
  /* *
   * | + - |
   * | - + |
   * */
  dest->_11 = m2->_11;
  dest->_12 = -m2->_12;
  dest->_21 = -m2->_21;
  dest->_22 = m2->_22;
}

void mat2_minor(mat2 *dest, const mat2 *m2) {
  dest->_11 = m2->_22;
  dest->_12 = m2->_12;
  dest->_21 = m2->_21;
  dest->_22 = m2->_11;
}

void mat2_adjugate(mat2 *dest, const mat2 *m2) {
  mat2 cof;
  mat2_cofactor(&cof ,m2);
  mat2_transpose(dest, &cof);
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

void mat2_rotate_z(mat2 *dest, float by) {
  /* *
   * |  cos sin |
   * | -sin cos |
   * */
  by = DEG2RAD(by);

  mat2 i = mat2_identity();
  i._11 = cosf(by);
  i._12 = sinf(by);
  i._21 = -sinf(by);
  i._22 = cosf(by);
  mat2_cpy(dest, &i);
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
  for (int i = 0; i < 9; ++i) {
    dest->arr[i] = m3->arr[i];
  }
}

void mat3_add(mat3 *dest, const mat3 *a, const mat3 *b) {
  for (int i = 0; i < 9; ++i) {
    dest->arr[i] = a->arr[i] + b->arr[i];
  }
}

void mat3_sub(mat3 *dest, const mat3 *a, const mat3 *b) {
  for (int i = 0; i < 9; ++i) {
    dest->arr[i] = a->arr[i] - b->arr[i];
  }
}

void mat3_scale(mat3 *dest, const mat3 *m3, float by) {
  for (int i = 0; i < 9; ++i) {
    dest->arr[i] = m3->arr[i] * by;
  }
}

void mat3_mul(mat3 *dest, const mat3 *a, const mat3 *b) {
  /* *
   * | a b c |   | a b c |   | (aa+bd+cg) (ab+be+ch) (ac+bf+ci) |
   * | d e f | x | d e f | = | (da+ed+fg) (db+ee+fh) (dc+ef+fi) |
   * | g h i |   | g h i |   | (ga+hd+ig) (gb+he+ih) (gc+hf+ii) |
   * */
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

void mat3_transpose(mat3 *dest, const struct mat3 *m3) {
  /* *
   * | a b c |    | a d g |
   * | d e f | -> | b e h |
   * | g h i |    | c f i |
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

void mat3_cofactor(mat3 *dest, const struct mat3 *m3) {
  /* *
   * | + - + |
   * | - + - |
   * | + - + |
   * */
  dest->_11 = m3->_11;
  dest->_12 = -m3->_12;
  dest->_13 = m3->_13;
  dest->_21 = -m3->_21;
  dest->_22 = m3->_22;
  dest->_23 = -m3->_23;
  dest->_31 = m3->_31;
  dest->_32 = -m3->_32;
  dest->_33 = m3->_33;
}

void mat3_cut(mat2 *dest, const mat3 *m3, int r, int c) {
  /* *
   * r = 0
   * c = 0
   *     0 1 2
   * 0 | a b c |    | -   -  - |
   * 1 | d e f | -> | - | e  f |
   * 2 | g h i |    | - | h  i |
   * */
  int n = 0;
  for (int i = 0; i < 3; ++i) {
    for (int j = 0; j < 3; ++j) {
      if (i == r or j == c) {
        continue;
      }
      int t = n;
      ++n;
      dest->arr[t] = m3->arr[3 * i + j];
    }
  }
}

void mat3_minor(mat3 *dest, const mat3 *m3) {
  /* *
   * | a b c |    | (ei - fh) (di - fg) (dh - eg) |
   * | d e f | -> | (bi - ch) (ai - cg) (ah - bg) |
   * | g h i |    | (bf - ce) (af - cd) (ae - bd) |
   * */
  for (int i = 0; i < 3; ++i) {
    for (int j = 0; j < 3; ++j) {
      mat2 m2;
      mat3_cut(&m2, m3, i, j);
      int id = 3 * i + j;
      dest->arr[id] = mat2_determinant(&m2);
    }
  }
}

void mat3_adjugate(mat3 *dest, const mat3 *m3) {
  mat3 cof;
  mat3_cofactor(&cof, m3);
  mat3_transpose(dest, &cof);
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

void mat3_scaling(mat3 *dest, float x, float y, float z) {
  /* *
   * | x  0  0 |
   * | 0  y  0 |
   * | 0  0  z |
   * */
  mat3 i = mat3_identity();
  i._11 = x;
  i._22 = y;
  i._33 = z;
  mat3_cpy(dest, &i);
}

void mat3_rotate_x(mat3 *dest, float by) {
  /* *
   * | 1    0    0 |
   * | 0  cos  sin |
   * | 0 -sin  cos |
   * */
  by = DEG2RAD(by);

  mat3 i = mat3_identity();
  i._22 = cosf(by);
  i._23 = sinf(by);
  i._32 = -sinf(by);
  i._33 = cosf(by);
  mat3_cpy(dest, &i);
}

void mat3_rotate_y(mat3 *dest, float by) {
  /* *
   * | cos  0  -sin |
   * |   0  1     0 |
   * | sin  0   cos |
   * */
  by = DEG2RAD(by);

  mat3 i = mat3_identity();
  i._11 = cosf(by);
  i._13 = -sinf(by);
  i._31 = sinf(by);
  i._33 = cosf(by);
  mat3_cpy(dest, &i);
}

void mat3_rotate_z(mat3 *dest, float by) {
  /* *
   * |  cos  sin    0 |
   * | -sin  cos    0 |
   * |    0    0    1 |
   * */
  by = DEG2RAD(by);

  mat3 i = mat3_identity();
  i._11 = cosf(by);
  i._12 = sinf(by);
  i._21 = -sinf(by);
  i._22 = cosf(by);
  mat3_cpy(dest, &i);
}

void mat3_rotation(mat3 *dest, float pitch, float yaw, float roll) {
  mat3 z;
  mat3_rotate_z(&z, roll);
  mat3 x;
  mat3_rotate_x(&x, pitch);
  mat3 y;
  mat3_rotate_y(&y, yaw);

  mat3 xy;
  mat3_mul(&xy, &x ,&y);

  mat3_mul(dest, &z, &xy);
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

void mat4_cpy(mat4 *dest, const mat4 *v4) {
  for (int i = 0; i < 16; ++i) {
    dest->arr[i] = v4->arr[i];
  }
}

void mat4_add(mat4 *dest, const mat4 *a, const mat4 *b) {
  for (int i = 0; i < 16; ++i) {
    dest->arr[i] = a->arr[i] + b->arr[i];
  }
}

void mat4_sub(mat4 *dest, const mat4 *a, const mat4 *b) {
  for (int i = 0; i < 16; ++i) {
    dest->arr[i] = a->arr[i] - b->arr[i];
  }
}

void mat4_scale(mat4 *dest, const mat4 *a, float by) {
  for (int i = 0; i < 16; ++i) {
    dest->arr[i] = a->arr[i] * by;
  }
}

void mat4_mul(mat4 *dest, const mat4 *a, const mat4 *b) {
  dest->_11 =
      a->_11 * b->_11 + a->_12 * b->_21 + a->_13 * b->_31 + a->_14 * b->_41;
  dest->_12 =
      a->_11 * b->_12 + a->_12 * b->_22 + a->_13 * b->_32 + a->_14 * b->_42;
  dest->_13 =
      a->_11 * b->_13 + a->_12 * b->_23 + a->_13 * b->_33 + a->_14 * b->_43;
  dest->_14 =
      a->_11 * b->_14 + a->_12 * b->_24 + a->_13 * b->_34 + a->_14 * b->_44;
  dest->_21 =
      a->_21 * b->_11 + a->_22 * b->_21 + a->_23 * b->_31 + a->_24 * b->_41;
  dest->_22 =
      a->_21 * b->_12 + a->_22 * b->_22 + a->_23 * b->_32 + a->_24 * b->_42;
  dest->_23 =
      a->_21 * b->_13 + a->_22 * b->_23 + a->_23 * b->_33 + a->_24 * b->_43;
  dest->_24 =
      a->_21 * b->_14 + a->_22 * b->_24 + a->_23 * b->_34 + a->_24 * b->_44;
  dest->_31 =
      a->_31 * b->_11 + a->_32 * b->_21 + a->_33 * b->_31 + a->_34 * b->_41;
  dest->_32 =
      a->_31 * b->_12 + a->_32 * b->_22 + a->_33 * b->_32 + a->_34 * b->_42;
  dest->_33 =
      a->_31 * b->_13 + a->_32 * b->_23 + a->_33 * b->_33 + a->_34 * b->_43;
  dest->_34 =
      a->_31 * b->_14 + a->_32 * b->_24 + a->_33 * b->_34 + a->_34 * b->_44;
  dest->_41 =
      a->_41 * b->_11 + a->_42 * b->_21 + a->_43 * b->_31 + a->_44 * b->_41;
  dest->_42 =
      a->_41 * b->_12 + a->_42 * b->_22 + a->_43 * b->_32 + a->_44 * b->_42;
  dest->_43 =
      a->_41 * b->_13 + a->_42 * b->_23 + a->_43 * b->_33 + a->_44 * b->_43;
  dest->_44 =
      a->_41 * b->_14 + a->_42 * b->_24 + a->_43 * b->_34 + a->_44 * b->_44;
}

void mat4_transpose(mat4 *dest, const mat4 *m4) {
  /* *
   * | a b c d |    | a e i m |
   * | e f g h | -> | b f j n |
   * | i j k l |    | c g k o |
   * | m n o p |    | d h l p |
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

void mat4_cofactor(mat4 *dest, const mat4 *m4) {
  /* *
   * |  +  -  +  - |
   * |  -  +  -  + |
   * |  +  -  +  - |
   * |  -  +  -  + |
   * */
  dest->_11 = m4->_11;
  dest->_12 = -m4->_12;
  dest->_13 = m4->_13;
  dest->_14 = -m4->_14;
  dest->_21 = -m4->_21;
  dest->_22 = m4->_22;
  dest->_23 = -m4->_23;
  dest->_24 = m4->_24;
  dest->_31 = m4->_31;
  dest->_32 = -m4->_32;
  dest->_33 = m4->_33;
  dest->_34 = -m4->_34;
  dest->_41 = -m4->_41;
  dest->_42 = m4->_42;
  dest->_43 = -m4->_43;
  dest->_44 = m4->_44;
}

void mat4_cut(mat3 *dest, const mat4 *m4, int row, int col) {
  /* *
   * r = 0
   * c = 0
   *     0 1 2 3
   * 0 | a b c d |    | -   -  -  - |
   * 1 | e f g h | -> | - | f  g  h |
   * 2 | i j k l |    | - | j  k  l |
   * 3 | m n o p |    | - | n  o  p |
   * */
  int n = 0;
  for (int i = 0; i < 4; ++i) {
    for (int j = 0; j < 4; ++j) {
      if (i == row or j == col) {
        continue;
      }
      int t = n;
      ++n;
      dest->arr[t] = m4->arr[4 * i + j];
    }
  }
}

void mat4_minor(mat4 *dest, const mat4 *m4) {
  for (int i = 0; i < 4; ++i) {
    for (int j = 0; j < 4; ++j) {
      mat3 m3;
      mat4_cut(&m3, m4, i, j);
      dest->arr[4 * i + j] = mat3_determinant(&m3);
    }
  }
}

void mat4_adjugate(mat4 *dest, const mat4 *m4) {
  mat4 cof;
  mat4_cofactor(&cof, m4);
  mat4_transpose(dest, &cof);
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

void mat4_translation(mat4 *dest, float x, float y, float z) {
  /* *
   * row order
   * |  1  0  0  0 |
   * |  0  1  0  0 |
   * |  0  0  1  0 |
   * ---------------
   * |  x  y  z  1 |
   * */

  mat4 i = mat4_identity();
  i._41 = x;
  i._42 = y;
  i._43 = z;
  mat4_cpy(dest, &i);
}

void mat4_scaling(mat4 *dest, float x, float y, float z) {
  /* *
   * | x  0  0  0 |
   * | 0  y  0  0 |
   * | 0  0  z  0 |
   * | 0  0  0  1 |
   * */

  mat4 i = mat4_identity();
  i._11 = x;
  i._22 = y;
  i._33 = z;
  mat4_cpy(dest, &i);
}

void mat4_rotate_x(mat4 *dest, float by) {
  /* *
   * | 1    0     0  0 |
   * | 0  cos  -sin  0 |
   * | 0  sin   cos  0 |
   * | 0    0     0  1 |
   * */
  by = DEG2RAD(by);

  mat4 i = mat4_identity();
  i._22 = cosf(by);
  i._23 = -sinf(by);
  i._32 = sinf(by);
  i._33 = cosf(by);
  mat4_cpy(dest, &i);
}

void mat4_rotate_y(mat4 *dest, float by) {
  /* *
   * |  cos  0  sin  0 |
   * |    0  1    0  0 |
   * | -sin  0  cos  0 |
   * |    0  0    0  1 |
   * */
  by = DEG2RAD(by);

  mat4 i = mat4_identity();
  i._11 = cosf(by);
  i._13 = sinf(by);
  i._31 = -sinf(by);
  i._33 = cosf(by);
  mat4_cpy(dest, &i);
}

void mat4_rotate_z(mat4 *dest, float by) {
  /* *
   * |  cos  sin  0  0 |
   * | -sin  cos  0  0 |
   * |    0    0  1  0 |
   * |    0    0  0  1 |
   * */
  by = DEG2RAD(by);

  mat4 i = mat4_identity();
  i._11 = cosf(by);
  i._12 = sinf(by);
  i._21 = -sinf(by);
  i._22 = cosf(by);
  mat4_cpy(dest, &i);
}

void mat4_rotation(mat4 *dest, float pitch, float yaw, float roll) {
  mat4 z;
  mat4_rotate_z(&z, roll);
  mat4 x;
  mat4_rotate_x(&x, pitch);
  mat4 y;
  mat4_rotate_y(&y, yaw);

  mat4 xy;
  mat4_mul(&xy, &x, &y);
  mat4_mul(dest, &z, &xy);
}

void mat4_print(const mat4 *m4) {
  printf("%.3f, %.3f, %.3f, %.3f \n", m4->_11, m4->_12, m4->_13, m4->_14);
  printf("%.3f, %.3f, %.3f, %.3f \n", m4->_21, m4->_22, m4->_23, m4->_24);
  printf("%.3f, %.3f, %.3f, %.3f \n", m4->_31, m4->_32, m4->_33, m4->_34);
  printf("%.3f, %.3f, %.3f, %.3f \n", m4->_41, m4->_42, m4->_43, m4->_44);
}
