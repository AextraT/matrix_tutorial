#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>

#include <assert.h>
#include "matrix.h"

#define PRINT_PRECISION "6"
#define PRINT_DECIMAL_PRECISION "2"

matrix matrix_create(unsigned n1, unsigned n2, scalar v)
{
  matrix m={n1,n2,true,NULL};
  if(n1==0 || n2==0)
    return m;
  else if(SIZE_MAX / n1 <= n2)
    m.ok = false;
  else if (!(m.data=calloc(((size_t) n1)*n2, sizeof(scalar))))
    m.ok = false;
  else {
    for(unsigned i=0; i<n1; ++i)
      for(unsigned j=0; j<n2; ++j)
        *matrix_get(m,i,j) = v;
  }
 
  return m;
}

matrix matrix_identity(unsigned n)
{
  matrix m = matrix_create(n, n, 0.);
  if(!m.ok)
    return m;

  for(unsigned i=0; i<n; ++i)
    *matrix_get(m, i, i) = 1.;
  return m;
}

void matrix_destroy(matrix m)
{
  if(m.ok) {
    m.ok = false;
    free(m.data);
  }
}

scalar *matrix_get(matrix m, unsigned i, unsigned j)
{
  if(!m.ok || i>m.n1 || i>m.n2)
    return NULL;

  return &m.data[i*m.n2+j];
}

matrix matrix_add(matrix m, matrix n)
{
  matrix res={0,0,false,NULL};

  if(m.n1!=n.n1 || m.n2!=n.n2 || !m.ok || !n.ok)
    return res;

  res=matrix_create(m.n1, m.n2, 0.);
  for(unsigned i=0; i<m.n1; ++i)
    for(unsigned j=0; j<m.n2; ++j)
      *matrix_get(res, i, j) = *matrix_get(m, i, j) + *matrix_get(n, i, j);

  return res;
}

matrix mult(matrix A, matrix B) {
    assert(A.n2 == B.n1);

    matrix M = matrix_create(A.n1, B.n2, 0);
    
    for (unsigned int i = 0; i < M.n1; i++)
        for (unsigned int j = 0; j < M.n2; j++)
            for (unsigned int k = 0; k < A.n2; k++)
                *matrix_get(M, i, j) += *matrix_get(A, i, k) * *matrix_get(B, k, j);

    return M;
}

matrix expo(matrix M, int pow) {
    if (pow == 0)
        return matrix_identity(M.n1);

    matrix ev = expo(M, pow / 2);

    if (pow % 2 == 0) {
        matrix ans = mult(ev, ev);
        matrix_destroy(ev);
        return ans;
    }

    matrix sq = mult(ev, ev);
    matrix_destroy(ev);
    matrix ans = mult(M, sq);
    matrix_destroy(sq);
    return ans;
}

void matrix_print(FILE *f, matrix m)
{
  if(!m.ok)
    fprintf(f, "Invalid matrix\n");
  else {
    for(unsigned i=0; i<m.n1; ++i) {
      for(unsigned j=0; j<m.n2; ++j)
        fprintf(
            f, 
            "%"PRINT_PRECISION"."PRINT_DECIMAL_PRECISION"f ",
            *matrix_get(m, i, j));
      fprintf(f, "\n");
    }
  }
}
