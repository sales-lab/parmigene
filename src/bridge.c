/*
Copyright 2010-2011 Gabriele Sales <gabriele.sales@unipd.it>


This file is part of parmigene.

knnmi is free software: you can redistribute it and/or modify
it under the terms of the GNU Affero General Public License 
version 3 as published by the Free Software Foundation.

knnmi is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU Affero General Public License for more details.

You should have received a copy of the GNU Affero General Public
License along with parmigene. If not, see <http://www.gnu.org/licenses/>.
*/

#include "mi.h"
#include "points.h"
#include <math.h>


unsigned int gen_seed(const double* const cs, const int n, const int k) {
  return n * k * ((int)cs[n/2]*100);
}


void mi_single(double* const x, double* const y, const int* const np, const int* const kp, const double* const noisep, double* const res) {
  const int n = *np;
  const int k = *kp;
  const double noise = *noisep;
  unsigned int seed = gen_seed(x, n, k);

  if (normalize(x, n) && normalize(y, n)) {
    add_noise(x, n, noise, &seed);
    add_noise(y, n, noise, &seed);

    mi_t mi;
    make_mi(&mi, n, k);
    *res = mutual_information(&mi, x, y);
    destroy_mi(&mi);
  }
  else
    *res = NAN;
}

void mi_cross(double* const xs, const int* const lp, double* const ys, const int* const mp, const int* const np, const int* const kp, const double* const noisep, double* res) {
  const int l = *lp;
  const int m = *mp;
  const int n = *np;
  const int k = *kp;
  const double noise = *noisep;
  int xnormed[l], ynormed[m];

#ifdef _OPENMP
  #pragma omp parallel
#endif
  {
    int i, j;
    unsigned int seed = gen_seed(xs, l*n, k);

#ifdef _OPENMP
    #pragma omp for nowait
#endif
    for (i = 0; i < l; i++) {
      double* const p = xs+(i*n);
      xnormed[i] = normalize(p, n);
      add_noise(p, n, noise, &seed);
    }

#ifdef _OPENMP
    #pragma omp for
#endif
    for (j = 0; j < m; j++) {
      double* const p = ys+(j*n);
      ynormed[j] = normalize(p, n);
      add_noise(p, n, noise, &seed);
    }

    mi_t mi;
    make_mi(&mi, n, k);

#ifdef _OPENMP
    #pragma omp for schedule(dynamic)
#endif
    for (i = 0; i < l; i++)
      for (j = 0; j < m; j++)
	res[i*m+j] = (xnormed[i] && ynormed[j]) ? mutual_information(&mi, xs+(i*n), ys+(j*n)) : NAN;

    destroy_mi(&mi);
  }
}

void mi_all(double* const xs, const int* const lp, const int* const np, const int* const kp, const double* const noisep, double* res) {
  const int l = *lp;
  const int n = *np;
  const int k = *kp;
  const double noise = *noisep;
  int xnormed[l];

#ifdef _OPENMP
  #pragma omp parallel
#endif
  {
    int i, j;
    unsigned int seed = gen_seed(xs, l*n, k);

#ifdef _OPENMP
    #pragma omp for nowait
#endif
    for (i = 0; i < l; i++) {
      double* const p = xs+(i*n);
      xnormed[i] = normalize(p, n);
      add_noise(p, n, noise, &seed);
    }

    #pragma omp for
    for (i = 0; i < l; i++)
      res[i*l+i] = 0.0;

    mi_t mi;
    make_mi(&mi, n, k);

#ifdef _OPENMP
    #pragma omp for schedule(dynamic)
#endif
    for (i = 1; i < l; i++)
      for (j = 0; j < i; j++)
    	res[i*l+j] = res[j*l+i] = (xnormed[i] && xnormed[j]) ? mutual_information(&mi, xs+(i*n), xs+(j*n)) : NAN;

    destroy_mi(&mi);
  }
}
