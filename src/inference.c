/*
Copyright 2011 Gabriele Sales <gabriele.sales@unipd.it>


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

#include <math.h>
#include <string.h>


#define ref(m,x,y) m[n*x+y]
#define dref(m,x,y) ref(m,x,y) = ref(m,y,x)


void aracne(const double* const mis, const int* const np, const double* const epsp, const double* const etap, double* const res) {
  const int n = *np;
  const double eps = *epsp;
  const double eta = *etap;
  int i, j, k;

  memcpy(res, mis, sizeof(double)*n*n);

#ifdef _OPENMP
  #pragma omp parallel for schedule(dynamic) private(j,k)
#endif
  for (i = 2; i < n; i++) {
    for (j = 1; j < i; j++) {
      const double x = ref(mis, i, j);

      for (k = 0; k < j; k++) {
	const double y = ref(mis, j, k);
	const double z = ref(mis, i, k);

	if (x <= y) {
	  if (z < x*eta-eps)
	    dref(res, i, k) = 0;
	  else if (x < y*eta-eps && x < z*eta-eps)
	    dref(res, i, j) = 0;
	} else {
	  if (x <= z) {
	    if (y < x*eta-eps)
	      dref(res, j, k) = 0;
	    else if (x < y*eta-eps && x < z*eta-eps)
	      dref(res, i, j) = 0;
	  } else {
	    if (y < z*eta-eps)
	      dref(res, j, k) = 0;
	    else if (z < x*eta-eps && z < y*eta-eps)
	      dref(res, i, k) = 0;
	  }
	}
      }
    }
  }
}

void clr(const double* const mis, const int* const np, double* const res) {
  const int n = *np;
  double ms[n], vs[n];

#ifdef _OPENMP
  #pragma omp parallel
#endif
  {
    int i, j;

#ifdef _OPENMP
    #pragma omp for private(j)
#endif
    for (i = 0; i < n; i++) {
      double acc;

      for (j = 0, acc = 0; j < n; j++)
	acc += ref(mis, i, j);
      ms[i] = acc / n;

      for (j = 0, acc = 0; j < n; j++) {
	const double v = ref(mis, i, j) - ms[i];
	acc += v*v;
      }
      vs[i] = acc / n;
    }

#ifdef _OPENMP
    #pragma omp for schedule(dynamic) private(j)
#endif
    for (i = 1; i < n; i++) {
      ref(res, i, i) = 0;

      for (j = 0; j < i; j++) {
	const double m = ref(mis, i, j);

	double zi = m - ms[i];
	zi = (zi < 0 || vs[i] == 0) ? 0 : zi*zi/vs[i];

	double zj = m - ms[j];
	zj = (zj < 0 || vs[j] == 0) ? 0 : zj*zj/vs[j];

	dref(res, i, j) = sqrt(zi*zi + zj*zj);
      }
    }
  }
}

void mrnet(const double* const mis, const int* const np, double* const res) {
  const int n = *np;
  int i, j;

  for (i = 0; i < n*n; i++)
    res[i] = 0;

#ifdef _OPENMP
  #pragma omp parallel for private(j)
#endif
  for (i = 0; i < n; i++) {
    int turn;
    char selected[n];
    double rs[n];

    turn = 0;
    memset(selected, 0, sizeof(char)*n);
    selected[i] = 1;
    for (j = 0; j < n; j++)
      rs[j] = 0;

    while (turn < n) {
      double max = 0;
      int max_idx;

      for (j = 0; j < n; j++) {
	if (selected[j])
	  continue;

	double s = ref(mis, i, j);
	if (turn)
	  s -= rs[j] / turn;

	if (s > max) {
	  max = s;
	  max_idx = j;
	}
      }

      if (max <= 0)
	break;

      ref(res, i, max_idx) = max;
      selected[max_idx] = 1;
      for (j = 0; j < n; j++) {
	if (!selected[j])
	  rs[j] += ref(mis, max_idx, j);
      }

      turn++;
    }
  }

  for (i = 1; i < n; i++) {
    for (j = 0; j < i; j++) {
      const double a = ref(res, i, j);
      const double b = ref(res, j, i);
      if (a > b)
	ref(res, j, i) = a;
      else
	ref(res, i, j) = b;
    }
  }
}
