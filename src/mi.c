/*
Copyright 2010-2024 Gabriele Sales <gabriele.sales@unipd.it>


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

#include "iqsort.h"
#include "mi.h"
#include <math.h>
#include <R.h>


static void init_psi(mi_t* const m) {
  static const coord_t c = 0.5772156649015328606065;

  m->psi = R_Calloc(m->n, coord_t);
  m->psi[0] = -c;

  int i;
  for (i = 1; i < m->n; i++)
    m->psi[i] = m->psi[i-1] + 1.0/i;
}

static coord_t get_psi(const mi_t* const m, const int i) {
  return m->psi[i-1];
}

#define compare_coords(a, b) ((**a) < (**b))

static void sort_coords(const coord_t* const cs, coord_t* const scs, int* const iis, const int n) {
  const coord_t* cps[n];

  int i;
  for (i = 0; i < n; i++)
    cps[i] = (coord_t*)(cs+i);

  QSORT(const coord_t*, cps, n, compare_coords);

  for (i = 0; i < n; i++) {
    scs[i] = *cps[i];
    iis[cps[i]-cs] = i;
  }
}

static dist_t find_range(const coord_t* const cs, const int center_idx, const int* const kis, const int k) {
  int i;
  dist_t md = 0;
  for (i = 0; i < k; i++) {
    const dist_t d = dist_abs(cs[center_idx] - cs[kis[i]]);
    if (d > md) md = d;
  }
  return md;
}

static int region_count(const coord_t* const scs, const int n, const int center_idx, const dist_t range) {
  const coord_t center = scs[center_idx];
  int c = 0;

  int i = center_idx-1;
  while (i >= 0 && center - scs[i] <= range) {
    c++;
    i--;
  }

  i = center_idx+1;
  while (i < n && scs[i] - center <= range) {
    c++;
    i++;
  }

  return c;
}

int make_mi(mi_t* const m, const int n, const int k) {
  if (n < k) return 0;

  m->k = k;
  m->n = n;
  init_psi(m);
  m->sxs  = R_Calloc(n, coord_t);
  m->xiis = R_Calloc(n, int);
  m->sys  = R_Calloc(n, coord_t);
  m->yiis = R_Calloc(n, int);

  return 1;
}

void destroy_mi(mi_t* const m) {
  R_Free(m->sxs);
  R_Free(m->xiis);
  R_Free(m->sys);
  R_Free(m->yiis);
  R_Free(m->psi);
}

coord_t mutual_information(mi_t* const m, const coord_t* const xs, const coord_t* const ys) {
  const coord_t* pxs;
  const coord_t* pys;
  make_grid(&m->grid, xs, ys, m->n, m->k);
  ordered_points(&m->grid, &pxs, &pys);

  sort_coords(pxs, m->sxs, m->xiis, m->n);
  sort_coords(pys, m->sys, m->yiis, m->n);

  int i;
  coord_t accum = 0;
  for (i = 0; i < m->n; i++) {
    int kis[m->k];
    search_knn(&m->grid, pxs[i], pys[i], kis);

    const dist_t mdx = find_range(pxs, i, kis, m->k);
    const int nx = region_count(m->sxs, m->n, m->xiis[i], mdx);

    const dist_t mdy = find_range(pys, i, kis, m->k);
    const int ny = region_count(m->sys, m->n, m->yiis[i], mdy);

    accum += get_psi(m, nx) + get_psi(m, ny);
  }

  destroy_grid(&m->grid);

  return get_psi(m, m->k) + get_psi(m, m->n) - (1.0/m->k) - (accum/m->n);
}
