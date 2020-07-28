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

#include <float.h>
#include "grid.h"
#include <math.h>
#include <R.h>


/*
 * The following two functions (distance and bounding_box) really belong to points.c.
 * Keeping them in a separate translation unit, however, means that they won't be inlined
 * unless link time optimizations are enabled.
 * As gcc only recently started to support LTO (and such extension is still considered
 * somewhat experimental), we manually move these functions here to be sure they are
 * considered for inlining.
 */
static dist_t distance(const coord_t x1, const coord_t y1, const coord_t x2, const coord_t y2) {
  const dist_t dx = dist_abs(x2 - x1);
  const dist_t dy = dist_abs(y2 - y1);
  return dx > dy ? dx : dy;
}

static void bounding_box(const coord_t* const xs, const coord_t* const ys, const int n, bbox_t* const bbox) {
  bbox->xmin = bbox->ymin = COORD_MAX;
  bbox->xmax = bbox->ymax = COORD_MIN;

  int i;
  for (i = 0; i < n; i++) {
    if (xs[i] < bbox->xmin) bbox->xmin = xs[i];
    if (xs[i] > bbox->xmax) bbox->xmax = xs[i];
    if (ys[i] < bbox->ymin) bbox->ymin = ys[i];
    if (ys[i] > bbox->ymax) bbox->ymax = ys[i];
  }
}


static int min(const int a, const int b) {
  return a < b ? a : b;
}

static int max(const int a, const int b) {
  return a > b ? a : b;
}

static dist_t dist_max(const dist_t a, const dist_t b) {
  return a > b ? a : b;
}

static int coords_to_idx(const grid_t* const g, const coord_t x, const coord_t y) {
  int i = (x - g->xmin) / g->size;
  if (i == g->cols) i--;

  int j = (y - g->ymin) / g->size;
  if (j == g->lines) j--;

  return j*g->cols + i;
}

static dist_t wall_distance(const grid_t* const g, const coord_t x, const coord_t y, const int ci, const int cj) {
  dist_t dmin, d;

  const coord_t wl = g->xmin + ci*g->size;
  dmin = x - wl;

  d = wl+g->size - x;
  if (d < dmin) dmin = d;

  const coord_t wb = g->ymin + cj*g->size;
  d = y - wb;
  if (d < dmin) dmin = d;

  d  = wb+g->size - y;
  if (d < dmin) dmin = d;

  return dmin;
}


static void append_index(cell_t* const cell, const int idx, const int load_factor) {
  if (cell->fill == cell->size) {
    cell->size = cell->size == 0 ? load_factor : cell->size*2;
    cell->idxs = Realloc(cell->idxs, cell->size, int);
  }
  cell->idxs[cell->fill++] = idx;
}

static void fill_cells(grid_t* const g, const coord_t* const xs, const coord_t* const ys, const int n) {
  const int load_factor = coord_ceil(coord_sqrt(n));

  int i;
  for (i = 0; i < n; i++) {
    const int cell_idx = coords_to_idx(g, xs[i], ys[i]);
    append_index(&g->cells[cell_idx], i, load_factor);
  }

  int c;
  coord_t *px, *py;
  px = g->xs = Calloc(n, coord_t);
  py = g->ys = Calloc(n, coord_t);
  for (c = 0; c < g->lines*g->cols; c++) {
    cell_t* const cell = &g->cells[c];
    cell->xs = px;
    cell->ys = py;
    for (i = 0; i < cell->fill; i++) {
      const int idx = cell->idxs[i];
      *px++ = xs[idx];
      *py++ = ys[idx];
    }
  }

  for (c = 0; c < g->lines*g->cols; c++) {
    if (g->cells[c].idxs)
      Free(g->cells[c].idxs);
  }
}


static void reset_candidates(candidates_t* const cs) {
  cs->used       = 0;
  cs->front.next = NULL;
  cs->max_dist   = DIST_MAX;
}

static void init_candidates(candidates_t* const cs, const int k) {
  cs->size = k+1;
  cs->candidates = (candidate_t*)Calloc(cs->size, candidate_t);
  reset_candidates(cs);
}

static void destroy_candidates(candidates_t* const cs) {
  Free(cs->candidates);
}

static void insert_candidate(const int idx, const dist_t dist, candidates_t* const cs) {
  if (dist >= cs->max_dist)
    return;

  candidate_t* insert_point = &cs->front;
  while (insert_point->next && dist < insert_point->next->dist)
    insert_point = insert_point->next;

  candidate_t* c;
  if (cs->used < cs->size) {
    c = &cs->candidates[cs->used++];

    c->idx = idx;
    c->dist = dist;
    c->next = insert_point->next;
    insert_point->next = c;

    if (cs->used == cs->size)
      cs->max_dist = cs->front.next->dist;
  } else {
    c = cs->front.next;
    c->idx = idx;
    c->dist = dist;

    if (c != insert_point) {
      cs->front.next = c->next;
      c->next = insert_point->next;
      insert_point->next = c;
    }

    cs->max_dist = cs->front.next->dist;
  }
}

static int candidate_last_dist(candidates_t* cs, dist_t* const d) {
  if (cs->used != cs->size)
    return 0;

  *d = cs->front.next->dist;
  return 1;
}

static void extract_candidate_idxs(candidates_t* const cs, const int k, int* const is) {
  const candidate_t* c = cs->front.next;
  int i;
  for (i = 0; i < k; i++) {
    is[k-i-1] = c->idx;
    c = c->next;
  }
}


void make_grid(grid_t* const g, const coord_t* const xs, const coord_t* const ys, const int n, const int k) {
  bbox_t bbox;
  bounding_box(xs, ys, n, &bbox);

  g->k = k;
  g->xmin = bbox.xmin;
  g->ymin = bbox.ymin;

  const dist_t w = bbox.xmax - bbox.xmin;
  const dist_t h = bbox.ymax - bbox.ymin;

  coord_t alpha = 1.23;
  while (1) {
    g->size  = alpha * dist_max(w/coord_sqrt(n), h/coord_sqrt(n));
    g->cols  = max(coord_ceil(w/g->size), 1);
    g->lines = max(coord_ceil(h/g->size), 1);

    if (w / g->size < g->cols && h / g->size < g->lines)
      break;

    alpha += 0.01;
  }

  g->cells = Calloc(g->cols*g->lines, cell_t);
  fill_cells(g, xs, ys, n);
  init_candidates(&g->candidates, k);
}

void destroy_grid(grid_t* const g) {
  destroy_candidates(&g->candidates);
  Free(g->cells);
  Free(g->xs);
  Free(g->ys);
}

void ordered_points(const grid_t* const g, const coord_t** xs, const coord_t** ys) {
  *xs = g->xs;
  *ys = g->ys;
}

void search_knn(grid_t* const g, const coord_t x, const coord_t y, int* const ris) {
  const int ci = (x - g->xmin) / g->size;
  const int cj = (y - g->ymin) / g->size;

  dist_t dsh = wall_distance(g, x, y, ci, cj);
  const int lmax = max(max(ci, g->cols-1-ci), max(cj, g->lines-1-cj));

  reset_candidates(&g->candidates);

  int l, j, c;
  for (l = 0; l <= lmax; l++, dsh += g->size) {
    const int il = ci-l;
    const int ih = ci+l;
    const int jl = cj-l;
    const int jh = cj+l;

    for (j = max(0, jl); j <= min(jh, g->lines-1); j++) {
      int i, istep;
      if (j == jl || j == jh) {
	istep = 1;
	i = max(0, il);
      } else {
	istep = ih - il;
	i = il >= 0 ? il : il+istep;
      }

      for (; i <= min(ih, g->cols-1); i += istep) {
	const cell_t* const cell = &g->cells[j*g->cols + i];
	const int base_idx = cell->xs - g->xs;
	for (c = 0; c < cell->fill; c++) {
	  const dist_t d = distance(x, y, cell->xs[c], cell->ys[c]);
	  insert_candidate(base_idx+c, d, &g->candidates);
	}
      }
    }

    dist_t dmin;
    if (candidate_last_dist(&g->candidates, &dmin) && dmin <= dsh)
      break;
  }

  extract_candidate_idxs(&g->candidates, g->k, ris);
}
