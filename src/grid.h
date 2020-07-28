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

#ifndef GRID_H
#define GRID_H

#include "points.h"

typedef struct {
  int fill;
  int size;
  int* idxs;
  coord_t *xs, *ys;
} cell_t;

typedef struct candidate {
  int idx;
  dist_t dist;
  struct candidate* next;
} candidate_t;

typedef struct {
  int size;
  int used;
  candidate_t* candidates;
  candidate_t front;
  dist_t max_dist;
} candidates_t;

typedef struct {
  coord_t* xs;
  coord_t* ys;
  int k;
  coord_t xmin, ymin;
  dist_t size;
  int cols, lines;
  cell_t* cells;
  candidates_t candidates;
} grid_t;

void make_grid(grid_t* const g, const coord_t* const xs, const coord_t* const ys, const int n, const int k);
void destroy_grid(grid_t* const g);
void ordered_points(const grid_t* const g, const coord_t** xs, const coord_t** ys);
void search_knn(grid_t* const g, const coord_t x, const coord_t y, int* const ris);

#endif //GRID_H
