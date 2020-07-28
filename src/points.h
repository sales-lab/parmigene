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

#ifndef POINTS_H
#define POINTS_H

#include <float.h>


typedef double coord_t;
typedef double dist_t;

typedef struct {
  coord_t xmin;
  coord_t xmax;
  coord_t ymin;
  coord_t ymax;
} bbox_t;

#define coord_sqrt(v) sqrt(v)
#define coord_ceil(v) ceil(v)
#define coord_read(p, e) strtod(p, e)
#define COORD_MIN DBL_MIN
#define COORD_MAX DBL_MAX
#define DIST_MIN  DBL_MIN
#define DIST_MAX  DBL_MAX
#define dist_abs(v) fabs(v)


int normalize(coord_t* const cs, const int n);
void add_noise(coord_t* const cs, const int n, const double noise, unsigned int* seed);

#endif //POINTS_H
