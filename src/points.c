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

#include "points.h"
#include <math.h>
#include <stdlib.h>

#ifndef isinf
#define isinf(x)						 \
  (sizeof (x) == sizeof (long double) ? isinf_ld (x)		 \
   : sizeof (x) == sizeof (double) ? isinf_d (x)		 \
   : isinf_f (x))
static inline int isinf_f  (float       x)
{ return !isnan (x) && isnan (x - x); }
static inline int isinf_d  (double      x)
{ return !isnan (x) && isnan (x - x); }
static inline int isinf_ld (long double x)
{ return !isnan (x) && isnan (x - x); }
#endif


int normalize(coord_t* const cs, const int n) {
  int i;
  coord_t m, s, l;

  for (i = 0, m = 0; i < n; i++)
    m += cs[i];
  m /= n;
  if (isinf(m)) return 0;

  for (i = 0, s = 0; i < n; i++)
    s += (cs[i]-m)*(cs[i]-m);
  s = coord_sqrt(s/(n-1));
  if (isinf(s)) return 0;

  if (s > 0) {
    for (i = 0, l = COORD_MAX; i < n; i++) {
      cs[i] = (cs[i]-m)/s;
      if (cs[i] < l) l = cs[i];
    }

    for (i = 0; i < n; i++)
      cs[i] -= l;
  }

  return 1;
}

#ifdef WIN32
// This version is in POSIX/C99 for [s]rand()
#define rand_r myrand_r
#undef RAND_MAX
#define RAND_MAX 31767
static int rand_r (unsigned int *seed)
{
    *seed = *seed * 1103515245 + 12345;
    return((unsigned)(*seed/65536) % 32768);
}
#endif

void add_noise(coord_t* const cs, const int n, const double noise, unsigned int* seed) {
  int i;
  for (i = 0; i < n; i++)
    cs[i] += (1.0*rand_r(seed)/RAND_MAX) * noise;
}
