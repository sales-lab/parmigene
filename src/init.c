/*
Copyright 2020 Gabriele Sales <gabriele.sales@unipd.it>


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

#include "bridge.h"
#include "inference.h"
#include <stdlib.h>
#include <R_ext/Rdynload.h>

static const R_CMethodDef entries[] = {
  {"aracne",    (DL_FUNC) &aracne,    5},
  {"clr",       (DL_FUNC) &clr,       3},
  {"mi_all",    (DL_FUNC) &mi_all,    6},
  {"mi_cross",  (DL_FUNC) &mi_cross,  8},
  {"mi_single", (DL_FUNC) &mi_single, 6},
  {"mrnet",     (DL_FUNC) &mrnet,     3},
  {NULL, NULL, 0}
};

void R_init_parmigene(DllInfo *dll) {
  R_registerRoutines(dll, entries, NULL, NULL, NULL);
  R_useDynamicSymbols(dll, FALSE);
  R_forceSymbols(dll, TRUE);
}
