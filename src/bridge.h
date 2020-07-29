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

void mi_single(double* const x, double* const y, const int* const np, const int* const kp, const double* const noisep, double* const res);
void mi_cross(double* const xs, const int* const lp, double* const ys, const int* const mp, const int* const np, const int* const kp, const double* const noisep, double* res);
void mi_all(double* const xs, const int* const lp, const int* const np, const int* const kp, const double* const noisep, double* res);
