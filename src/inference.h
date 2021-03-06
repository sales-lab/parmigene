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

void aracne(const double* const mis, const int* const np, const double* const epsp, const double* const etap, double* const res);
void clr(const double* const mis, const int* const np, double* const res);
void mrnet(const double* const mis, const int* const np, double* const res);
