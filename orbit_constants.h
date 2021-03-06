/*  Copyright 2019, 2020 Garrett Wright, Princeton Plasma Physic Lab

    This file is part of CuOrbit.

    CuOrbit is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    CuOrbit is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with CuOrbit.  If not, see <https://www.gnu.org/licenses/>.
*/

#ifndef SET_ORBIT_CONSTANTS_H_
#define SET_ORBIT_CONSTANTS_H_

#include <math.h>

/* note M_PI was dropped for c99 */
#ifndef M_PI
    #define M_PI 3.14159265358979323846
#endif
static const double pi2 = 2. * M_PI;
static const double pi2i = 0.5 / M_PI;

#endif
