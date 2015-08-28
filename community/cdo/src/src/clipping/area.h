/**
 * @file area.h
 *
 * @copyright Copyright  (C)  2013 Rene Redler <rene.redler@mpimet.mpg.de>
 *
 * @version 1.0
 * @author Rene Redler <rene.redler@mpimet.mpg.de>
 */
/*
 * Keywords:
 * Maintainer: Rene Redler <rene.redler@mpimet.mpg.de>
 * URL: https://redmine.dkrz.de/doc/YAC/html/index.html
 *
 * This file is part of YAC.
 *
 * YAC is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * YAC is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with YAC.  If not, see <http://www.gnu.org/licenses/gpl.txt>.
 */

#ifndef AREA_H
#define AREA_H

#include "grid.h"
#include "clipping.h"

/** \example test_area.c
 * A test for different area calculations.
 */

/** \file area.h
  * \brief Structs and interfaces for area calculations
  *
  **/

double yac_area_tol();

/** \brief Calculate the area of a cell on a unit sphere
  *
  * Taken from:
  *
  * http://www.psjournal.com/gis/Area-of-Polygon-on-a-Sphere-Algorithm-2534-.htm
  *
  * Bob (Robert G.) Chamberlain <rgc@jpl.nasa.gov>
  *
  * See also:
  * Robert.G. Chamberlain, William.H. Duquette
  * Jet Propulsion Laboratory
  * Some Algorithms for Polygons on a Sphere.
  * Association of American Geographers Annual Meeting San Francisco, California
  * 17 - 21 April 2007
  *
  * Given a cell with N edges, each of which is a segment of a great
  * circle. The cell must be simply connected (a single piece), without
  * holes, and no edge may cross another. Neither pole may be inside the
  * cell.
  *
  * The cell is described in a counterclockwise direction by a succession
  * of vertices, numbered from 0 to N-1. A vertex N, if given, is identical
  * to vertex 0. Note that "inside" and "outside" are defined by the
  * requirement that the cell be counterclockwise. (If the cell is
  * described in a clockwise direction, the area computed will be negative.)
  *
  * The location of each vertex is given by latitude and longitude. In this
  * algorithm, latitude and longitude are expressed in radian.
  *
  * Latitude is zero at the equator; north is positive, south is negative.
  * The latitude of point i is denoted by lat.i or lat.(i)
  *
  * Longitude is zero at Greenwich: east is positive, west is negative. The
  * longitude of point i is denoted by lon.i or lon.(i)
  *
  * With the radius of the Earth denoted by R, the area of the cell,
  * expressed in the square of the units used for R, is given by:
  *
  * Area = - ( R^2 / 2 ) * ( SUM )
  *
  * SUM = ( lon.1     - lon.(N-1) ) * sin( lat.0 )
  *     + ( lon.2     - lon.0     ) * sin( lat.1 )
  *     + ( lon.3     - lon.1     ) * sin( lat.2 )
  *     + ( lon.4     - lon.2     ) * sin( lat.3 )
  *     + ...
  *     + ( lon.(N-1) - lon.(N-3) ) * sin( lat.(N-2) )
  *     + ( lon.0     - lon.(N-2) ) * sin( lat.(N-1) ) 
  *
  * \remark This does not work for cells that include the pole.
  *
 **/

double yac_cell_approx_area ( struct grid_cell cell );

/** \brief Calculate the area of a triangle on a unit sphere
  *
  * Adopted from the ICON code, mo_base_geometry.f90
  * provided by Luis Kornblueh, MPI-M. 
  *
  * The algorithm is based on Girards' theorem and formula.
  *
  * Converted to c by Rene Redler, MPI-M.
  *
  * Vertex coordinates need to be provided as cartesian coordinates
  *
  *  The algorithm is based on Girards' theorem and formula.
  *
  *  R:  Earth radius
  *  n:  number of vertices
  *  pi: guess what
  *  Theta: Sum of inner angle of the element (in rad)
  *
  *  The Formula reads as
  *
  *  S = [ Theta - (n-2) * pi ] * R*R
  *
  *  Ad with n=3 for triangles simplifies to
  *
  *  S = [ Theta - pi ] * R*R
  *
  *  @param[in] cell cell for which he area shal be calculated
  *
  *  @return approximate area of the cell
  *
  **/

double yac_triangle_area ( struct grid_cell cell );

/** \brief Calculate the area of a cell on a unit sphere
  *
  * Generalized version of triangle_area
  *
  * The algorithm is based on Girards' theorem and formula.
  *
  * Converted to c by Rene Redler, MPI-M.
  *
  * Vertex coordinates need to be provided as cartesian coordinates
  *
  *  The algorithm is based on Girards' theorem and formula.
  *
  *  R:  Earth radius
  *  n:  number of vertices
  *  pi: guess what
  *  Theta: Sum of inner angle of the element (in rad)
  *
  *  The Formula reads as
  *
  *  S = [ Theta - (n-2) * pi ] * R*R
  *
  *  @param[in] cell cell for which he area shall be calculated
  *
  *  @return area of the triangle
  *
  **/

double yac_cell_area ( struct grid_cell cell );

/** \brief Calculate the area of a cell on a unit sphere
  *
  * Following M. Bevis and G. Cambareri, 1987: Computing the Area
  *   of a Spherical Polygon of Arbitray Shape. Math. Geology,
  *   Vol. 19, No 4, 335 - 346.
  *
  * The algorithm is based on Girards' theorem and formula.
  *
  * Converted to c and expanded for cells by Rene Redler, MPI-M.
  *
  * In contrast to the above simple algorithm by Chamberlain and Duquette
  * this approach seems to work over the pole as well.
  *
  *  R:  Earth radius
  *  n:  number of vertices
  *  pi: guess what
  *  Theta: Sum of inner angle of the element (in rad)
  *
  *  The Formula reads as
  *
  *  S = [ Theta - (n-2) * pi ] * R*R
  *
  *  This algorithm does not work for very small elements. In this case negative
  *  areas are delivered which are set to zero currently. An alternative would be
  *  to transform the coordinates into the s space, ignore the curvature
  *  of the edges (as they can be neglected for small elements) and calculate the
  *  area of a plane triangle.
  *
  *  http://www.mathopenref.com/coordpolygonarea2.html \n
  *  http://geomalgorithms.com/a01-_area.html \n
  *  http://stackoverflow.com/questions/2350604/get-the-area-of-a-3d-surface
  *  
  *  @param[in] cell cell for which the area shall be calculated
  *
  *  @return area of the cell
  *
  **/
double yac_girards_area ( struct grid_cell cell );

/** \brief Calculate the area of a cell in a 3d plane on a unit sphere
  *
  *  see http://geomalgorithms.com/a01-_area.html
  *
  *  other references:
  *
  *  http://www.mathopenref.com/coordpolygonarea2.html \n
  *  http://geomalgorithms.com/a01-_area.html \n
  *  http://stackoverflow.com/questions/2350604/get-the-area-of-a-3d-surface
  *
  *  @param[in] cell cell for which the area shall be calculated
  *
  *  @return area of the cell
  *
  **/

double yac_pole_area ( struct grid_cell cell );

/**
  * \brief Area calculation on a unit sphere of a planar polygon in 3D
  *
  * (http://gaim.umbc.edu/2010/06/03/polygon-area)\n
  *
  * This area calculation works for any planar polygon (concave or convex)
  * with non-intersecting edges in 3D. It requires vertex coordinates in
  * Carthesian space. In our case this is applicable for very
  * small elements on the sphere.
  *
  */

double yac_planar_3dcell_area (struct grid_cell cell);

/**
  * \brief Area calculation on a unit sphere taken from ESMF based on L'Huilier's Theorem
  *
  * (http://mathworld.wolfram.com/LHuiliersTheorem.html)\n
  * (http://mathforum.org/library/drmath/view/65316.html)\n
  * (http://math.stackexchange.com/questions/9819/area-of-a-spherical-triangle)
  *
  * The cell in split up into triangles that all have one corner in common,
  * then the area for each of the triangle is computed and summed up to build
  * the area of the cell. L'Huilier's Theorem is used to compute the area of
  * the triangles. This seems to be sufficiently accurate for elements on the
  * Earth surface with edge lengths of approx. 100 m and gives results comparable
  * to our implementation of Huilier's algorithm for edge lengths of up to 1 km.
  */
double yac_huiliers_area(struct grid_cell cell);

#endif // AREA_H
