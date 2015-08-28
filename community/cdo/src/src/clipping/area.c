/**
 * @file area.c
 *
 * @copyright Copyright  (C)  2013 Moritz Hanke <hanke@dkrz.de>
 *                                 Rene Redler <rene.redler@mpimet.mpg.de>
 *
 * @version 1.0
 * @author Moritz Hanke <hanke@dkrz.de>
 *         Rene Redler <rene.redler@mpimet.mpg.de>
 */
/*
 * Keywords:
 * Maintainer: Moritz Hanke <hanke@dkrz.de>
 *             Rene Redler <rene.redler@mpimet.mpg.de>
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

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <float.h>

#include "area.h"
#include "clipping.h"
#include "geometry.h"
#include "utils.h"
#include "ensure_array_size.h"

// area tolerance (10 m^2 = 0.0001 km^2)
// const double area_tol = 0.02*0.02; // 20m*20m

// an area of 20m x 20m on the Earth Surface is equivalent to an area on the unit sphere:

double yac_area_tol () {
  return 0.02 / EarthRadius;
}

static inline double scalar_product(double a[], double b[]);

static inline double inner_angle ( double plat, double plon,
                                   double qlon, double qlat );

static double yac_partial_area ( double a_lon, double a_lat,
                             double b_lon, double b_lat,
                             double c_lon, double c_lat );

/* ----------------------------------- */

double yac_cell_approx_area ( struct grid_cell cell ) {

  /* adopted from Robert.G. Chamberlain and William.H. Duquette */

  int m, M;
  // double SUM = 0.0;
  double area = 0.0;

  M = cell.num_corners;

  for ( m = 0; m < M; m++ )
    area += ( cell.coordinates_x[m] - cell.coordinates_x[(m+2)%M] ) *
              sin(cell.coordinates_y[(m+1)%M]);

  // area = EarthRadius2 * SUM;

  return area;
}

/* ----------------------------------- */

double yac_triangle_area ( struct grid_cell cell ) {

  /* taken from the ICON code, mo_base_geometry.f90
     provided by Luis Kornblueh, MPI-M. */

  double s01, s12, s20;
  double ca1, ca2, ca3;
  double a1, a2, a3;

  double * triangle[3];
  double u01[3], u12[3], u20[3];

  if ( cell.num_corners != 3 ) {
    printf ("Only for triangles!\n");
    return -1;
  }

  triangle[0] = cell.coordinates_xyz + 0*3;
  triangle[1] = cell.coordinates_xyz + 1*3;
  triangle[2] = cell.coordinates_xyz + 2*3;

  /* First, compute cross products Uij = Vi x Vj. */

  crossproduct_ld(triangle[0], triangle[1], u01);
  crossproduct_ld(triangle[1], triangle[2], u12);
  crossproduct_ld(triangle[2], triangle[0], u20);

  /*  Normalize Uij to unit vectors. */

  s01 = scalar_product(u01, u01);
  s12 = scalar_product(u12, u12);
  s20 = scalar_product(u20, u20);

  /* Test for a degenerated triangle associated with collinear vertices. */

  if ( s01 == 0.0 &&
       s12 == 0.0 &&
       s20 == 0.0 )
    return 0.0;

  s01 = sqrt(s01);
  s12 = sqrt(s12);
  s20 = sqrt(s20);

  for (int m = 0; m < 3; m++ ) {
    u01[m] = u01[m]/s01;
    u12[m] = u12[m]/s12;
    u20[m] = u20[m]/s20;
  }

  /*  Compute interior angles Ai as the dihedral angles between planes:

      CA1 = cos(A1) = -<U01,U20>
      CA2 = cos(A2) = -<U12,U01>
      CA3 = cos(A3) = -<U20,U12>

  */

  ca1 = -u01[0]*u20[0]-u01[1]*u20[1]-u01[2]*u20[2];
  ca2 = -u12[0]*u01[0]-u12[1]*u01[1]-u12[2]*u01[2];
  ca3 = -u20[0]*u12[0]-u20[1]*u12[1]-u20[2]*u12[2];

  if ( ca1 < -1.0 ) ca1 = -1.0;
  if ( ca1 >  1.0 ) ca1 =  1.0;
  if ( ca2 < -1.0 ) ca2 = -1.0;
  if ( ca2 >  1.0 ) ca2 =  1.0;
  if ( ca3 < -1.0 ) ca3 = -1.0;
  if ( ca3 >  1.0 ) ca3 =  1.0;

  a1 = acos(ca1);
  a2 = acos(ca2);
  a3 = acos(ca3);

  /*  Compute areas = a1 + a2 + a3 - pi.

      here for a unit sphere: */

  // return MAX(( a1+a2+a3-M_PI ) * EarthRadius * EarthRadius, 0.0);
  return MAX( (a1+a2+a3-M_PI) , 0.0 );
}

/* ----------------------------------- */

double yac_cell_area ( struct grid_cell cell ) {

  /* generalised version based on the ICON code, mo_base_geometry.f90
     provided by Luis Kornblueh, MPI-M. */

  int const M = cell.num_corners; // number of vertices

  double area;
  double s[cell.num_corners];
  double ca[cell.num_corners];
  double a[cell.num_corners];

  double * p[cell.num_corners];
  double u[cell.num_corners][3];

  for (int i = 0; i < cell.num_corners; ++i)
    p[i] = cell.coordinates_xyz + i * 3;

  /* First, compute cross products Uij = Vi x Vj. */

  for (int m = 0; m < M; m++ )
    crossproduct_ld (p[m], p[(m+1)%M], u[m]);

  /*  Normalize Uij to unit vectors. */

  area = 0.0;

  for (int m = 0; m < M; m++ ) {
    s[m] = scalar_product(u[m], u[m]);
    area += s[m];
  }

  /* Test for a degenerated cells associated with collinear vertices. */

  if ( area != 0.0 ) {

    for (int m = 0; m < M; m++ )
      s[m] = sqrt(s[m]);

    for (int m = 0; m < M; m++ )
      for (int i = 0; i < 3; i++ )
        u[m][i] = u[m][i]/s[m];

    /*  Compute interior angles Ai as the dihedral angles between planes
        by using the definition of the scalar product

                    ab = |a| |b| cos (phi)

        As a and b are already normalised this reduces to

                    ab = cos (phi)

        There is no explanation so far for the - in the loop below.
        But otherwise we don't get the correct results for triangles
        and cells. Must have something to do with the theorem.

     */

    for (int m = 0; m < M; m++ ) {
      ca[m] = - scalar_product(u[m], u[(m+1)%M]);
      if ( ca[m] < -1.0 ) ca[m] = -1.0;
      if ( ca[m] >  1.0 ) ca[m] =  1.0;
      a[m] = acos(ca[m]);
    }

    /*  Compute areas = a1 + a2 + a3 - (M-2) * pi.

        here for a unit sphere: */

    area = - (double) (M-2) * M_PI;

    for (int m = 0; m < M; m++ )
      area += a[m];
  }

  // return MAX(area * EarthRadius * EarthRadius, 0.0);
  return MAX(area,0.0);
}

/* ----------------------------------- */

double yac_girards_area ( struct grid_cell cell  ) {

  /* Bevis and Cambareri, 1987

     this algorithm provides wrong results if one
     of the inner angles becomes less than 1e-18 rad

     (R. Redler, M. Hanke 2013)
  */

  int m;
  const double tol = 1e-18;
  double area = 0.0;

  int M = cell.num_corners;
  if (M < 3) return area;  // a degenerate cell

  double * theta = malloc ( M * sizeof(theta[0]) );

  for ( m = 0; m < M; m++ ) {
     theta[m] = yac_partial_area(cell.coordinates_x[(m+1)%M], cell.coordinates_y[(m+1)%M],
                             cell.coordinates_x[(m+2)%M], cell.coordinates_y[(m+2)%M],
                             cell.coordinates_x[m%M], cell.coordinates_y[m%M]);
     if ( theta[m] < tol ) {
       area = yac_planar_3dcell_area ( cell );
       return area;
     }
  }

  /*  Sum up partial areas */

  area = - (double) (M-2) * M_PI;

  for ( m = 0; m < M; m++ )
    area += theta[m];

  /* Area on Sphere with radius EarthRadius */

  // area *= EarthRadius * EarthRadius;

  free ( theta );

  if ( area <= 0.0 ) area = yac_planar_3dcell_area ( cell );

  return area;
}

/* ----------------------------------- */

/** area of a spherical triangle based on L'Huilier's Theorem
  *
  * source code is taken from code by Robert Oehmke of Earth System Modeling
  * Framework (www.earthsystemmodeling.org)
  *
  * it has been extended by a more accurate computation of vector angles
  *
  * the license statement for this routine is as follows:
  * Earth System Modeling Framework
  * Copyright 2002-2013, University Corporation for Atmospheric Research,
  * Massachusetts Institute of Technology, Geophysical Fluid Dynamics
  * Laboratory, University of Michigan, National Centers for Environmental
  * Prediction, Los Alamos National Laboratory, Argonne National Laboratory,
  * NASA Goddard Space Flight Center.
  * Licensed under the University of Illinois-NCSA License.
  *
  * \remark all edges are on great circle
  */
static double
tri_area(double u[3], double v[3], double w[3]) {

  double a_ = get_vector_angle(u,v);
  double b_ = get_vector_angle(u,w);
  double c_ = get_vector_angle(w,v);

  double a, b, c;

  if (a_ < b_) {
    if (a_ < c_) {
      if (b_ < c_) {
        a = a_, b = b_, c = c_;
      } else {
        a = a_, b = c_, c = b_;
      }
    } else {
      a = c_, b = a_, c = b_;
    }
  } else {
    if (b_ < c_) {
      if (a_ < c_) {
        a = b_, b = a_, c = c_;
      } else {
        a = b_, b = c_, c = a_;
      }
    } else {
      a = c_, b = b_, c = a_;
    }
  }

  // see: http://en.wikipedia.org/wiki/Heron%27s_formula#Numerical_stability
  // see also: http://www.eecs.berkeley.edu/~wkahan/Triangle.pdf

  // the tolerance value is determined empirically
  if ((a + (b - c)) < c * 1e-14) return 0.0;

  double t = tan(0.25 * (a + (b + c))) * tan(0.25 * (c - (a - b))) *
             tan(0.25 * (c + (a - b))) * tan(0.25 * (a + (b - c)));

  return fabs(4.0 * atan(sqrt(fabs(t))));
}

/* ----------------------------------- */

static inline int compute_norm_vector(double a[], double b[], double norm[]) {

  crossproduct_ld(a, b, norm);

  double scale = sqrt(norm[0] * norm[0] + norm[1] * norm[1] + norm[2] * norm[2]);

  if (scale == 0) return 1;

  scale = 1.0 / scale;

  norm[0] *= scale;
  norm[1] *= scale;
  norm[2] *= scale;

  return 0;
}

static double
lat_edge_correction(double base_point[3], double a[3], double b[3],
                    double lon_a, double lon_b) {

  double const tol = 1e-8;

  if (fabs(a[2] - b[2]) > tol)
    yac_internal_abort_message("ERROR: latitude of both corners is not identical\n",
                                __FILE__, __LINE__);

  double h = fabs(a[2]);

  // if we are at the equator or at a pole
  if (h < tol || fabs(1.0 - h) < tol)
    return 0.0;

  double lat_area = fabs((1.0 - h) * get_angle(lon_a, lon_b));

  double pole[3] = {0, 0, (a[2] >= 0.0)?1:-1};
  double gc_area = tri_area(a, b, pole);

  double correction = MAX(lat_area - gc_area, 0.0);

  double middle_lat[3] = {a[0]+b[0], a[1]+b[1], a[2]};
  double scale = sqrt(middle_lat[0]*middle_lat[0]+middle_lat[1]*middle_lat[1]);

  if (scale == 0) yac_internal_abort_message("internal error", __FILE__, __LINE__);

  scale = sqrt(1.0 - a[2]*a[2]) / scale;

  middle_lat[0] *= scale;
  middle_lat[1] *= scale;

  double norm_ab[3];

  // if the angle between a and b is to small to compute a norm vector
  if (compute_norm_vector(a, b, norm_ab)) return 0.0;

  double scalar_base = scalar_product(norm_ab, base_point);
  double scalar_middle_lat = scalar_product(norm_ab, middle_lat);

  // if the base point is on the same plane as a and b
  if (fabs(scalar_base) < 1e-11) {

    double norm_middle[3];
    double pole[3] = {0,0,((a[2]>0)?1:-1)};

    if (compute_norm_vector(middle_lat, pole, norm_middle)) return 0.0;

    double scalar_a = scalar_product(norm_middle, a);

    if (scalar_a > 0)
      return correction;
    else
      return - correction;

  } else {

    if (scalar_middle_lat < 0)
      return correction;
    else
      return - correction;
  }
}

double yac_pole_area ( struct grid_cell cell ) {

  double area = 0.0;

  int M = cell.num_corners;

  if (M < 2) return 0.0;

  int closer_to_south_pole = cell.coordinates_y[0] < 0;

  double pole_vec[3] = {0, 0, (closer_to_south_pole)?-1:1};

  // it would also be possible to use the equator instead
  // of the poles as the baseline
  // this could be used as special case for cell close
  // the equator (were the other method is probably most
  // inaccurate)

  for (int i = 0; i < M; ++i) {

    // if one of the points it at the pole
    if (fabs(fabs(cell.coordinates_y[i]) - M_PI_2) < 1e-12) continue;
    if (fabs(fabs(cell.coordinates_y[(i+1)%M]) - M_PI_2) < 1e-12) {
      ++i; // we can skip the next edge
      continue;
    }

    if (cell.edge_type[i] == GREAT_CIRCLE || cell.edge_type[i] == LON_CIRCLE) {

      double * a;
      double * b;

      a = cell.coordinates_xyz + i * 3;
      b = cell.coordinates_xyz + ((i+1)%M) * 3;

      double edge_direction = a[0]*b[1]-a[1]*b[0]; // 3. component of cross product

      // if the edge is nearly on a circle of longitude
      if (fabs(edge_direction) < 1e-12) continue;

      double tmp_area = tri_area(a, b, pole_vec);

      // or the other way round
      if (edge_direction > 0) area -= tmp_area;
      else                    area += tmp_area;

    } else if (cell.edge_type[i] == LAT_CIRCLE) {

      // the area of a sphere cap is:
      // A = 2 * PI * r * h (where r == 1 and h == 1 - cos(d_lat))
      // scaled with the longitude angle this is:
      // A' = (d_lon / (2 * PI)) * A
      // => A' = d_lon * (1 - cos(d_lat))

      double d_lon = get_angle(cell.coordinates_x[i], cell.coordinates_x[(i+1)%M]);
      double d_lat = M_PI_2;

      if (closer_to_south_pole)
        d_lat += cell.coordinates_y[i];
      else
        d_lat -= cell.coordinates_y[i];

      double h = 1 - cos(d_lat);

      area += d_lon * h;

    } else {

      yac_internal_abort_message("ERROR: unsupported edge type\n", __FILE__, __LINE__);
    }
  }
  // return fabs(area) * EarthRadius * EarthRadius;
  return fabs(area);
}

double yac_planar_3dcell_area (struct grid_cell cell) {

 /*
  * source code is originally based on http://gaim.umbc.edu/2010/06/03/polygon-area/
  *
  */

  double area = 0.0;
  double norm[3] = {0,0,0};

  if (cell.num_corners < 3) return area;

  for ( int i0=cell.num_corners-1, i1=0;  i1<cell.num_corners; i0=i1, ++i1) {
    norm[0] += cell.coordinates_xyz[1+i0*3]*cell.coordinates_xyz[2+i1*3] - cell.coordinates_xyz[1+i1*3]*cell.coordinates_xyz[2+i0*3];
    norm[1] += cell.coordinates_xyz[2+i0*3]*cell.coordinates_xyz[0+i1*3] - cell.coordinates_xyz[2+i1*3]*cell.coordinates_xyz[0+i0*3];
    norm[2] += cell.coordinates_xyz[0+i0*3]*cell.coordinates_xyz[1+i1*3] - cell.coordinates_xyz[0+i1*3]*cell.coordinates_xyz[1+i0*3];
  };

  return area = 0.5 * sqrt(norm[0]*norm[0] + norm[1]*norm[1] + norm[2]*norm[2] );

}

 /*
  * source code is originally based on code by Robert Oehmke of Earth System Modeling
  * Framework (www.earthsystemmodeling.org)
  *
  * it has be extended to support YAC data structures and two types of
  * grid cell edges (great circle and circle of latitude)
  *
  */
double yac_huiliers_area (struct grid_cell cell) {

  if (cell.num_corners < 2) return 0;

  int lat_flag = 0;

  for (int i = 0; i < cell.num_corners; i++)
    lat_flag |= cell.edge_type[i] == LAT_CIRCLE;

  if (cell.num_corners == 3 && !lat_flag)
    return fabs(tri_area(cell.coordinates_xyz + 0*3,
                         cell.coordinates_xyz + 1*3,
                         cell.coordinates_xyz + 2*3));
    //         * EarthRadius * EarthRadius;

  // sum areas around cell
  double area = 0.0;

  for (int i = 2; i < cell.num_corners; ++i) {

    double tmp_area = tri_area(cell.coordinates_xyz + 0*3,
                               cell.coordinates_xyz + (i-1)*3,
                               cell.coordinates_xyz + i * 3);

    double norm[3];

    if (compute_norm_vector(cell.coordinates_xyz + (i-1) * 3,
                            cell.coordinates_xyz + i * 3, norm)) continue;

    double scalar_base = scalar_product(norm, cell.coordinates_xyz + 0*3);

    if (scalar_base > 0)
      area += tmp_area;
    else
      area -= tmp_area;
  }

  // if there is at least one latitude circle edge
  if (lat_flag) {

    for (int i = 0; i < cell.num_corners; ++i) {

      if (cell.edge_type[i] == LAT_CIRCLE) {

        int i_ = (i+1)%cell.num_corners;

        area += lat_edge_correction(cell.coordinates_xyz + 0 * 3,
                                    cell.coordinates_xyz + i * 3,
                                    cell.coordinates_xyz + i_ * 3,
                                    cell.coordinates_x[i],
                                    cell.coordinates_x[i_]);
      }
    }
  }

  return fabs(area);
  // return fabs(area) * EarthRadius * EarthRadius;
}

/* ----------------------------------- */

double yac_partial_area ( double a_lon, double a_lat,
                          double b_lon, double b_lat,
                          double c_lon, double c_lat ) {

  double theta;
  double angle_f;
  double angle_b;

  angle_f = inner_angle ( a_lat, a_lon, b_lat, b_lon );
  angle_b = inner_angle ( a_lat, a_lon, c_lat, c_lon );

  theta = angle_b - angle_f;

  if ( theta < 0.0 ) theta = theta + M_PI + M_PI;

  return theta;
}

/* ----------------------------------- */

static inline double inner_angle ( double plat, double plon,
                                   double qlat, double qlon ) {

  double t = sin((qlon-plon))*cos(qlat);

  double b = sin(qlat)*cos(plat)
           - cos(qlat)*sin(plat) * cos((qlon-plon));

  return atan2(b,t);
}

/* ----------------------------------- */

static inline double scalar_product(double a[], double b[]) {
  return (a[0] * b[0] + a[1] * b[1] + a[2] * b[2]);
}

/* ----------------------------------- */
