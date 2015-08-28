/**
 * @file intersection.c
 * @brief Set of functions to determine the intersection between two edges
 *
 * very interesting literature:
 * - http://geospatialmethods.org/spheres/GCIntersect.html
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

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include <stdlib.h>
#include <math.h>
#include <stdio.h>

#include "utils.h"
#include "geometry.h"

#if defined(CDO) && !defined(HAVE_SQRTL)
#define sqrtl(x) (sqrt((double)x))
#endif

// angle tolerance
const double yac_angle_tol = 1e-9;
static double const tol = 1.0e-10;

static int vector_is_between (double a[], double b[], double p[],
                              double * angle_ab, double dot_ab) {

   double dot_ap = a[0]*p[0] + a[1]*p[1] + a[2]*p[2];
   double dot_pb = p[0]*b[0] + p[1]*b[1] + p[2]*b[2];

   // catches most obvious false-cases
   if ((dot_ap < dot_ab - 0.1) ||
       (dot_pb < dot_ab - 0.1))
   return 0;

   if (*angle_ab < 0)
      *angle_ab = get_vector_angle(a, b);

/* determines whether p is between a and b
   (a, b, p are in the same plane AB)
   angle_ab is the angle between a and b */

   return fabs(get_vector_angle(a, p) +
               get_vector_angle(b, p) -
               *angle_ab) < yac_angle_tol;
}

static int vector_is_between_lat (double a[], double b[], double p[]) {

/* determines whether p is between a and b
   (a, b, p have the same latitude)*/

/*  I. a_0 * alpha + b_0 * beta = p_0
   II. a_1 * alpha + b_1 * beta = p_1
   
   if alpha > 0 and beta > 0 -> p is between a and b */

   if (fabs(fabs(a[2]) - 1.0) < tol) return 1;

   int flag = fabs(a[0]) > fabs(a[1]);
   double a_0 = a[flag], a_1 = a[flag^1];
   double b_0 = b[flag], b_1 = b[flag^1];
   double p_0 = p[flag], p_1 = p[flag^1];

   double temp = b_0 - (b_1 * a_0) / a_1;

   // if a and b are nearly identical
   if (fabs(temp) < tol) return (fabs(a_0 - p_0) < tol) &&
                                (fabs(a_1 - p_1) < tol);

   double beta = (p_0 - (p_1 * a_0) / a_1) / temp;

   if (beta < -tol) return 0;

   double alpha = (p_1 - b_1 * beta) / a_1;

   return alpha > -tol;
}

/** \brief compute the intersection points of two great circles
  *
  * if p and q are != NULL they contain the intersection points
  *
  * the return value is :
  *    -  0 if the intersection points are neither between (a and b) or (c and d)
  *    - -1 if an error occurred
  *    - 1st bit will be set if p is between a and b
  *    - 2nd bit will be set if q is between a and b
  *    - 3rd bit will be set if p is between c and d
  *    - 4th bit will be set if q is between c and d
  *    - 5th bit will be set if both great circles are identically
  *
  * based on
  * - http://www.geoclub.de/viewtopic.php?f=54&t=29689
  **/

 int yac_gcxgc (struct edge edge_a, struct edge edge_b,
                struct point * p, struct point * q) {

   double a[3], b[3], c[3], d[3], p_[3], q_[3];

   LLtoXYZ( edge_a.points[0].lon, edge_a.points[0].lat, a);
   LLtoXYZ( edge_a.points[1].lon, edge_a.points[1].lat, b);
   LLtoXYZ( edge_b.points[0].lon, edge_b.points[0].lat, c);
   LLtoXYZ( edge_b.points[1].lon, edge_b.points[1].lat, d);

   int ret_val =  yac_gcxgc_vec(a, b, c, d, p_, q_);

   if (p != NULL) XYZtoLL(p_, &(p->lon), &(p->lat));
   if (q != NULL) XYZtoLL(q_, &(q->lon), &(q->lat));

   return ret_val;
}

/** \brief compute the intersection points of two great circles
  *
  * if p and q are != NULL they contain the intersection points
  *
  * the return value is :
  *    -  0 if the intersection points are neither between (a and b) or (c and d)
  *    - -1 if an error occurred
  *    - 1st bit will be set if p is between a and b
  *    - 2nd bit will be set if q is between a and b
  *    - 3rd bit will be set if p is between c and d
  *    - 4th bit will be set if q is between c and d
  *    - 5th bit will be set if both great circles are identically
  *
  * based on
  * - http://www.geoclub.de/viewtopic.php?f=54&t=29689
  **/
 int yac_gcxgc_vec (double a[3], double b[3], double c[3], double d[3],
                    double p[3], double q[3]) {

   double length_cross_ab, length_cross_cd, length_cross_abxcd;
   long double cross_ab[3], cross_cd[3], cross_abxcd[3];

   // p' = (a X b) X (c X d)
   // p  = p' / length(p')
   {
      long double a_[3] = {a[0], a[1], a[2]};
      long double b_[3] = {b[0], b[1], b[2]};
      long double c_[3] = {c[0], c[1], c[2]};
      long double d_[3] = {d[0], d[1], d[2]};

      cross_ab[0] = a_[1] * b_[2] - a_[2] * b_[1];
      cross_ab[1] = a_[2] * b_[0] - a_[0] * b_[2];
      cross_ab[2] = a_[0] * b_[1] - a_[1] * b_[0];
      cross_cd[0] = c_[1] * d_[2] - c_[2] * d_[1];
      cross_cd[1] = c_[2] * d_[0] - c_[0] * d_[2];
      cross_cd[2] = c_[0] * d_[1] - c_[1] * d_[0];
      cross_abxcd[0] = cross_ab[1] * cross_cd[2] - cross_ab[2] * cross_cd[1];
      cross_abxcd[1] = cross_ab[2] * cross_cd[0] - cross_ab[0] * cross_cd[2];
      cross_abxcd[2] = cross_ab[0] * cross_cd[1] - cross_ab[1] * cross_cd[0];
   }

   length_cross_ab = sqrtl(cross_ab[0] * cross_ab[0] +
                           cross_ab[1] * cross_ab[1] +
                           cross_ab[2] * cross_ab[2]);
   length_cross_cd = sqrtl(cross_cd[0] * cross_cd[0] +
                           cross_cd[1] * cross_cd[1] +
                           cross_cd[2] * cross_cd[2]);
   length_cross_abxcd = sqrtl(cross_abxcd[0] * cross_abxcd[0] +
                              cross_abxcd[1] * cross_abxcd[1] +
                              cross_abxcd[2] * cross_abxcd[2]);

   if (length_cross_ab < tol) {

      if (p != NULL) p[0] = a[0], p[1] = a[1], p[2] = a[2];
      if (q != NULL) q[0] = -a[0], q[1] = -a[1], q[2] = -a[2];

      if (length_cross_cd < tol) {

         double cross_ac[3], angle_ac;

         // we only care about accuracy of the angle in case the angle is very
         // small or nearly 180°, therefore we use directly the cross product
         // instead of the routine get_vector_angle
         crossproduct_ld(a, c, cross_ac);
         angle_ac = asin(cross_ac[0] * cross_ac[0] +
                         cross_ac[1] * cross_ac[1] +
                         cross_ac[2] * cross_ac[2]);

         // if points are identically or have an angle of 180°
         if (fabs(angle_ac) < tol) {

            double dot_ac = a[0]*c[0] + a[1]*c[1] + a[2]*c[2];

            return (1 + ((dot_ac > 0)?4:8));

         } else

            return -1;

      } else {

         // dot product is most accurate for angles around half PI
         // (for angles very close to half PI: alpha = fabs(acos(alpha))
         double angle = fabs(cross_cd[0] * a[0] +
                             cross_cd[1] * a[1] +
                             cross_cd[2] * a[2]) / length_cross_cd ;

         // if ab is on the plane of cd
         if (fabs(angle) < tol) {

            int result = 1;
            double angle_cd = -1;
            double dot_cd = c[0]*d[0] + c[1]*d[1] + c[2]*d[2];

            if (vector_is_between(c, d, a, &angle_cd, dot_cd))
               result |= 1 << 2;
            double tvec[] = {-a[0], -a[1], -a[2]};
            if (vector_is_between(c, d, tvec, &angle_cd, dot_cd))
	       result |= 1 << 3;

            return result;
         }

         return -1;
      }

   } else {

      if (length_cross_cd < tol) {

         if (p != NULL) p[0] = c[0], p[1] = c[1], p[2] = c[2];
         if (q != NULL) q[0] = -c[0], q[1] = -c[1], q[2] = -c[2];

         // dot product is most accurate for angles around half PI
         // (for angles very close to half PI: alpha = fabs(acos(alpha))
         double angle = fabs(cross_ab[0] * c[0] +
                             cross_ab[1] * c[1] +
                             cross_ab[2] * c[2]) / length_cross_ab;

         // if cd is on the plane of ab
         if (fabs(angle) < tol) {

            int result = 4;

            double angle_ab = -1;
            double dot_ab = a[0]*b[0] + a[1]*b[1] + a[2]*b[2];

            if (vector_is_between(a, b, c, &angle_ab, dot_ab))
               result |= 1 << 0;
	    double tvec[] = {-c[0], -c[1], -c[2]};
            if (vector_is_between(a, b, tvec, &angle_ab, dot_ab))
	       result |= 1 << 1;

            return result;
         }

         return -1;
      }
   }

   // if both great circles are nearly identically
   if (length_cross_abxcd / (length_cross_ab * length_cross_cd) < tol) {

      int ret_value = 1 << 4;

      int a_between_cd, b_between_cd, c_between_ab, d_between_ab;
      double angle_ab = -1;
      double angle_cd = -1;
      double dot_ab = a[0]*b[0] + a[1]*b[1] + a[2]*b[2];
      double dot_cd = c[0]*d[0] + c[1]*d[1] + c[2]*d[2];

      a_between_cd = vector_is_between(c, d, a, &angle_cd, dot_cd) << 0;
      b_between_cd = vector_is_between(c, d, b, &angle_cd, dot_cd) << 1;
      c_between_ab = vector_is_between(a, b, c, &angle_ab, dot_ab) << 2;
      d_between_ab = vector_is_between(a, b, d, &angle_ab, dot_ab) << 3;

      switch (a_between_cd + b_between_cd + c_between_ab + d_between_ab) {

         case (0):
            p[0] = a[0], p[1] = a[1], p[2] = a[2];
            q[0] = b[0], q[1] = b[1], q[2] = b[2];
            ret_value |= 1 + 2;
            return ret_value;
         case (1+2):
         case (1+2+4):
         case (1+2+8):
            p[0] = a[0], p[1] = a[1], p[2] = a[2];
            q[0] = b[0], q[1] = b[1], q[2] = b[2];
            break;
         case (1+4):
         case (1+4+8):
            p[0] = a[0], p[1] = a[1], p[2] = a[2];
            q[0] = c[0], q[1] = c[1], q[2] = c[2];
            break;
         case (1+8):
            p[0] = a[0], p[1] = a[1], p[2] = a[2];
            q[0] = d[0], q[1] = d[1], q[2] = d[2];
            break;
         case (2+4):
         case (2+4+8):
            p[0] = b[0], p[1] = b[1], p[2] = b[2];
            q[0] = c[0], q[1] = c[1], q[2] = c[2];
            break;
         case (2+8):
            p[0] = b[0], p[1] = b[1], p[2] = b[2];
            q[0] = d[0], q[1] = d[1], q[2] = d[2];
            break;
         case (4+8):
            p[0] = c[0], p[1] = c[1], p[2] = c[2];
            q[0] = d[0], q[1] = d[1], q[2] = d[2];
            break;
         case (1+2+4+8):
            p[0] = a[0], p[1] = a[1], p[2] = a[2];
            q[0] = a[0], q[1] = a[1], q[2] = a[2];
            break;
         default:
            yac_internal_abort_message("internal error", __FILE__, __LINE__);
      }

      ret_value |= 1 + 2 + 4 + 8;

      return ret_value;
   }

    long double scale = 1.0l / length_cross_abxcd;
    // determine p and q
    double p_[3], q_[3];

    p_[0]= cross_abxcd[0] * scale;
    p_[1]= cross_abxcd[1] * scale;
    p_[2]= cross_abxcd[2] * scale;

    q_[0]=-p_[0];
    q_[1]=-p_[1];
    q_[2]=-p_[2];

    // set p and q
    if (p != 0) {
       p[0] = p_[0];
       p[1] = p_[1];
       p[2] = p_[2];
    }
    if (q != 0) {
       q[0] = q_[0];
       q[1] = q_[1];
       q[2] = q_[2];
    }

    int result;
    double angle_ab = -1;
    double angle_cd = -1;
    double dot_ab = a[0]*b[0] + a[1]*b[1] + a[2]*b[2];
    double dot_cd = c[0]*d[0] + c[1]*d[1] + c[2]*d[2];

    result = 0;
    if (vector_is_between(a, b, p_, &angle_ab, dot_ab)) result |= 1 << 0;
    if (vector_is_between(a, b, q_, &angle_ab, dot_ab)) result |= 1 << 1;
    if (vector_is_between(c, d, p_, &angle_cd, dot_cd)) result |= 1 << 2;
    if (vector_is_between(c, d, q_, &angle_cd, dot_cd)) result |= 1 << 3;

    return result;
}

static
int gcxgc_vec_ (double a[3], double b[3], double c[3], double d[3]) {

   double length_cross_ab, length_cross_cd, length_cross_abxcd;
   long double cross_ab[3], cross_cd[3], cross_abxcd[3];

   // p' = (a X b) X (c X d)
   // p  = p' / length(p')
   {
      long double a_[3] = {a[0], a[1], a[2]};
      long double b_[3] = {b[0], b[1], b[2]};
      long double c_[3] = {c[0], c[1], c[2]};
      long double d_[3] = {d[0], d[1], d[2]};

      cross_ab[0] = a_[1] * b_[2] - a_[2] * b_[1];
      cross_ab[1] = a_[2] * b_[0] - a_[0] * b_[2];
      cross_ab[2] = a_[0] * b_[1] - a_[1] * b_[0];
      cross_cd[0] = c_[1] * d_[2] - c_[2] * d_[1];
      cross_cd[1] = c_[2] * d_[0] - c_[0] * d_[2];
      cross_cd[2] = c_[0] * d_[1] - c_[1] * d_[0];
      cross_abxcd[0] = cross_ab[1] * cross_cd[2] - cross_ab[2] * cross_cd[1];
      cross_abxcd[1] = cross_ab[2] * cross_cd[0] - cross_ab[0] * cross_cd[2];
      cross_abxcd[2] = cross_ab[0] * cross_cd[1] - cross_ab[1] * cross_cd[0];
   }

   length_cross_ab = sqrtl(cross_ab[0] * cross_ab[0] +
                           cross_ab[1] * cross_ab[1] +
                           cross_ab[2] * cross_ab[2]);
   length_cross_cd = sqrtl(cross_cd[0] * cross_cd[0] +
                           cross_cd[1] * cross_cd[1] +
                           cross_cd[2] * cross_cd[2]);
   length_cross_abxcd = sqrtl(cross_abxcd[0] * cross_abxcd[0] +
                              cross_abxcd[1] * cross_abxcd[1] +
                              cross_abxcd[2] * cross_abxcd[2]);

   if (length_cross_ab < tol) {

      if (length_cross_cd < tol) {

         double dot_ac = a[0]*c[0] + a[1]*c[1] + a[2]*c[2];

         // if the angle between a and c is bigger than ~45°
         if (dot_ac < 0.75) return 0;

         double cross_ac[3], angle_ac;

         // we only care about accuracy of the angle in case the angle is very
         // small or nearly 180°, therefore we use directly the cross product
         // instead of the routine get_vector_angle
         crossproduct_ld(a, c, cross_ac);
         angle_ac = asin(cross_ac[0] * cross_ac[0] +
                         cross_ac[1] * cross_ac[1] +
                         cross_ac[2] * cross_ac[2]);

         // if points are identically or have an angle of 180°
         return fabs(angle_ac) < tol;

      } else {

         // dot product is most accurate for angles around half PI
         // (for angles very close to half PI: alpha = fabs(acos(alpha))
         double angle = fabs(cross_cd[0] * a[0] +
                             cross_cd[1] * a[1] +
                             cross_cd[2] * a[2]) / length_cross_cd ;

         // if ab is not on the plane of cd
         if (fabs(angle) >= tol)
            return 0;

         double angle_cd = -1;
         double dot_cd = c[0]*d[0] + c[1]*d[1] + c[2]*d[2];

         return vector_is_between(c, d, a, &angle_cd, dot_cd);
      }

   } else {

      if (length_cross_cd < tol) {

         // dot product is most accurate for angles around half PI
         // (for angles very close to half PI: alpha = fabs(acos(alpha))
         double angle = fabs(cross_ab[0] * c[0] +
                             cross_ab[1] * c[1] +
                             cross_ab[2] * c[2]) / length_cross_ab;

         // if cd is not on the plane of ab
         if (fabs(angle) >= tol)
            return 0;

         double angle_ab = -1;
         double dot_ab = a[0]*b[0] + a[1]*b[1] + a[2]*b[2];

         return vector_is_between(a, b, c, &angle_ab, dot_ab);
      }
   }

   // if both great circles are nearly identically
   if (length_cross_abxcd / (length_cross_ab * length_cross_cd) < tol) {

      double angle_ab = -1;
      double angle_cd = -1;
      double dot_ab = a[0]*b[0] + a[1]*b[1] + a[2]*b[2];
      double dot_cd = c[0]*d[0] + c[1]*d[1] + c[2]*d[2];

      return vector_is_between(c, d, a, &angle_cd, dot_cd) ||
             vector_is_between(c, d, b, &angle_cd, dot_cd) ||
             vector_is_between(a, b, c, &angle_ab, dot_ab) ||
             vector_is_between(a, b, d, &angle_ab, dot_ab);
   }

    long double scale = 1.0l / length_cross_abxcd;
    // determine p and q
    double p_[3];

    p_[0]= cross_abxcd[0] * scale;
    p_[1]= cross_abxcd[1] * scale;
    p_[2]= cross_abxcd[2] * scale;

    double angle_ab = -1;
    double angle_cd = -1;
    double dot_ab = a[0]*b[0] + a[1]*b[1] + a[2]*b[2];
    double dot_cd = c[0]*d[0] + c[1]*d[1] + c[2]*d[2];

    if (vector_is_between(a, b, p_, &angle_ab, dot_ab) &&
        vector_is_between(c, d, p_, &angle_cd, dot_cd))
       return 1;

    double q_[3];
    q_[0]=-p_[0];
    q_[1]=-p_[1];
    q_[2]=-p_[2];

    return vector_is_between(a, b, q_, &angle_ab, dot_ab) &&
           vector_is_between(c, d, q_, &angle_cd, dot_cd);

}

/** \brief compute the intersection point two circles of latitude
 *
 * compute the intersection points of two circle of latitude
 * if p and q are != NULL they contain the intersection points
 * the return value is:
 *      - 0 if the intersection points are neither between (a and b) or (c and d)
 *      - -1 if an error occurred
 *      - 1st bit will be set if p is between a and b
 *      - 2nd bit will be set if q is between a and b
 *      - 3rd bit will be set if p is between c and d
 *      - 4th bit will be set if q is between c and d
 *      - 5th bit will be set if both edges are on the same circle of latitude
 **/
int yac_latcxlatc (struct edge edge_a, struct edge edge_b,
                   struct point * p, struct point * q) {

   // two circles of latitude can only intersect if they are on the same latitude
   if (fabs(edge_a.points[0].lat - edge_b.points[0].lat) > tol)
      return -1;

   int ret_value = 1 << 4;

   // if the circles are on a pole
   if (fabs(fabs(edge_a.points[0].lat) - M_PI_2) < tol) {

      if (p != NULL) *p = edge_a.points[0];
      if (q != NULL) *q = edge_a.points[0];

      ret_value |= 1 + 4;

      return ret_value;
   }

   // check whether the two circles overlap
   
   double angle_ab;
   double angle_ac;
   double angle_ad;
   double angle_bc;
   double angle_bd;
   double angle_cd;
   double angle_cb;

   angle_ab = fabs(get_angle(edge_a.points[0].lon, edge_a.points[1].lon));
   angle_ac = fabs(get_angle(edge_a.points[0].lon, edge_b.points[0].lon));
   angle_ad = fabs(get_angle(edge_a.points[0].lon, edge_b.points[1].lon));
   angle_bc = fabs(get_angle(edge_a.points[1].lon, edge_b.points[0].lon));
   angle_bd = fabs(get_angle(edge_a.points[1].lon, edge_b.points[1].lon));
   angle_cd = fabs(get_angle(edge_b.points[0].lon, edge_b.points[1].lon));
   angle_cb = fabs(get_angle(edge_b.points[0].lon, edge_a.points[1].lon));

   int a_between_cd, b_between_cd, c_between_ab, d_between_ab;

   a_between_cd = ((angle_cd + tol) > (angle_ac + angle_ad)) << 0;
   b_between_cd = ((angle_cd + tol) > (angle_bc + angle_bd)) << 1;
   c_between_ab = ((angle_ab + tol) > (angle_ac + angle_cb)) << 2;
   d_between_ab = ((angle_ab + tol) > (angle_ad + angle_bd)) << 3;

   switch (a_between_cd + b_between_cd + c_between_ab + d_between_ab) {
      case (0):
      {
         if (p != NULL) *p = edge_a.points[0];
         if (q != NULL) *q = edge_a.points[1];
         ret_value |= 1 + 2;
         return ret_value;
      }
      case (1+2):
      case (1+2+4):
      case (1+2+8):
         if (p != NULL) *p = edge_a.points[0];
         if (q != NULL) *q = edge_a.points[1];
         break;
      case (1+4):
      case (1+4+8):
         if (p != NULL) *p = edge_a.points[0];
         if (q != NULL) *q = edge_b.points[0];
         break;
      case (1+8):
         if (p != NULL) *p = edge_a.points[0];
         if (q != NULL) *q = edge_b.points[1];
         break;
      case (2+4):
      case (2+4+8):
         if (p != NULL) *p = edge_a.points[1];
         if (q != NULL) *q = edge_b.points[0];
         break;
      case (2+8):
         if (p != NULL) *p = edge_a.points[1];
         if (q != NULL) *q = edge_b.points[1];
         break;
      case (4+8):
         if (p != NULL) *p = edge_b.points[0];
         if (q != NULL) *q = edge_b.points[1];
         break;
      case (1+2+4+8):
         if (p != NULL) *p = edge_a.points[0];
         if (q != NULL) *q = edge_a.points[0];
         break;
      default:
         yac_internal_abort_message("internal error", __FILE__, __LINE__);
   }

   if (angle_ab < tol || angle_cd < tol)
      ret_value |= 1 + 4;
   else
      ret_value |= 1 + 2 + 4 + 8;

   return ret_value;
}

/** \brief compute the intersection point two circles of latitude
 *
 * compute the intersection points of two circle of latitude
 * if p and q are != NULL they contain the intersection points
 * the return value is:
 *      - 0 if the intersection points are neither between (a and b) or (c and d)
 *      - -1 if an error occurred
 *      - 1st bit will be set if p is between a and b
 *      - 2nd bit will be set if q is between a and b
 *      - 3rd bit will be set if p is between c and d
 *      - 4th bit will be set if q is between c and d
 *      - 5th bit will be set if both edges are on the same circle of latitude
 **/
int yac_latcxlatc_vec (double a[3], double b[3], double c[3], double d[3],
                       double p[3], double q[3]) {

   // two circles of latitude can only intersect if they are on the same latitude
   if (fabs(a[2] - c[2]) > tol)
      return -1;

   int result = 16;

   int a_between_cd, b_between_cd, c_between_ab, d_between_ab;

   a_between_cd = vector_is_between_lat(c, d, a);
   b_between_cd = vector_is_between_lat(c, d, b);
   c_between_ab = vector_is_between_lat(a, b, c);
   d_between_ab = vector_is_between_lat(a, b, d);

   if (a_between_cd && b_between_cd && c_between_ab && d_between_ab) {

      p[0] = a[0], p[1] = a[1], p[2] = a[2];
      q[0] = a[0], q[1] = a[1], q[2] = a[2];

      result |= 1 + 4;

   } else if (a_between_cd) {

      p[0] = a[0], p[1] = a[1], p[2] = a[2];

      result |= 1 + 2 + 4 + 8;

      if (b_between_cd) q[0] = b[0], q[1] = b[1], q[2] = b[2];
      else if (c_between_ab) q[0] = c[0], q[1] = c[1], q[2] = c[2];
      else if (d_between_ab) q[0] = d[0], q[1] = d[1], q[2] = d[2];
      else yac_internal_abort_message("internal error", __FILE__, __LINE__);

   } else if (b_between_cd) {

      p[0] = b[0], p[1] = b[1], p[2] = b[2];

      result |= 1 + 2 + 4 + 8;

      if (c_between_ab) q[0] = c[0], q[1] = c[1], q[2] = c[2];
      else if (d_between_ab) q[0] = d[0], q[1] = d[1], q[2] = d[2];
      else yac_internal_abort_message("internal error", __FILE__, __LINE__);

   } else if (c_between_ab && d_between_ab) {

      p[0] = c[0], p[1] = c[1], p[2] = c[2];
      q[0] = d[0], q[1] = d[1], q[2] = d[2];

      result |= 1 + 2 + 4 + 8;

   } else {

      p[0] = a[0], p[1] = a[1], p[2] = a[2];
      q[0] = b[0], q[1] = b[1], q[2] = b[2];

      result |= 1 + 2;
   }

   return result;
}

static int latcxlatc_(struct edge edge_a, struct edge edge_b) {

   // if both edges are on the same circle of latitude
   if (fabs(edge_a.points[0].lat - edge_b.points[0].lat) > tol)
      return 0;

   // if both edges are on the pole
   if (fabs(fabs(edge_a.points[0].lat) - M_PI_2) < tol)
      return 1;

   double angle_ab = fabs(get_angle(edge_a.points[0].lon,
                                    edge_a.points[1].lon)) + tol;
   double angle_cd = fabs(get_angle(edge_b.points[0].lon,
                                    edge_b.points[1].lon)) + tol;
   double angle_ac = fabs(get_angle(edge_a.points[0].lon,
                                    edge_b.points[0].lon));
   double angle_ad = fabs(get_angle(edge_a.points[0].lon,
                                    edge_b.points[1].lon));
   double angle_bc = fabs(get_angle(edge_a.points[1].lon,
                                    edge_b.points[0].lon));
   double angle_bd = fabs(get_angle(edge_a.points[1].lon,
                                    edge_b.points[1].lon));

   return ((angle_ac + angle_bc) < angle_ab) ||
          ((angle_ad + angle_bd) < angle_ab) ||
          ((angle_ac + angle_ad) < angle_cd) ||
          ((angle_bc + angle_bd) < angle_cd);
}

/** \brief compute the intersection point two circles of longitude
 *
 * compute the intersection points of two circle of longitude
 * if p and q are != NULL they contain the intersection points
 * the return value is:
 *      - 0 if the intersection points are neither between (a and b) or (c and d)
 *      - -1 if an error occurred
 *      - 1st bit will be set if p is between a and b
 *      - 2nd bit will be set if q is between a and b
 *      - 3rd bit will be set if p is between c and d
 *      - 4th bit will be set if q is between c and d
 *      - 5th bit will be set if both edges are on the same circle of longitude
 **/
int yac_loncxlonc (struct edge edge_a, struct edge edge_b,
                   struct point * p, struct point * q) {

   double angle_ab = fabs(get_angle(edge_a.points[0].lon, edge_a.points[1].lon));
   double angle_ac = fabs(get_angle(edge_a.points[0].lon, edge_b.points[0].lon));
   double angle_ad = fabs(get_angle(edge_a.points[0].lon, edge_b.points[1].lon));
   double angle_cd = fabs(get_angle(edge_b.points[0].lon, edge_b.points[1].lon));

   int ret_value = 0;

   // if both edges are on identical circles of longitude
   if ((angle_ac < tol) || (fabs(M_PI - angle_ac) < tol)) {

      ret_value |= 16;

      // rotate circles to the equator (lat == 0)

      double lon_a, lon_b, lon_c, lon_d;

      lon_a = edge_a.points[0].lat;
      lon_b = (angle_ab < tol)?(edge_a.points[1].lat):(M_PI - edge_a.points[1].lat);
      lon_c = (angle_ac < tol)?(edge_b.points[0].lat):(M_PI - edge_b.points[0].lat);
      lon_d = (angle_ad < tol)?(edge_b.points[1].lat):(M_PI - edge_b.points[1].lat);

      double angle_ab_ = fabs(get_angle(lon_a, lon_b));
      double angle_ac_ = fabs(get_angle(lon_a, lon_c));
      double angle_ad_ = fabs(get_angle(lon_a, lon_d));
      double angle_cb_ = fabs(get_angle(lon_c, lon_b));
      double angle_cd_ = fabs(get_angle(lon_c, lon_d));
      double angle_bc_ = fabs(get_angle(lon_b, lon_c));
      double angle_bd_ = fabs(get_angle(lon_b, lon_d));

      int a_between_cd, b_between_cd, c_between_ab, d_between_ab;

      a_between_cd = ((angle_cd_ + tol) > (angle_ac_ + angle_ad_)) << 0;
      b_between_cd = ((angle_cd_ + tol) > (angle_bc_ + angle_bd_)) << 1;
      c_between_ab = ((angle_ab_ + tol) > (angle_ac_ + angle_cb_)) << 2;
      d_between_ab = ((angle_ab_ + tol) > (angle_ad_ + angle_bd_)) << 3;

      switch (a_between_cd + b_between_cd + c_between_ab + d_between_ab) {

         case (0):
         {
            double temp = (edge_a.points[0].lat > 0)?M_PI_2:-M_PI_2;
            if (p != NULL) *p = (struct point){.lon = 0, .lat = temp};
            if (q != NULL) *q = (struct point){.lon = 0, .lat = -temp};
            // if edge a goes across a pole
            if ((angle_ab > tol) || (fabs(edge_a.points[0].lat) > (M_PI_2 - tol)) ||
                                    (fabs(edge_a.points[1].lat) > (M_PI_2 - tol))) {

               ret_value |= 1 << 0;
            }
            // if edge b goes across a pole
            if ((angle_cd > tol) || (fabs(edge_b.points[0].lat) > (M_PI_2 - tol)) ||
                                    (fabs(edge_b.points[1].lat) > (M_PI_2 - tol))) {

               ret_value |= ((edge_b.points[0].lat > 0) ^ (temp < 0))?(1 << 2):(1 << 3);
            }
            return ret_value;
         }
         case (1+2):
         case (1+2+4):
         case (1+2+8):
            if (p != NULL) *p = edge_a.points[0];
            if (q != NULL) *q = edge_a.points[1];
            break;
         case (1+4):
         case (1+4+8):
            if (p != NULL) *p = edge_a.points[0];
            if (q != NULL) *q = edge_b.points[0];
            break;
         case (1+8):
            if (p != NULL) *p = edge_a.points[0];
            if (q != NULL) *q = edge_b.points[1];
            break;
         case (2+4):
         case (2+4+8):
            if (p != NULL) *p = edge_a.points[1];
            if (q != NULL) *q = edge_b.points[0];
            break;
         case (2+8):
            if (p != NULL) *p = edge_a.points[1];
            if (q != NULL) *q = edge_b.points[1];
            break;
         case (4+8):
            if (p != NULL) *p = edge_b.points[0];
            if (q != NULL) *q = edge_b.points[1];
            break;
         case (1+2+4+8):
            if (p != NULL) *p = edge_a.points[0];
            if (q != NULL) *q = edge_a.points[0];
            break;
         default:
            yac_internal_abort_message("internal error", __FILE__, __LINE__);
      }
      ret_value |= 1 + 2 + 4 + 8;

   } else {

      double sign = (edge_a.points[0].lat > 0)?1.0:-1.0;

      if (p != NULL)
         p->lon = edge_a.points[0].lon, p->lat = M_PI_2 * sign;
      if (q != NULL)
         q->lon = edge_a.points[0].lon+M_PI, q->lat = -M_PI_2 * sign;

      // if edge a goes across a pole
      if ((angle_ab > tol) || (fabs(edge_a.points[0].lat) > (M_PI_2 - tol)) ||
                              (fabs(edge_a.points[1].lat) > (M_PI_2 - tol))) {

         ret_value |= 1 << 0;
      }

      // if edge b goes across a pole
      if ((angle_cd > tol) || (fabs(edge_b.points[0].lat) > (M_PI_2 - tol)) ||
                              (fabs(edge_b.points[1].lat) > (M_PI_2 - tol))) {

         ret_value |= ((edge_b.points[0].lat > 0) ^ (sign < 0))?(1 << 2):(1 << 3);
      }

      // if both edges are only points at the poles
      if ((fabs(edge_a.points[0].lat) > (M_PI_2 - tol)) &&
          (fabs(edge_a.points[1].lat) > (M_PI_2 - tol)) &&
          (fabs(edge_b.points[0].lat) > (M_PI_2 - tol)) &&
          (fabs(edge_b.points[1].lat) > (M_PI_2 - tol)))
         ret_value |= 1 << 4;
   }

   return ret_value;
}

/** \brief compute the intersection point two circles of longitude
 *
 * compute the intersection points of two circle of longitude
 * if p and q are != NULL they contain the intersection points
 * the return value is:
 *      - 0 if the intersection points are neither between (a and b) or (c and d)
 *      - -1 if an error occurred
 *      - 1st bit will be set if p is between a and b
 *      - 2nd bit will be set if q is between a and b
 *      - 3rd bit will be set if p is between c and d
 *      - 4th bit will be set if q is between c and d
 *      - 5th bit will be set if both edges are on the same circle of longitude
 **/
int yac_loncxlonc_vec (double a[3], double b[3], double c[3], double d[3],
                       double p[3], double q[3]) {

   int ret_value = 0;

   double cross_ab[3], cross_cd[3];

   crossproduct_ld(a, b, cross_ab);
   crossproduct_ld(c, d, cross_cd);

   double abs_norm_cross_ab[2], abs_norm_cross_cd[2];

   // abs_norm_cross_ab[2] = 0;
   // abs_norm_cross_cd[2] = 0;

   int edge_is_pole_point = 0;

   double * ref_point;

   ref_point = (fabs(a[2]) > fabs(b[2]))?b:a;
   // if both points are at the pole
   if (fabs(ref_point[2]) > 1.0-tol) {
      abs_norm_cross_ab[0] = 1;
      abs_norm_cross_ab[1] = 0;
      edge_is_pole_point = 1;
   } else {
      double scale = 1.0 / sqrt(ref_point[0]*ref_point[0] +
                                ref_point[1]*ref_point[1]);
      abs_norm_cross_ab[0] = ref_point[1] * scale;
      abs_norm_cross_ab[1] = ref_point[0] * scale;
      double max_abs_val = (fabs(abs_norm_cross_ab[0]) >
                            fabs(abs_norm_cross_ab[1]))?(abs_norm_cross_ab[0]):
                                                        (abs_norm_cross_ab[1]);
      if (max_abs_val < 0) {
         abs_norm_cross_ab[0] *= -1.0;
         abs_norm_cross_ab[1] *= -1.0;
      }
   }

   ref_point = (fabs(c[2]) > fabs(d[2]))?d:c;
   // if both points are at the pole
   if (fabs(ref_point[2]) > 1.0-tol) {
      abs_norm_cross_cd[0] = 0;
      abs_norm_cross_cd[1] = 1;
      edge_is_pole_point = 1;
   } else {
      double scale = 1.0 / sqrt(ref_point[0]*ref_point[0] +
                                ref_point[1]*ref_point[1]);
      abs_norm_cross_cd[0] = ref_point[1] * scale;
      abs_norm_cross_cd[1] = ref_point[0] * scale;
      double max_abs_val = (fabs(abs_norm_cross_cd[0]) >
                            fabs(abs_norm_cross_cd[1]))?(abs_norm_cross_cd[0]):
                                                        (abs_norm_cross_cd[1]);
      if (max_abs_val < 0) {
         abs_norm_cross_cd[0] *= -1.0;
         abs_norm_cross_cd[1] *= -1.0;
      }
   }

   double angle_ab = -1;
   double angle_cd = -1;

   // if both edges are on the same circle of longitude
   if (fabs(abs_norm_cross_ab[0] - abs_norm_cross_cd[0]) < tol &&
       fabs(abs_norm_cross_ab[1] - abs_norm_cross_cd[1]) < tol) {

      if (!edge_is_pole_point)
         ret_value |= 16;

      int a_between_cd, b_between_cd, c_between_ab, d_between_ab;
      double dot_ab = a[0]*b[0] + a[1]*b[1] + a[2]*b[2];
      double dot_cd = c[0]*d[0] + c[1]*d[1] + c[2]*d[2];

      a_between_cd = vector_is_between(c, d, a, &angle_cd, dot_cd) << 0;
      b_between_cd = vector_is_between(c, d, b, &angle_cd, dot_cd) << 1;
      c_between_ab = vector_is_between(a, b, c, &angle_ab, dot_ab) << 2;
      d_between_ab = vector_is_between(a, b, d, &angle_ab, dot_ab) << 3;

      switch (a_between_cd + b_between_cd + c_between_ab + d_between_ab) {

         case (0):
            p[0] = 0, p[1] = 0, p[2] = 1;
            q[0] = 0, q[1] = 0, q[2] = -1;
            break;
         case (1+2):
         case (1+2+4):
         case (1+2+8):
            p[0] = a[0], p[1] = a[1], p[2] = a[2];
            q[0] = b[0], q[1] = b[1], q[2] = b[2];
            break;
         case (1+4):
         case (1+4+8):
            p[0] = a[0], p[1] = a[1], p[2] = a[2];
            q[0] = c[0], q[1] = c[1], q[2] = c[2];
            break;
         case (1+8):
            p[0] = a[0], p[1] = a[1], p[2] = a[2];
            q[0] = d[0], q[1] = d[1], q[2] = d[2];
            break;
         case (2+4):
         case (2+4+8):
            p[0] = b[0], p[1] = b[1], p[2] = b[2];
            q[0] = c[0], q[1] = c[1], q[2] = c[2];
            break;
         case (2+8):
            p[0] = b[0], p[1] = b[1], p[2] = b[2];
            q[0] = d[0], q[1] = d[1], q[2] = d[2];
            break;
         case (4+8):
            p[0] = c[0], p[1] = c[1], p[2] = c[2];
            q[0] = d[0], q[1] = d[1], q[2] = d[2];
            break;
         case (1+2+4+8):
            p[0] = a[0], p[1] = a[1], p[2] = a[2];
            q[0] = a[0], q[1] = a[1], q[2] = a[2];
            break;
         default:
            yac_internal_abort_message("internal error", __FILE__, __LINE__);
      }

   } else {

      p[0] = 0, p[1] = 0; p[2] = 1;
      q[0] = 0, q[1] = 0; q[2] = -1;
   }

   double dot_ab = a[0]*b[0] + a[1]*b[1] + a[2]*b[2];
   double dot_cd = c[0]*d[0] + c[1]*d[1] + c[2]*d[2];

   if (ret_value & 16) {
      if (vector_is_between(a, b, p, &angle_ab, dot_ab)) ret_value |= 1;
      if (vector_is_between(a, b, q, &angle_ab, dot_ab)) ret_value |= 2;
      if (vector_is_between(c, d, p, &angle_cd, dot_cd)) ret_value |= 4;
      if (vector_is_between(c, d, q, &angle_cd, dot_cd)) ret_value |= 8;
   } else {
      if (vector_is_between(a, b, p, &angle_ab, dot_ab)) ret_value |= 1;
      else if (vector_is_between(a, b, q, &angle_ab, dot_ab)) ret_value |= 2;
      if (vector_is_between(c, d, p, &angle_cd, dot_cd)) ret_value |= 4;
      else if (vector_is_between(c, d, q, &angle_cd, dot_cd)) ret_value |= 8;
   }

   return ret_value;
}

static
int loncxlonc_ (struct edge edge_a, struct edge edge_b) {

   // if edge goes across pole
   if (fabs(get_angle(edge_a.points[0].lon, edge_a.points[1].lon)) > tol) {

      double lat = (edge_a.points[0].lat > 0)?M_PI_2:-M_PI_2;

      return
         loncxlonc_(
            (struct edge){
               .edge_type = LON_CIRCLE,
               .points = {{.lon = edge_a.points[0].lon,
                           .lat = edge_a.points[0].lat},
                          {.lon = edge_a.points[0].lon,
                           .lat = lat}}}, edge_b)
         ||
         loncxlonc_(
            (struct edge){
               .edge_type = LON_CIRCLE,
               .points = {{.lon = edge_a.points[1].lon,
                           .lat = edge_a.points[1].lat},
                          {.lon = edge_a.points[1].lon,
                           .lat = lat}}}, edge_b);

   } else if (fabs(get_angle(edge_b.points[0].lon,
                             edge_b.points[1].lon)) > tol) {

      double lat = (edge_b.points[0].lat > 0)?M_PI_2:-M_PI_2;

      return
         loncxlonc_(edge_a,
            (struct edge){
               .edge_type = LON_CIRCLE,
               .points = {{.lon = edge_b.points[0].lon,
                           .lat = edge_b.points[0].lat},
                          {.lon = edge_b.points[0].lon,
                           .lat = lat}}})
         ||
         loncxlonc_(edge_a,
            (struct edge){
               .edge_type = LON_CIRCLE,
               .points = {{.lon = edge_b.points[1].lon,
                           .lat = edge_b.points[1].lat},
                          {.lon = edge_b.points[1].lon,
                           .lat = lat}}});
   }

   // if both edges touch a pole
   if (((fabs(fabs(edge_a.points[0].lat) - M_PI_2) < tol) ||
        (fabs(fabs(edge_a.points[1].lat) - M_PI_2) < tol)) &&
       ((fabs(fabs(edge_b.points[0].lat) - M_PI_2) < tol) ||
        (fabs(fabs(edge_b.points[1].lat) - M_PI_2) < tol))) {

      return !((edge_a.points[0].lat > 0) ^ (edge_b.points[0].lat > 0));
   }

   if (fabs(get_angle(edge_a.points[0].lon, edge_b.points[0].lon)) > tol)
      return 0;

   double angle_ab = fabs(edge_a.points[0].lat - edge_a.points[1].lat) + tol;
   double angle_cd = fabs(edge_b.points[0].lat - edge_b.points[1].lat) + tol;
   double angle_ac = fabs(edge_a.points[0].lat - edge_b.points[0].lat);
   double angle_ad = fabs(edge_a.points[0].lat - edge_b.points[1].lat);
   double angle_bc = fabs(edge_a.points[1].lat - edge_b.points[0].lat);
   double angle_bd = fabs(edge_a.points[1].lat - edge_b.points[1].lat);

   return ((angle_ac + angle_bc) < angle_ab) ||
          ((angle_ad + angle_bd) < angle_ab) ||
          ((angle_ac + angle_ad) < angle_cd) ||
          ((angle_bc + angle_bd) < angle_cd);
}

/** \brief compute the intersection point of a meridian and a parallel
 *
 * compute the intersection points of a circle of longitude (defined by a and b)
 * and a circle of latitude (defined by c and d)
 * if p and q are != NULL they contain the intersection points
 * the return value is:
 *      - 0 if the intersection points are neither between (a and b) or (c and d)
 *      - -1 if an error occurred
 *      - 1st bit will be set if p is between a and b
 *      - 2nd bit will be set if q is between a and b
 *      - 3rd bit will be set if p is between c and d
 *      - 4th bit will be set if q is between c and d
 **/
int yac_loncxlatc (struct edge edge_a, struct edge edge_b,
                   struct point * p, struct point * q) {

   double lon_a, lat_a[2];
   double lon_b[2], lat_b;

   unsigned ret_value;

   ret_value = 0;

   lon_a = edge_a.points[0].lon;

   if (edge_a.points[0].lat > edge_a.points[1].lat) {

      lat_a[0] = edge_a.points[1].lat;
      lat_a[1] = edge_a.points[0].lat;
   } else {
      lat_a[0] = edge_a.points[0].lat;
      lat_a[1] = edge_a.points[1].lat;
   }

   lat_b = edge_b.points[0].lat;

   if (edge_b.points[0].lon > edge_b.points[1].lon) {

      lon_b[0] = edge_b.points[1].lon;
      lon_b[1] = edge_b.points[0].lon;
   } else {
      lon_b[0] = edge_b.points[0].lon;
      lon_b[1] = edge_b.points[1].lon;
   }

   unsigned a_goes_across_pole;
   unsigned b_is_on_pole;

   a_goes_across_pole = fabs(get_angle(edge_a.points[0].lon,
                                           edge_a.points[1].lon)) > tol;
   b_is_on_pole = fabs(M_PI_2 - fabs(lat_b)) < tol;

   if (b_is_on_pole) {

      if (p != NULL) *p = edge_b.points[0];
      if (q != NULL) *q = edge_b.points[1];

      if (((a_goes_across_pole) && (fabs(lat_b - lat_a[0]) < M_PI_2)) ||
          ((fabs(lat_b - lat_a[0]) < tol) ||
           (fabs(lat_b - lat_a[1]) < tol))){

         ret_value |= 1;
      }

      ret_value |= 4;

   } else {

      if (p != NULL) p->lon = lon_a, p->lat = lat_b;
      if (q != NULL) q->lon = lon_a + M_PI, q->lat = lat_b;

      double angle_cd, angle_cp, angle_cq;

      angle_cd = get_angle(lon_b[0], lon_b[1]);
      angle_cp = get_angle(lon_b[0], lon_a);
      angle_cq = get_angle(lon_b[0], lon_a + M_PI);

      if (angle_cd > 0) {
         if (angle_cp <= angle_cd + tol && angle_cp > -tol) ret_value |= 1 << 2;
         if (angle_cq <= angle_cd + tol && angle_cq > -tol) ret_value |= 1 << 3;
      } else {
         if (angle_cp >= angle_cd - tol && angle_cp < tol) ret_value |= 1 << 2;
         if (angle_cq >= angle_cd - tol && angle_cq < tol) ret_value |= 1 << 3;
      }

      if (a_goes_across_pole) {

         if (lat_a[0] > 0.0) {

            if (lat_b > lat_a[0]) ret_value |= 1 << 0;
            if (lat_b > lat_a[1]) ret_value |= 1 << 1;

         } else {

            if (lat_b < lat_a[0]) ret_value |= 1 << 0;
            if (lat_b < lat_a[1]) ret_value |= 1 << 1;
         }
      } else if ((lat_b >= lat_a[0]) && (lat_b <= lat_a[1])) ret_value |= 1 << 0;
   }

   return ret_value;
}

/** \brief compute the intersection point of a meridian and a parallel
 *
 * compute the intersection points of a circle of longitude (defined by a and b)
 * and a circle of latitude (defined by c and d)
 * if p and q are != NULL they contain the intersection points
 * the return value is:
 *      - 0 if the intersection points are neither between (a and b) or (c and d)
 *      - -1 if an error occurred
 *      - 1st bit will be set if p is between a and b
 *      - 2nd bit will be set if q is between a and b
 *      - 3rd bit will be set if p is between c and d
 *      - 4th bit will be set if q is between c and d
 **/
int yac_loncxlatc_vec (double a[3], double b[3], double c[3], double d[3],
                       double p[3], double q[3]) {

   unsigned ret_value;

   ret_value = 0;

   // this test is not very accurate but should catch the most obvious cases
   // the accuracy of this test is not allowed to be higher than the one in then
   // routine is_inside_gc
   if ((fabs(a[0] * b[1] - a[1] * b[0]) > 1e-7) &&
       (fabs(fabs(a[2]) - 1.0) > tol) &&
       (fabs(fabs(b[2]) - 1.0) > tol)) {

      yac_internal_abort_message("edge is not a circle of longitude", __FILE__, __LINE__);
   }

   unsigned ab_goes_across_pole;
   unsigned cd_is_on_pole;
   unsigned ab_is_point;
   unsigned cd_is_point;
   double angle_ab = get_vector_angle(a, b);
   double angle_cd = get_vector_angle(c, d);

   ab_goes_across_pole =
      ((fabs(1.0 - fabs(a[2])) < tol || fabs(1.0 - fabs(b[2])) < tol) ||
       (((a[0] > 0.0) ^ (b[0] > 0.0)) && ((a[1] > 0.0) ^ (b[1] > 0.0))));

   cd_is_on_pole = fabs(1.0 - fabs(c[2])) < tol;

   ab_is_point = angle_ab < tol;

   cd_is_point = angle_cd < tol;

   if (cd_is_on_pole) {

      if (((ab_goes_across_pole) && (fabs(a[2] - c[2]) < 1.0)) ||
          ((fabs(a[2] - c[2]) < tol)) || (fabs(b[2] - c[2]) < tol))
         ret_value |= 1;

         if (p != NULL)
            p[0] = c[0], p[1] = c[1], p[2] = c[2];
         if (q != NULL)
            q[0] = c[0], q[1] = c[1], q[2] = c[2];

      ret_value |= 4;

   } else if (ab_goes_across_pole && ab_is_point) {

      if (p != NULL)
         p[0] = c[0], p[1] = c[1], p[2] = c[2];
      if (q != NULL)
         q[0] = -p[0], q[1] = -p[1], q[2] = p[2];

      ret_value |= 4;

   } else {

      /*
      // the cos is too inaccurate close to the equator
      {
         if (fabs(a[2]) < fabs(b[2])) {

            double scale = cos(c[2]) / cos(a[2]);

            if (p != NULL)
               p[0] = a[0] * scale, p[1] = a[1] * scale, p[2] = c[2];
            if (q != NULL)
               q[0] = -a[0] * scale, q[1] = -a[1] * scale, q[2] = c[2];
         } else {d

            double scale = cos(c[2]) / cos(b[2]);

            if (p != NULL)
               p[0] = b[0] * scale, p[1] = b[1] * scale, p[2] = c[2];
            if (q != NULL)
               q[0] = -b[0] * scale, q[1] = -b[1] * scale, q[2] = c[2];
         }
      }
      */
      {

        double tmp_scale_a = a[0] * a[0] + a[1] * a[1];
        double tmp_scale_b = b[0] * b[0] + b[1] * b[1];

         if (tmp_scale_a > tmp_scale_b) {

            double scale = sqrt((1.0 - c[2] * c[2])/tmp_scale_a);

            if (p != NULL)
               p[0] = a[0] * scale, p[1] = a[1] * scale, p[2] = c[2];
            if (q != NULL)
               q[0] = -a[0] * scale, q[1] = -a[1] * scale, q[2] = c[2];

         } else {

            double scale = sqrt((1.0 - c[2] * c[2])/tmp_scale_b);

            if (p != NULL)
               p[0] = b[0] * scale, p[1] = b[1] * scale, p[2] = c[2];
            if (q != NULL)
               q[0] = -b[0] * scale, q[1] = -b[1] * scale, q[2] = c[2];
         }
      }

      

      if (cd_is_point) {

         if (fabs(c[0]-p[0]) < tol &&
             fabs(c[1]-p[1]) < tol &&
             fabs(c[2]-p[2]) < tol) ret_value |= 1 << 2;
         if (fabs(c[0]-q[0]) < tol &&
             fabs(c[1]-q[1]) < tol &&
             fabs(c[2]-q[2]) < tol) ret_value |= 1 << 3;

      } else {

         if (vector_is_between_lat(c, d, p))
            ret_value |= 1 << 2;
         if (vector_is_between_lat(c, d, q))
            ret_value |= 1 << 3;
      }

      if (ab_is_point) {
         if (fabs(a[0]-p[0]) < tol &&
             fabs(a[1]-p[1]) < tol &&
             fabs(a[2]-p[2]) < tol) ret_value |= 1 << 0;
         if (fabs(a[0]-q[0]) < tol &&
             fabs(a[1]-q[1]) < tol &&
             fabs(a[2]-q[2]) < tol) ret_value |= 1 << 1;
      } else {

         double dot_ab = a[0]*b[0] + a[1]*b[1] + a[2]*b[2];

         if (vector_is_between(a, b, p, &angle_ab, dot_ab)) ret_value |= 1 << 0;
         if (vector_is_between(a, b, q, &angle_ab, dot_ab)) ret_value |= 1 << 1;
      }
   }

   return ret_value;
}

static int loncxlatc_ (struct edge edge_a, struct edge edge_b) {

   // if edge goes across pole
   if (fabs(get_angle(edge_a.points[0].lon, edge_a.points[1].lon)) > tol) {

      double lat = (edge_a.points[0].lat > 0)?M_PI_2:-M_PI_2;

      return
         loncxlatc_(
            (struct edge){
               .edge_type = LON_CIRCLE,
               .points = {{.lon = edge_a.points[0].lon,
                           .lat = edge_a.points[0].lat},
                          {.lon = edge_a.points[0].lon,
                           .lat = lat}}}, edge_b)
         ||
         loncxlatc_(
            (struct edge){
               .edge_type = LON_CIRCLE,
               .points = {{.lon = edge_a.points[1].lon,
                           .lat = edge_a.points[1].lat},
                          {.lon = edge_a.points[1].lon,
                           .lat = lat}}}, edge_b);
   }

   // if edge b is at the pole
   if (fabs(fabs(edge_b.points[0].lat) - M_PI_2) < tol) {

      return (fabs(edge_a.points[0].lat - edge_b.points[0].lat) < tol) ||
             (fabs(edge_a.points[1].lat - edge_b.points[0].lat) < tol);
   }

   double lat_diff = fabs(edge_a.points[0].lat - edge_a.points[1].lat) + tol;

   if ((fabs(edge_a.points[0].lat - edge_b.points[0].lat) > lat_diff) ||
       (fabs(edge_a.points[1].lat - edge_b.points[0].lat) > lat_diff))
      return 0;

   double lon_diff = fabs(get_angle(edge_b.points[0].lon,
                                    edge_b.points[1].lon)) + tol;

   return (fabs(edge_a.points[0].lon - edge_b.points[0].lon) < lon_diff) &&
          (fabs(edge_a.points[0].lon - edge_b.points[1].lon) < lon_diff);
}

/** \brief compute the intersection of a great circle with the parallel
 *
 *  compute the intersection points of a great circle (defined by a and b)
 * and a circle of latitude (defined by c and d)
 * if p and q are != NULL they contain the intersection points
 * the return value is:
 *    - 0 if the intersection points are neither between (a and b) or (c and d)
 *    - -1 if the two circles do not intersect or an error occurred
 *    - 1st bit will be set if p is between a and b
 *    - 2nd bit will be set if q is between a and b
 *    - 3rd bit will be set if p is between c and d
 *    - 4th bit will be set if q is between c and d
 *    - 5th bit will be set if both circles are identically
 *
 *   based on
 *   - http://geospatialmethods.org/spheres/GCIntersect.html
 **/
int yac_gcxlatc(struct edge edge_a, struct edge edge_b,
                struct point * p, struct point * q) {

   // if the great circle is nearly a lon circle, then the accuracy of the normal
   // computation gets messy, therefore we handle it as a lon circle...

   if (fabs(edge_a.points[0].lon - edge_a.points[1].lon) < tol) {

      edge_a.points[0].lon = edge_a.points[1].lon = (edge_a.points[0].lon + edge_a.points[1].lon) / 2.0;
      edge_a.edge_type = LON_CIRCLE;

      return yac_loncxlatc(edge_a, edge_b, p, q);
   }

   // if the great circle is the equator, we handle the great circle as a circle of latitude
   if (fabs(edge_a.points[0].lat) < tol && fabs(edge_a.points[1].lat) < tol) {

      edge_a.points[0].lat = edge_a.points[1].lat = 0;
      edge_a.edge_type = LAT_CIRCLE;
      return yac_latcxlatc(edge_a, edge_b, p, q);
   }

   double a[3], b[3], c[3], d[3], p_[3], q_[3];

   LLtoXYZ(edge_a.points[0].lon, edge_a.points[0].lat, a);
   LLtoXYZ(edge_a.points[1].lon, edge_a.points[1].lat, b);
   LLtoXYZ(edge_b.points[0].lon, edge_b.points[0].lat, c);
   LLtoXYZ(edge_b.points[1].lon, edge_b.points[1].lat, d);

   int ret_value = yac_gcxlatc_vec(a, b, c, d, p_, q_);

   if (ret_value == -1) return -1;

   if (p != NULL) XYZtoLL(p_, &p->lon, &p->lat);
   if (q != NULL) XYZtoLL(q_, &q->lon, &q->lat);

   return ret_value;
}

/** \brief compute the intersection of a great circle with the parallel
 *
 *  compute the intersection points of a great circle (defined by a and b)
 * and a circle of latitude (defined by c and d)
 * if p and q are != NULL they contain the intersection points
 * the return value is:
 *    - 0 if the intersection points are neither between (a and b) or (c and d)
 *    - -1 if the two circles do not intersect or an error occurred
 *    - 1st bit will be set if p is between a and b
 *    - 2nd bit will be set if q is between a and b
 *    - 3rd bit will be set if p is between c and d
 *    - 4th bit will be set if q is between c and d
 *    - 5th bit will be set if both circles are identically
 * \remarks if -1 is returned neither p or q is set
 * \remarks if the two circles only have one intersection point,
 *          p and q will be identically, but only the p bits will be set
 **/

int yac_gcxlatc_vec(double a[3], double b[3], double c[3], double d[3],
                    double p[3], double q[3]) {

   unsigned result = 0;

   double cross_ab[3], scale;
   
   crossproduct_ld(a, b, cross_ab);
   scale = sqrt(cross_ab[0]*cross_ab[0]+
                cross_ab[1]*cross_ab[1]+
                cross_ab[2]*cross_ab[2]);

   // if the great circle is the equator
   if (fabs(a[2]) < tol && fabs(b[2]) < tol) {

      return yac_latcxlatc_vec(a, b, c, d, p, q);

   // if the great circle is  a circle of longitude
   } else if (scale < tol || fabs(cross_ab[2]/scale) < tol ||
              fabs(fabs(a[2])-1.0) < 1e-13 ||
              fabs(fabs(b[2])-1.0) < 1e-13) {

      return yac_loncxlatc_vec(a, b, c, d, p, q);
   }

   double t[3], s[3];

   if (fabs(a[2]) > fabs(b[2])) {

      double scale = c[2] / a[2];

      t[0] = scale * a[0];
      t[1] = scale * a[1];
      
   } else {

      double scale = c[2] / b[2];

      t[0] = scale * b[0];
      t[1] = scale * b[1];
   }

   t[2] = c[2];

   s[2] = 0;

   if (fabs(a[2]) < tol)
      s[0] = a[0], s[1] = a[1];
   else if (fabs(b[2]) < tol)
      s[0] = b[0], s[1] = b[1];
   else if (fabs(a[2]) > fabs(b[2])) {
      double scale = b[2] / a[2];
      s[0] = b[0] - scale * a[0];
      s[1] = b[1] - scale * a[1];
   } else {
      double scale = a[2] / b[2];
      s[0] = a[0] - scale * b[0];
      s[1] = a[1] - scale * b[1];
   }

   {
      // the intersection of the planes of both circles is defined by:
      // x = t + n * s

      // x_0^2 + x_1^2 + x_2^2 = 1
      // x_2 = c_2

      double a_ = s[0] * s[0] + s[1] * s[1];
      double b_ = 2.0 * (t[0] * s[0] + t[1] * s[1]);
      double c_ = t[0] * t[0] + t[1] * t[1] + c[2] * c[2] - 1.0;

      if (a_ == 0.0)
         yac_internal_abort_message("internal error", __FILE__, __LINE__);

      double temp = b_ * b_ - 4.0 * a_ * c_;

      // no intersection possible
      if (temp < 0.0) {
        if (temp < -tol)
          return -1;
        else
          temp = 0;
      }

      double n[2];

      n[0] = - (b_ + sqrt(temp)) / (2.0 * a_);
      n[1] = - (b_ - sqrt(temp)) / (2.0 * a_);

      p[0] = t[0] + n[0] * s[0];
      p[1] = t[1] + n[0] * s[1];
      p[2] = t[2] + n[0] * s[2];

      double angle_ab = -1;
      double dot_ab = a[0]*b[0] + a[1]*b[1] + a[2]*b[2];

      if (vector_is_between(a, b, p, &angle_ab, dot_ab)) result |= 1;
      if (vector_is_between_lat(c, d, p)) result |= 4;

      if (fabs(n[0] - n[1]) >= tol) {

         q[0] = t[0] + n[1] * s[0];
         q[1] = t[1] + n[1] * s[1];
         q[2] = t[2] + n[1] * s[2];

         if (vector_is_between(a, b, q, &angle_ab, dot_ab)) result |= 2;
         if (vector_is_between_lat(c, d, q)) result |= 8;
      } else
         q[0] = p[0], q[1] = p[1], q[2] = p[2];
   }

   return result;
}

static
int gcxlatc_vec_(double a[3], double b[3], double c[3], double d[3]) {

   double angle_ab = get_vector_angle(a, b);
   double angle_cd = get_vector_angle(c, d);

   // if ab is a point
   if (angle_ab < tol) {

      return (fabs(a[2] - c[2]) < tol) && vector_is_between_lat(c, d, a);

   } else if (angle_cd < tol) {

      double cross_ab[3];

      crossproduct_ld(a, b, cross_ab);

      double n = 1.0 / sqrt(cross_ab[0] * cross_ab[0] +
                            cross_ab[1] * cross_ab[1] +
                            cross_ab[2] * cross_ab[2]);
      cross_ab[0] *= n;
      cross_ab[1] *= n;
      cross_ab[2] *= n;

      double dot_ab = a[0]*b[0] + a[1]*b[1] + a[2]*b[2];

      return (fabs(get_vector_angle(cross_ab, c) - M_PI_2) < tol) &&
             vector_is_between(a, b, c, &angle_ab, dot_ab);
   }

   double t[3], s[3];

   if (fabs(a[2]) > fabs(b[2])) {

      double scale = c[2] / a[2];

      t[0] = scale * a[0];
      t[1] = scale * a[1];
      
   } else {

      double scale = c[2] / b[2];

      t[0] = scale * b[0];
      t[1] = scale * b[1];
   }

   t[2] = c[2];

   s[2] = 0;

   if (fabs(a[2]) < tol)
      s[0] = a[0], s[1] = a[1];
   else if (fabs(b[2]) < tol)
      s[0] = b[0], s[1] = b[1];
   else if (fabs(a[2]) > fabs(b[2])) {
      double scale = b[2] / a[2];
      s[0] = b[0] - scale * a[0];
      s[1] = b[1] - scale * a[1];
   } else {
      double scale = a[2] / b[2];
      s[0] = a[0] - scale * b[0];
      s[1] = a[1] - scale * b[1];
   }

   if (fabs(s[0]) < tol && fabs(s[1]) < tol)
      yac_internal_abort_message("internal error", __FILE__, __LINE__);

   {
      // the intersection of the planes of both circles is defined by:
      // x = t + n * s

      // x_0^2 + x_1^2 + x_2^2 = 1
      // x_2 = c_2

      double a_ = s[0] * s[0] + s[1] * s[1];
      double b_ = 2.0 * (t[0] * s[0] + t[1] * s[1]);
      double c_ = t[0] * t[0] + t[1] * t[1] + c[2] * c[2] - 1.0;

      double temp = b_ * b_ - 4.0 * a_ * c_;

      // no intersection possible
      if (temp < 0.0) {
        if (temp < -tol)
          return 0;
        else
          temp = 0;
      }

      double n[2];

      n[0] = - (b_ + sqrt(temp)) / (2.0 * a_);
      n[1] = - (b_ - sqrt(temp)) / (2.0 * a_);

      double p[3];

      p[0] = t[0] + n[0] * s[0];
      p[1] = t[1] + n[0] * s[1];
      p[2] = t[2] + n[0] * s[2];

      double dot_ab = a[0]*b[0] + a[1]*b[1] + a[2]*b[2];

      if (vector_is_between(a, b, p, &angle_ab, dot_ab) &&
          vector_is_between_lat(c, d, p))
         return 1;

      if (fabs(n[0] - n[1]) >= tol) {

         p[0] = t[0] + n[1] * s[0];
         p[1] = t[1] + n[1] * s[1];
         p[2] = t[2] + n[1] * s[2];

         return vector_is_between(a, b, p, &angle_ab, dot_ab) &&
                vector_is_between_lat(c, d, p);
      }  else
         return 0;
   }
}

int yac_intersect (struct edge const edge_a, struct edge const edge_b,
                   struct point * intersection) {

   int switch_edges;

   switch_edges = 0;

   int (*intersect_func)(struct edge, struct edge, struct point *, struct point *);

   // if both edges are on circles of latitude
   if (edge_a.edge_type == LAT_CIRCLE && edge_b.edge_type == LAT_CIRCLE) {

      intersect_func = yac_latcxlatc;

   // if both edges are on circle of longitude
   } else if (edge_a.edge_type == LON_CIRCLE && edge_b.edge_type == LON_CIRCLE) {

      intersect_func = yac_loncxlonc;

   // if both edges are on great circles
   } else if ((edge_a.edge_type == GREAT_CIRCLE &&
        edge_b.edge_type == GREAT_CIRCLE) ||
       (edge_a.edge_type == LON_CIRCLE   &&
        edge_b.edge_type == GREAT_CIRCLE) ||
       (edge_a.edge_type == GREAT_CIRCLE &&
        edge_b.edge_type == LON_CIRCLE)) {

      intersect_func = yac_gcxgc;

   // if one edge a is on a great circle and edge b on a circle of latitude
   } else if (edge_a.edge_type == GREAT_CIRCLE &&
              edge_b.edge_type == LAT_CIRCLE) {

      intersect_func = yac_gcxlatc;

   // if one edge a is on a circle of latitude and edge b on a great circle
   } else if (edge_a.edge_type == LAT_CIRCLE &&
              edge_b.edge_type == GREAT_CIRCLE ) {

      switch_edges = 1;
      intersect_func = yac_gcxlatc;

   // if one edge a is on a circle of longitude and edge b on a circle of latitude
   } else if (edge_a.edge_type == LON_CIRCLE &&
              edge_b.edge_type == LAT_CIRCLE) {

      intersect_func = yac_loncxlatc;

   // if one edge a is on a circle of latitude and edge b on a circle of longitude
   } else if (edge_a.edge_type == LAT_CIRCLE &&
              edge_b.edge_type == LON_CIRCLE ) {

      switch_edges = 1;
      intersect_func = yac_loncxlatc;

   } else {

      yac_internal_abort_message ( "ERROR: unknown edge type.", __FILE__, __LINE__ );
      exit(EXIT_FAILURE);
   }

   unsigned const p_on_a = 1 << 0;
   unsigned const q_on_a = 1 << 1;
   unsigned const p_on_b = 1 << 2;
   unsigned const q_on_b = 1 << 3;

   int ret_value;

   struct point p, q;

   // compute the intersection between both circles
   if (switch_edges) ret_value = intersect_func(edge_b, edge_a, &p, &q);
   else              ret_value = intersect_func(edge_a, edge_b, &p, &q);

   // check whether the circles defined by the edges intersect
   if (ret_value > 0) {

      // check for intersection of the edges
      if ( (ret_value & p_on_a) && (ret_value & p_on_b) ) {

         if (intersection != NULL) *intersection = p;
         return 1;

      } else if ( (ret_value & q_on_a) && (ret_value & q_on_b) ) {

         if (intersection != NULL) *intersection = q;
         return 1;
      }
   }
   return 0;
}

int yac_intersect_vec (enum yac_edge_type edge_type_a, double a[3], double b[3],
                       enum yac_edge_type edge_type_b, double c[3], double d[3],
                       double p[3], double q[3]) {

   int switch_edges;

   switch_edges = 0;

   int (*intersect_func)(double *, double *, double *, double *, double *, double *);

   // if both edges are on circles of latitude
   if (edge_type_a == LAT_CIRCLE &&
       edge_type_b == LAT_CIRCLE) {

      intersect_func = yac_latcxlatc_vec;

   // if both edges are on circle of longitude
   } else if (edge_type_a == LON_CIRCLE &&
              edge_type_b == LON_CIRCLE) {

      intersect_func = yac_loncxlonc_vec;

   // if both edges are on great circles
   } else if ((edge_type_a == GREAT_CIRCLE &&
               edge_type_b == GREAT_CIRCLE) ||
              (edge_type_a == LON_CIRCLE   &&
               edge_type_b == GREAT_CIRCLE) ||
              (edge_type_a == GREAT_CIRCLE &&
               edge_type_b == LON_CIRCLE)) {

      intersect_func = yac_gcxgc_vec;

   // if one edge a is on a great circle and edge b on a circle of latitude
   } else if (edge_type_a == GREAT_CIRCLE &&
              edge_type_b == LAT_CIRCLE) {

      intersect_func = yac_gcxlatc_vec;

   // if one edge a is on a circle of latitude and edge b on a great circle
   } else if (edge_type_a == LAT_CIRCLE &&
              edge_type_b == GREAT_CIRCLE ) {

      switch_edges = 1;
      intersect_func = yac_gcxlatc_vec;

   // if one edge a is on a circle of longitude and edge b on a circle of latitude
   } else if (edge_type_a == LON_CIRCLE &&
              edge_type_b == LAT_CIRCLE) {

      intersect_func = yac_loncxlatc_vec;

   // if one edge a is on a circle of latitude and edge b on a circle of longitude
   } else if (edge_type_a == LAT_CIRCLE &&
              edge_type_b == LON_CIRCLE ) {

      switch_edges = 1;
      intersect_func = yac_loncxlatc_vec;

   } else {

      yac_internal_abort_message ( "ERROR: unknown edge type.", __FILE__, __LINE__ );
      exit(EXIT_FAILURE);
   }

   int ret_value;

   // compute the intersection between both circles
   if (switch_edges) ret_value = intersect_func(c, d, a, b, p, q);
   else              ret_value = intersect_func(a, b, c, d, p, q);

   if (switch_edges)
      ret_value = (ret_value & (~(1 + 2 + 4 + 8))) +
                  ((ret_value & (1 + 2)) << 2) +
                  ((ret_value & (4 + 8)) >> 2);

   return ret_value;
}

int yac_do_intersect (struct edge edge_a, double a[3], double b[3],
                      struct edge edge_b, double c[3], double d[3]) {

   int flag = ((edge_a.edge_type == LAT_CIRCLE)   << 0) |
              ((edge_a.edge_type == LON_CIRCLE)   << 1) |
              ((edge_a.edge_type == GREAT_CIRCLE) << 2) |
              ((edge_b.edge_type == LAT_CIRCLE)   << 3) |
              ((edge_b.edge_type == LON_CIRCLE)   << 4) |
              ((edge_b.edge_type == GREAT_CIRCLE) << 5);

   switch (flag) {

      case ((1 << 0) | (1 << 3)):
         return latcxlatc_(edge_a, edge_b);
      case ((1 << 0) | (1 << 4)):
         return loncxlatc_(edge_b, edge_a);
      case ((1 << 0) | (1 << 5)):
         return gcxlatc_vec_(c, d, a, b);
      case ((1 << 1) | (1 << 3)):
         return loncxlatc_(edge_a, edge_b);
      case ((1 << 1) | (1 << 4)):
         return loncxlonc_(edge_a, edge_b);
      case ((1 << 1) | (1 << 5)):
      case ((1 << 2) | (1 << 4)):
      case ((1 << 2) | (1 << 5)):
         return gcxgc_vec_(a, b, c, d);
      case ((1 << 2) | (1 << 3)):
         return gcxlatc_vec_(a, b, c, d);
      default:
         yac_internal_abort_message ( "ERROR: unknown edge type.", __FILE__, __LINE__ );
         exit(EXIT_FAILURE);
   };
}
