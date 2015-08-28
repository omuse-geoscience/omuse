#include <stdio.h>
#include <math.h>
#include "grid.h"

#ifndef IS_EQUAL
#  define IS_NOT_EQUAL(x,y) (x < y || y < x)
#  define IS_EQUAL(x,y)     (!IS_NOT_EQUAL(x,y))
#endif

#ifndef  M_PI
#define  M_PI		3.14159265358979323846	/* pi */
#endif

const double pi = M_PI;
#define deg_per_rad  (180./pi)
#define rad_per_deg  (pi / 180.)

// Mean Earth Radius in m.  The value below is consistent with NCEP's routines and grids.
const double earth_radius_m = 6371200.;
// const double earth_radius_m = 6370000.; //  consistent with rest of MM5 system
#define radians_per_degree  (pi / 180.)
 

static
void map_init(proj_info_t *proj)
{
  // Initializes the map projection structure to missing values

  proj->lat1     = -999.9;
  proj->lon1     = -999.9;
  proj->dx       = -999.9;
  proj->stdlon   = -999.9;
  proj->truelat1 = -999.9;
  proj->truelat2 = -999.9;
  proj->cone     = -999.9;
  proj->polei    = -999.9;
  proj->polej    = -999.9;
  proj->rsw      = -999.9;
  proj->hemi     = 0;
  proj->init     = 0;
}

static
double lc_cone(double truelat1, double truelat2)
{
  // Subroutine to compute the cone factor of a Lambert Conformal projection
    
  // Input Args
  //    truelat1   (-90 -> 90 degrees N)
  //    truelat2     "   "  "   "     "

  // First, see if this is a secant or tangent projection.  For tangent
  // projections, truelat1 = truelat2 and the cone is tangent to the 
  // Earth's surface at this latitude.  For secant projections, the cone
  // intersects the Earth's surface at each of the distinctly different
  // latitudes

  double cone;

  if (fabs(truelat1-truelat2) > 0.1)
    {
      cone = log10(cos(truelat1*rad_per_deg)) - 
  	     log10(cos(truelat2*rad_per_deg));
      cone = cone /(log10(tan((45.0 - fabs(truelat1)/2.0) * rad_per_deg)) -
		    log10(tan((45.0 - fabs(truelat2)/2.0) * rad_per_deg)));   
    }
  else
    {
      cone = sin(fabs(truelat1)*rad_per_deg);
    }

  return cone;
}

static
void set_lc(proj_info_t *proj)
{
  // Initialize the remaining items in the proj structure for a
  // lambert conformal grid.

  double  arg;
  double  deltalon1;
  double  tl1r;
  double  ctl1r;

  // Compute cone factor
  proj->cone = lc_cone(proj->truelat1, proj->truelat2);
  // fprintf(stdout, "Computed cone factor: %g\n", proj->cone);

  // Compute longitude differences and ensure we stay out of the
  // forbidden "cut zone"
  deltalon1 = proj->lon1 - proj->stdlon;
  if (deltalon1 > +180.) deltalon1 = deltalon1 - 360.;
  if (deltalon1 < -180.) deltalon1 = deltalon1 + 360.;

  // Convert truelat1 to radian and compute COS for later use
  tl1r = proj->truelat1 * rad_per_deg;
  ctl1r = cos(tl1r);

  // Compute the radius to our known lower-left (SW) corner
  proj->rsw = proj->rebydx * ctl1r/proj->cone *
    pow((tan((90.*proj->hemi-proj->lat1)*rad_per_deg/2.) /
	 tan((90.*proj->hemi-proj->truelat1)*rad_per_deg/2.)), proj->cone);

  // Find pole point
  arg = proj->cone*(deltalon1*rad_per_deg);
  proj->polei = 1. - proj->hemi * proj->rsw * sin(arg);
  proj->polej = 1. + proj->rsw * cos(arg)  ;
  // fprintf(stdout, "Computed pole (x,y) = %g %g\n", proj->polei, proj->polej);
}


void map_set(int proj_code, double lat1, double lon1, double dx, double stdlon,
	     double truelat1, double truelat2, proj_info_t *proj)
{
  // Given a partially filled proj_info structure, this routine computes
  // polei, polej, rsw, and cone (if LC projection) to complete the 
  // structure.  This allows us to eliminate redundant calculations when
  // calling the coordinate conversion routines multiple times for the
  // same map.
  // This will generally be the first routine called when a user wants
  // to be able to use the coordinate conversion routines, and it
  // will call the appropriate subroutines based on the 
  // proj->code which indicates which projection type  this is.
    
  // First, check for validity of mandatory variables in proj
  /*
  if ( fabs(lat1) > 90. ) {
    PRINT '(A)', 'Latitude of origin corner required as follows:'
      PRINT '(A)', '    -90N <= lat1 < = 90.N'
      STOP 'MAP_INIT'
    }
    if ( fabs(lon1) > 180.) {
      PRINT '(A)', 'Longitude of origin required as follows:'
      PRINT '(A)', '   -180E <= lon1 <= 180W'
      STOP 'MAP_INIT'
    }
    if ((dx <= 0.).AND.(proj_code .NE. PROJ_LATLON)) {
      PRINT '(A)', 'Require grid spacing (dx) in meters be positive//'
      STOP 'MAP_INIT'
    }
    if ((fabs(stdlon) > 180.).AND.(proj_code != PROJ_MERC)) {
      PRINT '(A)', 'Need orientation longitude (stdlon) as: '
      PRINT '(A)', '   -180E <= lon1 <= 180W' 
      STOP 'MAP_INIT'
    }
    if (fabs(truelat1)>90.) {
      PRINT '(A)', 'Set true latitude 1 for all projections//'
      STOP 'MAP_INIT'
    }
  */
  map_init(proj);

  proj->code     = proj_code;
  proj->lat1     = lat1;
  proj->lon1     = lon1;
  proj->dx       = dx;
  proj->stdlon   = stdlon;
  proj->truelat1 = truelat1;
  proj->truelat2 = truelat2;

  if ( proj->code != PROJ_LATLON )
    {
      proj->dx = dx;
      if (truelat1 < 0.)
	proj->hemi = -1;
      else
	proj->hemi =  1;

      proj->rebydx = earth_radius_m / dx;
    }
  /*
 pick_proj: SELECT CASE(proj->code)

      CASE(PROJ_PS)
        PRINT '(A)', 'Setting up POLAR STEREOGRAPHIC map...'
        CALL set_ps(proj)

      CASE(PROJ_LC)
  */
  // fprintf(stdout, "Setting up LAMBERT CONFORMAL map...\n");
  if (fabs(proj->truelat2) > 90.)
    {
      // fprintf(stdout, "Second true latitude not set, assuming a tangent\n");
      // fprintf(stdout, "projection at truelat1: %g\n", proj->truelat1);
      proj->truelat2 = proj->truelat1;
    }

  set_lc(proj);
  /*
      CASE (PROJ_MERC)
        PRINT '(A)', 'Setting up MERCATOR map...'
        CALL set_merc(proj)
   
      CASE (PROJ_LATLON)
        PRINT '(A)', 'Setting up CYLINDRICAL EQUIDISTANT LATLON map...'
        // Convert lon1 to 0->360 notation
        if (proj->lon1 < 0.) proj->lon1 = proj->lon1 + 360.
   
      CASE DEFAULT
        PRINT '(A,I2)', 'Unknown projection code: ', proj->code
        STOP 'MAP_INIT'
    
    END SELECT pick_proj
  */

  proj->init = 1;
}


void ijll_lc(double i, double j, proj_info_t proj, double *lat, double *lon)
{
  // Subroutine to convert from the (i,j) cartesian coordinate to the 
  // geographical latitude and longitude for a Lambert Conformal projection.

  // History:
  // 25 Jul 01: Corrected by B. Shaw, NOAA/FSL
  // 10 Sep 09: Converted to ANSI C: Uwe Schulzweida, MPIMET

  // Input Args
  // double, INTENT(IN)              :: i        // Cartesian X coordinate
  // double, INTENT(IN)              :: j        // Cartesian Y coordinate
  // TYPE(proj_info),INTENT(IN)      :: proj     // Projection info structure

  // Output Args                 
  // double, INTENT(OUT)             :: lat      // Latitude (-90->90 deg N)
  // double, INTENT(OUT)             :: lon      // Longitude (-180->180 E)

  // Locals 
  double  inew;
  double  jnew;
  double  r;
  double  chi,chi1,chi2;
  double  r2;
  double  xx;
  double  yy;

  chi1 = (90. - proj.hemi*proj.truelat1)*rad_per_deg;
  chi2 = (90. - proj.hemi*proj.truelat2)*rad_per_deg;
    
  // See if we are in the southern hemispere and flip the indices if we are. 
  if ( proj.hemi == -1 )
    { 
      inew = -i + 2.;
      jnew = -j + 2.;
    }
  else
    {
      inew = i;
      jnew = j;
    }

  // Compute radius**2 to i/j location
  xx = inew - proj.polei;
  yy = proj.polej - jnew;
  r2 = (xx*xx + yy*yy);
  r  = sqrt(r2)/proj.rebydx;
   
  // Convert to lat/lon
  if ( IS_EQUAL(r2, 0.) )
    {
      *lat = proj.hemi * 90.;
      *lon = proj.stdlon;
    }
  else
    {
      // Longitude
      *lon = proj.stdlon + deg_per_rad * atan2(proj.hemi*xx,yy)/proj.cone;
      *lon = fmod(*lon+360., 360.);

      // Latitude.  Latitude determined by solving an equation adapted 
      // from:
      //  Maling, D.H., 1973: Coordinate Systems and Map Projections
      // Equations #20 in Appendix I.  
        
      if ( IS_EQUAL(chi1, chi2) )
	chi = 2.0*atan( pow( r/tan(chi1), (1./proj.cone) ) * tan(chi1*0.5) );
      else
	chi = 2.0*atan( pow( r*proj.cone/sin(chi1), (1./proj.cone) ) * tan(chi1*0.5)) ;

      *lat = (90.0-chi*deg_per_rad)*proj.hemi;
    }

  if ( *lon > +180. ) *lon = *lon - 360.;
  if ( *lon < -180. ) *lon = *lon + 360.;
}

/*
int main(void)
{
  int    status = 0;
  int    nlon = 245;
  int    nlat = 277;
  double xi, xj;
  double lat_ll_p   =  47.806;
  double lon_ll_p   = -10.063;
  double lat_tan_p  =  59.2;
  double dx_p       =  11000.0;
  double lon_xx_p   = -10.0;
  double zlat, zlon;
  proj_info_t proj;

  map_set(PROJ_LC, lat_ll_p, 360.+lon_ll_p, dx_p, 360.+lon_xx_p, lat_tan_p, lat_tan_p, &proj);

  xi = 1;
  xj = 1;
  ijll_lc(xi, xj, proj, &zlat, &zlon);
  printf("1 1 47.806 349.937 0\n");
  printf("%g %g %g %g %d\n", xj, xi, zlat, zlon, status);

  xi = nlon;
  xj = nlat;
  ijll_lc(xi, xj, proj, &zlat, &zlon);
  printf("277 245 63.086 51.4192 0\n");
  printf("%g %g %g %g %d\n", xj, xi, zlat, zlon, status);
			  

  {
    int i, j;
    for ( j = 1; j <= nlat; j++ )
      for ( i = 1; i <= nlon; i++ )
	{
	  xi = i;
	  xj = j;
	  ijll_lc(xi, xj, proj, &zlat, &zlon);
	  if ( zlon < 0 ) zlon+=360;
	}
  }
			
  return (0);
}
*/
