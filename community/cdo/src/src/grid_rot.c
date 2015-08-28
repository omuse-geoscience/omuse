#include <stdio.h>
#include <math.h>

#include "grid.h"


double lamrot_to_lam(double phirot, double lamrot, double polphi, double pollam, double polgam)
{
  /*
    This function converts lambda from one rotated system to lambda in another system. 
    If the optional argument polgam is present, the other system can also be a rotated one, 
    where polgam is the angle between the two north poles.
    If polgam is not present, the other system is the real geographical system.

    phirot : latitude in the rotated system
    lamrot : longitude in the rotated system (E>0)
    polphi : latitude of the rotated north pole
    pollam : longitude of the rotated north pole

    result : longitude in the geographical system
  */
  double zsinpol, zcospol, zlampol;
  double zphirot, zlamrot, zarg1, zarg2;
  double zgam;
  double result = 0;

  zsinpol = sin(DEG2RAD*polphi);
  zcospol = cos(DEG2RAD*polphi);

  zlampol = DEG2RAD*pollam;
  zphirot = DEG2RAD*phirot;
  if ( lamrot > 180.0 ) lamrot -= 360.0;
  zlamrot = DEG2RAD*lamrot;

  if ( polgam > 0 )
    {
      zgam  = DEG2RAD*polgam;
      zarg1 = sin(zlampol) *                                               
 	    (- zsinpol*cos(zphirot) * (cos(zlamrot)*cos(zgam) - sin(zlamrot)*sin(zgam)) 
 	     + zcospol*sin(zphirot))                                              
	- cos(zlampol)*cos(zphirot) * (sin(zlamrot)*cos(zgam) + cos(zlamrot)*sin(zgam));

      zarg2 = cos(zlampol) *                                               
 	    (- zsinpol*cos(zphirot) * (cos(zlamrot)*cos(zgam) - sin(zlamrot)*sin(zgam)) 
	     + zcospol*sin(zphirot))                                              
	+ sin(zlampol)*cos(zphirot) * (sin(zlamrot)*cos(zgam) + cos(zlamrot)*sin(zgam));
      }
  else
    {
      zarg1 = sin(zlampol)*(- zsinpol*cos(zlamrot)*cos(zphirot)  +
      		              zcospol*             sin(zphirot)) -
	      cos(zlampol)*           sin(zlamrot)*cos(zphirot);
      zarg2 = cos(zlampol)*(- zsinpol*cos(zlamrot)*cos(zphirot)  +
                              zcospol*             sin(zphirot)) +
              sin(zlampol)*           sin(zlamrot)*cos(zphirot);
    }

  if ( fabs(zarg2) > 0 ) result = RAD2DEG*atan2(zarg1, zarg2);
  if ( fabs(result) < 9.e-14 ) result = 0;

  return (result);
}


double phirot_to_phi(double phirot, double lamrot, double polphi, double polgam)
{
  /*
    This function converts phi from one rotated system to phi in another
    system. If the optional argument polgam is present, the other system
    can also be a rotated one, where polgam is the angle between the two
    north poles.
    If polgam is not present, the other system is the real geographical
    system.

    phirot : latitude in the rotated system
    lamrot : longitude in the rotated system (E>0)
    polphi : latitude of the rotated north pole
    polgam : angle between the north poles of the systems

    result : latitude in the geographical system
  */
  double zsinpol, zcospol;
  double zphirot, zlamrot, zarg;
  double zgam;

  zsinpol = sin(DEG2RAD*polphi);
  zcospol = cos(DEG2RAD*polphi);

  zphirot   = DEG2RAD*phirot;
  if ( lamrot > 180.0 ) lamrot -= 360.0;
  zlamrot   = DEG2RAD*lamrot;

  if ( polgam > 0 )
    {
      zgam = DEG2RAD*polgam;
      zarg = zsinpol*sin(zphirot) +
             zcospol*cos(zphirot)*(cos(zlamrot)*cos(zgam) - sin(zgam)*sin(zlamrot));
    }
  else
    zarg   = zcospol*cos(zphirot)*cos(zlamrot) + zsinpol*sin(zphirot);

  return (RAD2DEG*asin(zarg));
}

static
double rl_to_rls(double phi, double rla, double polphi, double pollam)
{
  /*
    Umrechnung von rla (geo. System) auf rlas (rot. System)

    phi    : Breite im geographischen System (N>0)
    rla    : Laenge im geographischen System (E>0)
    polphi : Geographische Breite des Nordpols des rot. Systems
    pollam : Geographische Laenge des Nordpols des rot. Systems

    result : Rotierte Laenge
  */
  double zsinpol, zcospol, zlampol;
  double zphi, zrla, zarg1, zarg2;

  zsinpol = sin(DEG2RAD*polphi);
  zcospol = cos(DEG2RAD*polphi);
  zlampol =     DEG2RAD*pollam;

  if ( rla > 180.0 ) rla -= 360.0;

  zrla = DEG2RAD*rla;
  zphi = DEG2RAD*phi;

  zarg1  = - sin(zrla-zlampol)*cos(zphi);
  zarg2  = - zsinpol*cos(zphi)*cos(zrla-zlampol)+zcospol*sin(zphi);

  if ( fabs(zarg2) < 1.0e-20 ) zarg2 = 1.0e-20;

  return (RAD2DEG*atan2(zarg1,zarg2));
}

static
double ph_to_phs(double phi, double rla, double polphi, double pollam)
{
  /*
    Umrechnung von phi (geo. System) auf phis (rot. System)

    phi    : Breite im geographischen System (N>0)
    rla    : Laenge im geographischen System (E>0)
    polphi : Geographische Breite des Nordpols des rot. Systems
    pollam : Geographische Laenge des Nordpols des rot. Systems

    result : Rotierte Breite
  */
  double zsinpol, zcospol, zlampol;
  double zphi, zrla, zarg;

  zsinpol = sin(DEG2RAD*polphi);
  zcospol = cos(DEG2RAD*polphi);
  zlampol =     DEG2RAD*pollam;

  zphi = DEG2RAD*phi;
  if ( rla > 180.0 ) rla -= 360.0;
  zrla = DEG2RAD*rla;

  zarg = zcospol*cos(zphi)*cos(zrla-zlampol) + zsinpol*sin(zphi);

  return (RAD2DEG*asin(zarg));
}


void usvs_to_uv(double us, double vs, double phi, double rla,
		double polphi, double pollam, double *u, double *v)
{
  /*
    Umrechnen der windkomponente us, vs im rotierten sphaerischen
    system in die windkomponenten u, v, im geographischen system

    us     : 'zonaler wind im rotierten system
    vs     : 'merid. wind im rotierten  system
    phi    : breite im geographischen system (n>0)
    rla    : laenge im geographischen system (e>0)
    polphi : geographische breite des n-pols des rot. sys.
    pollam : geographische laenge des n-pols des rot. sys.
 
    u      : zonaler wind im geographischen system
    v      : merid. wind im geographischen system
  */
  double zpolphi, zpollam, zrla, zphi, pollamd, zrlas, zarg, zbeta;

  /* umrechnung von grad in bogenmass */
  zpolphi = polphi*DEG2RAD;
  zpollam = pollam*DEG2RAD;
  zrla    = rla   *DEG2RAD;
  zphi    = phi   *DEG2RAD;
  pollamd = pollam;
  if ( pollamd < 0.0 ) pollamd += 360.0;

  /* laenge im rotierten system berechnen */
  zrlas = rl_to_rls(phi, rla, polphi, pollam)*DEG2RAD;

  /* winkel zbeta berechen (schnittwinkel der breitenkreise) */
  zarg = - sin(zpolphi)*sin(zrla-zpollam)*sin(zrlas) - cos(zrla-zpollam)*cos(zrlas);
  if ( zarg >  1.0 ) zarg =  1.0;
  if ( zarg < -1.0 ) zarg = -1.0;
  /*
  zbeta = acos(zarg);
  zbeta = sign(zbeta, -(rla - (pollamd-180.0)));
  */
  zbeta = fabs(acos(zarg));
  /*  if ( -(rla - (pollamd-180.0)) < 0 ) zbeta = -zbeta;*/
  if ( (-(rla - (pollamd-180.0)) < 0) && (-(rla - (pollamd-180.0)) >= -180) ) zbeta = -zbeta;

  /* us - wind transformieren */
  *u = us*cos(zbeta) - vs*sin(zbeta);
  
  /* vs - wind transformieren */
  *v = us*sin(zbeta) + vs*cos(zbeta);
}

/*
int main(void)
{
  double polphi, pollam;
  double x0, y0, x1, y1, x2, y2;
  int i;

  polphi = 90.0;
  pollam = 0.0;

  polphi = 32.5;
  pollam = -170.0;

  x0 = -20.0;
  y0 = 0.0;

  for ( i = 0; i < 10; i++ )
    {
      x0 = i *20.0;
      printf("%g %g\n", x0, y0);

      x1 = rls_to_rl(y0, x0, polphi, pollam);
      y1 = phs_to_ph(y0, x0, polphi);

      printf("%g %g\n", x1, y1);

      x2 = rl_to_rls(y1, x1, polphi, pollam);
      y2 = ph_to_phs(y1, x1, polphi, pollam);

      printf("%g %g\n", x2, y2);
    }

  usvs_to_uv(30.0, 20.0, 30.0, 0.0, polphi, pollam, &x1, &x2);
  printf("usvs_to_uv %g %g %g %g\n", polphi, pollam, x1, x2);

  return (0);
}
*/
