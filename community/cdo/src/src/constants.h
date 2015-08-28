#ifndef _CONSTANTS_H
#define _CONSTANTS_H

/* Thermodynamical constants adopted from ECMWF IFS-Code */

#define  C_RKBOL         (1.380658e-23)     /* Boltzmann constant in J/K   */
#define  C_RNAVO         (6.0221367e+23)    /* Avogadro constant in 1/mol  */
#define  C_RMD           (28.9644)          /* molecular weight of dry air */
#define  C_RMV           (18.0153)          /* molecular weight of water vapor */
#define  C_R             (C_RKBOL * C_RNAVO)
#define  C_RV            (1000. * C_R / C_RMV)

#define  C_EARTH_RD      (1000. * C_R / C_RMD)
#define  C_EARTH_RADIUS  (6371000.0)        /* radius of the Earth in m */
#define  C_EARTH_GRAV    (9.80665)

#define  C_MARS_RD       (189.0 )
#define  C_MARS_RADIUS   (3400000.0)        /* radius of the Mars in m */
#define  C_MARS_GRAV     (3.7)

#define  C_RG            (1.0 / PlanetGrav)

#define  C_RCPV          (4.0 * C_RV)
#define  C_RETV          (C_RV / PlanetRD - 1.)
#define  C_RCW           (4218.)            /* specific water heat capacity ?? */
#define  C_RCS           (2106.)            /* specific ice heat capacity ?? */
#define  C_RTT           (273.16)           /* melting temperature of ice/snow */
#define  C_RLVTT         (2.5008e+6)        /* latent heat for vaporisation in J/kg */
#define  C_RLSTT         (2.8345e+6)        /* latent heat for sublimation in J/kg */
#define  C_RESTT         (611.14)
#define  C_RCPD          (3.5 * PlanetRD)

#define  C_TIMES_RHOH2O  (-333700000.0)

extern double PlanetRD;
extern double PlanetRadius;
extern double PlanetGrav;

#endif  /* _CONSTANTS_H */
