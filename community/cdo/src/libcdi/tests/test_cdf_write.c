#ifdef HAVE_CONFIG_H
#  include "config.h"
#endif

#include <math.h>
#include <stdlib.h>
#include <stdio.h>

#include "cdi.h"

struct cart_coord {
  double lat, lon;
};

static double
my_gamma_dist(double x);

static void
compute_curvilinear(double *coords_,
                    const struct cart_coord a[4], size_t sizex, size_t sizey);

#define DEG2RAD(phi) ((M_PI / 180.0) * (phi))
#define RAD2DEG(phi) ((180./M_PI) * (phi))
#define missValue (-50.0)

int main(int argc, const char **argv)
{
  /* todo: handle optional arguments here to increase test coverage */
  const char *fname = "test.nc";
  if (argc > 1)
    fname = argv[1];
  int streamID = streamOpenWrite(fname, FILETYPE_NC);

  if ( streamID < 0 )
    {
      fprintf(stderr, "Open failed on %s: %s\n", fname,
              cdiStringError(streamID));
      return EXIT_FAILURE;
    }

  enum {
    sizey = 40,
    sizex = 2 * sizey,
  };

  size_t datasize = (size_t)sizex * (size_t)sizey;
  int gridID = gridCreate(GRID_CURVILINEAR, (int)datasize);
  gridDefXsize(gridID, sizex);
  gridDefYsize(gridID, sizey);
  {
    /* anti-clockwise coordinates around Amazonia */
    static struct cart_coord region[4]
      = { { .lon = DEG2RAD(-85.0), .lat = DEG2RAD(-25.0) },
          { .lon = DEG2RAD(-44.0), .lat = DEG2RAD(-18.0) },
          { .lon = DEG2RAD(-50.0), .lat = DEG2RAD(7.0) },
          { .lon = DEG2RAD(-80.0), .lat = DEG2RAD(10.0) } };
    double (*gridCoords)[sizey][sizex]
      = (double (*)[sizey][sizex])
      malloc(sizeof (*gridCoords) * sizey * sizex * 2);
    if (gridCoords == NULL)
      {
        perror("grid coordinate memory allocation failed");
        return EXIT_FAILURE;
      }
    compute_curvilinear((double *)gridCoords, region, sizex, sizey);
    gridDefXvals(gridID, (double *)(gridCoords[1]));
    gridDefYvals(gridID, (double *)(gridCoords[0]));
    free(gridCoords);
  }

  int zaxisID = zaxisCreate(ZAXIS_SURFACE, 1);

  int vlistID = vlistCreate();
  int varID = vlistDefVar(vlistID, gridID, zaxisID, TSTEP_INSTANT);
  vlistDefVarMissval(vlistID, varID, missValue);
  {
    static const char creatorText[] = "CDI test_cdf_write";
    vlistDefAttTxt(vlistID, varID, "CDI Text Attribute test, created by",
                   sizeof (creatorText) - 1, creatorText);
  }

  int taxisID = taxisCreate(TAXIS_ABSOLUTE);
  vlistDefTaxis(vlistID, taxisID);

  streamDefVlist(streamID, vlistID);

  (void)streamDefTimestep(streamID, 0);

  {
    double (*data)[sizex] = malloc(sizeof (**data) * sizex * sizey);
    if (!data)
      {
        perror("data values memory allocation failed");
        return EXIT_FAILURE;
      }
    for (size_t j = 0; j < sizey; ++j)
      for (size_t i = 0; i < sizex; ++i)
        {
          data[j][i] = my_gamma_dist((double)i/(double)(sizex - 1));
        }
    data[sizey/3][sizex/2] = missValue;
    streamWriteVar(streamID, 0, (const double *)data, 1);
    free(data);
  }

  streamClose(streamID);

  return EXIT_SUCCESS;
}

static inline double
cart_distance(struct cart_coord p1, struct cart_coord p2)
{
  double d_lat = sin((p1.lat-p2.lat)/2.0),
    d_lon = sin((p1.lon - p2.lon)/2.0),
    d = 2.0 * asin(sqrt(d_lat * d_lat +
                        cos(p1.lat)*cos(p2.lat) * (d_lon * d_lon)));
  return d;
}

static inline struct cart_coord
intermediate_coord(struct cart_coord p1, struct cart_coord p2, double f)
{
  double d = cart_distance(p1, p2),
    sine_of_d = sin(d),
    A = sin((1 - f) * d) / sine_of_d,
    B = sin(f * d) / sine_of_d,
    x = A * cos(p1.lat) * cos(p1.lon) + B * cos(p2.lat) * cos(p2.lon),
    y = A * cos(p1.lat) * sin(p1.lon) + B * cos(p2.lat) * sin(p2.lon),
    z = A * sin(p1.lat) + B * sin(p2.lat);
  struct cart_coord ic = { .lat = atan2(z, sqrt(x * x + y * y)),
                           .lon = atan2(y, x) };
  return ic;
}

static void
compute_curvilinear(double *coords_,
                    const struct cart_coord a[4], size_t sizex, size_t sizey)
{
  double (*coords)[sizey][sizex] = (double (*)[sizey][sizex])coords_;
  for (size_t j = 0; j < sizey; ++j)
    {
      double g = (double)j / (double)(sizey - 1);
      /* compute start/end coordinates of great circle in x direction */
      struct cart_coord gc_left = intermediate_coord(a[0], a[3], g),
        gc_right = intermediate_coord(a[1], a[2], g);
      for (size_t i = 0; i < sizex; ++i)
        {
          double f = (double)i / (double)(sizex - 1);
          struct cart_coord pij = intermediate_coord(gc_left, gc_right, f);
          coords[0][j][i] = RAD2DEG(pij.lat);
          coords[1][j][i] = RAD2DEG(pij.lon);
        }
    }
}

static double
my_gamma_dist(double x)
{
  enum {
    k = 9,
  };
  const double theta = 0.5;
  x *= 20.0;
  double pdf_x = 1.0 / ( tgamma(k) * pow(theta, k) ) * pow(x, k - 1)
    * exp(-x/theta);
  return pdf_x;
}

/*
 * Local Variables:
 * c-file-style: "Java"
 * c-basic-offset: 2
 * indent-tabs-mode: nil
 * show-trailing-whitespace: t
 * require-trailing-newline: t
 * End:
 */
