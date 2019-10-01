import numpy
from omuse.units import units
from amuse.units import trigo

sin=trigo.sin
cos=trigo.cos
tan=trigo.tan
atan2=trigo.arctan2
atan=trigo.arctan
sqrt=numpy.sqrt

def distance(lat1, lon1, lat2, lon2, R=1):
    """
    returns the spherical distance between two coordinates in given in lat,lon
    result should be multiplied by R for the actual distance in meters
    """

    dlat = lat2-lat1
    dlon = lon2-lon1

    #using haversine formula
    a = sin(dlat/2.0) * sin(dlat/2.0) + \
        cos(lat1) * cos(lat2) * \
        sin(dlon/2.0) * sin(dlon/2)

    d = 2 * atan2(sqrt(a), sqrt(1.0-a))

    return R*d

def triangle_area(a, b, c, R=1):
    """
    returns the area of the triangle enclosed by three lines of distances a, b, and c
    result should be multiplied with R**2 for the actual surface area in meters
    """

    #compute the semiperimeter
    s = (a + b + c)/2

    #compute spherical excess using L'Huilier's Theorem
    tans = tan(s/2)
    tana = tan((s-a)/2)
    tanb = tan((s-b)/2)
    tanc = tan((s-c)/2)

    tanE4 = sqrt(tans * tana * tanb * tanc)

    E = 4 * atan(tanE4)

    return R**2*E

if __name__=="__main__":
    R=1| units.Rearth
    lon1=0. | units.deg
    lat1=0. | units.deg
    lon2=1. | units.deg
    lat2=0. | units.deg
    print(distance(lat1,lon1,lat2,lon2, R=R).in_(units.km))

    lon1=0. | units.deg
    lat1=0. | units.deg
    lon2=90. | units.deg
    lat2=0. | units.deg
    lon3=0. | units.deg
    lat3=90. | units.deg    
    d1=distance(lat1,lon1,lat2,lon2)
    d2=distance(lat2,lon2,lat3,lon3)
    d3=distance(lat3,lon3,lat1,lon1)
    print(d1,d2,d3)
    print(triangle_area(d1,d2,d3,R=R)/4/trigo.pi/R**2)
