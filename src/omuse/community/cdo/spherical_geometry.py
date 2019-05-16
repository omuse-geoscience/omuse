import numpy
from math import *

#radius of Earth in meters
R = 6371000

#returns the spherical distance between two coordinates in given in lat,lon
#result should be multiplied by R for the actual distance in meters
def distance(lat1, lon1, lat2, lon2):
    dlat = lat2-lat1
    dlon = lon2-lon1

    #using haversine formula
    a = sin(dlat/2.0) * sin(dlat/2.0) + \
        cos(lat1) * cos(lat2) * \
        sin(dlon/2.0) * sin(dlon/2)

    d = 2 * atan2(sqrt(a), sqrt(1.0-a))

    return d

#returns the area of the triangle enclosed by three lines of distances a, b, and c
#result should be multiplied with R**2 for the actual surface area in meters
def triangle_area(a, b, c):
    #compute the semiperimeter
    s = (a + b + c)/2

    #compute spherical excess using L'Huilier's Theorem
    tans = tan(s/2)
    tana = tan((s-a)/2)
    tanb = tan((s-b)/2)
    tanc = tan((s-c)/2)

    tanE4 = sqrt(tans * tana * tanb * tanc)

    E = 4 * atan(tanE4)

    return E


def triangle_area_points(lon1, lat1, lon2, lat2, lon3, lat3):
    a = distance(lat1, lon1, lat2, lon2)
    b = distance(lat2, lon2, lat3, lon3)
    c = distance(lat3, lon3, lat1, lon1)
    return triangle_area(a, b, c)
