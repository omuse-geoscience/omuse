import numpy as np


def get_haversine(lon1, lat1, lon2, lat2, d2r):

    dlat = d2r * (lat2 - lat1)
    dlon = d2r * (lon2 - lon1)
    lt1 = d2r * lat1
    lt2 = d2r * lat2
    
#    a = np.sin(0.5 * dlon) * np.sin(0.5 * dlon)
#    a = a * np.cos(lt1) * np.cos(lt2)
#    a = a + (np.sin(0.5 * dlat) * np.sin(0.5 * dlat))
#    thedist = 2.0 * np.arctan2(np.sqrt(a), np.sqrt(1.0 - a))

    a = np.sin(0.5 * dlon) * np.sin(0.5 * dlon) * \
            np.cos(lt1) * np.cos(lt2) + \
            (np.sin(0.5 * dlat) * np.sin(0.5 * dlat))
    thedist = 2.0 * np.arctan2(np.sqrt(a), np.sqrt(1.0 - a))
    
    return thedist
    


def distance_vector(lon1, lat1, lon2, lat2):

    erad = 6371315.0
    d2r = np.arctan2(0.,-1.)/ 180. # atan2(0.,-1.) == pi
    
    m = np.shape(lon1)[0]
    
    dist = np.zeros((m))
    for i in range(m):
        dist[i] = get_haversine(lon1[i], lat1[i], lon2[i], lat2[i], d2r)
        
    dist = dist * erad
    
    return dist


    
def distance_matrix(xa, xb):

    erad = 6371315.0
    d2r = np.arctan2(0.,-1.)/ 180. # atan2(0.,-1.) == pi
    
    m = np.shape(xa)[1] #changed by Ben, used to be 0
    n = np.shape(xb)[1] #changed by Ben, used to be 0

    dist = np.zeros((m,n))

    for j in range(m):
        #for i in range(n):
            #dist[j,i] = get_haversine(xa[j,0], xa[j,1], xb[i,0], xb[i,1],d2r)  #fiona probably had xa and xb zipped
            #dist[j,i] = get_haversine(xa[0,j], xa[1,j], xb[0,i], xb[1,i], d2r)  #first edit by Ben
        dist[j,:] = get_haversine(xa[0,j], xa[1,j], xb[0,:], xb[1,:], d2r)  #Ben's attempt to do entire row at once

        

            
    dist = dist * erad
    
    return dist


def distance(lon1, lat1, lon2, lat2):

    erad = 6371315.0
    d2r = np.arctan2(0.,-1.)/ 180. # atan2(0.,-1.) == pi
    
    thedist = get_haversine(lon1, lat1, lon2, lat2, d2r)
    thedist = thedist * erad
    
    return thedist
