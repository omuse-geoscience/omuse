#!/usr/bin/env python

import numpy

from spherical_geometry import *

class adcirc_grid_converter(object):

    def __init__(self, filename, coordinates):

        from omuse.community.adcirc.read_grid import adcirc_grid_reader

        r = adcirc_grid_reader(filename, coordinates)
        r.read_grid()
        self.num_elems = num_elems = r.parameters['NE']
        self.num_nodes = num_nodes = r.parameters['NP']

        self.size = num_elems
        self.corners = num_corners = 3  #assuming only triangles are used for now

        # r.p contains the lines of node definitions, strip last column
        self.lon = lon = r.p[:,0]
        self.lat = lat = r.p[:,1]

        # r.t contains the connectivity of nodes to form triangles, strip fist column 
        self.triangles = triangles = r.t[:,1:] -1  #minus one to offset Fortran indexes starting at 1

        self.grid_center_lat = grid_center_lat = numpy.zeros(num_elems, dtype=numpy.double)
        self.grid_center_lon = grid_center_lon = numpy.zeros(num_elems, dtype=numpy.double)

        self.grid_corner_lat = grid_corner_lat = numpy.zeros(num_elems*num_corners, dtype=numpy.double)
        self.grid_corner_lon = grid_corner_lon = numpy.zeros(num_elems*num_corners, dtype=numpy.double)

        self.inverse_mapping = inverse_mapping = [[] for i in range(num_nodes)]

        self.elem_area = elem_area = numpy.zeros(num_elems, dtype=numpy.double)
        self.elem_area_spherical = elem_area_spherical = numpy.zeros(num_elems, dtype=numpy.double)

        i=0
        for triangle in triangles:
            n1,n2,n3 = triangle
            lat1 = lat[n1]
            lat2 = lat[n2]
            lat3 = lat[n3]
            lon1 = lon[n1]
            lon2 = lon[n2]
            lon3 = lon[n3]

            #add element i as neighbors of nodes n1, n2, and n3
            inverse_mapping[n1].append(i)
            inverse_mapping[n2].append(i)
            inverse_mapping[n3].append(i)

            grid_center_lat[i] = (lat1 + lat2 + lat3)/3.0
            grid_center_lon[i] = (lon1 + lon2 + lon3)/3.0

            grid_corner_lat[i*3+0] = lat1
            grid_corner_lat[i*3+1] = lat2
            grid_corner_lat[i*3+2] = lat3

            grid_corner_lon[i*3+0] = lon1
            grid_corner_lon[i*3+1] = lon2
            grid_corner_lon[i*3+2] = lon3

            #triangle area, uses counter-clockwise sorting of nodes in trangle
            #assumes 2D euclidian geometry, not entirely accurate for lat-lon
            elem_area[i] = (((lon2-lon1)*(lat3-lat1))-((lon3-lon1)*(lat2-lat1)))/2.0

            #compute triangle area using spherical geometry
            a = distance(lat1, lon1, lat2, lon2)
            b = distance(lat2, lon2, lat3, lon3)
            c = distance(lat3, lon3, lat1, lon1)
            elem_area_spherical[i] = triangle_area(a, b, c)

            i+=1

        self.adjacency_list = adjacency_list = []
        self.boundary_node = boundary_node = [False for i in range(num_nodes)]

        for node in range(num_nodes):

            #create a list for each cell center connected to this corner
            cell_neighbors = [[] for i in inverse_mapping[node]]
            #use celli as index into cell_neighbors array
            celli = 0

            #for each neighboring cell/triangle
            for cell in inverse_mapping[node]:
                
                corners = triangles[cell][:] #copy corner list
                corner_list = list(corners)
                corner_list.remove(node) #remove the shared corner, look for another one
                for i in inverse_mapping[node]:
                    matches = 0
                    for c in triangles[i]:
                        if c in corner_list:
                            matches += 1
                    if matches == 1: #matches should be one, if it's 2 then it is the cell itself
                        cell_neighbors[celli].append(i)

                #if there is a cell with less than two neighboring cells for this corner, then
                #this corner must be on a boundary
                if len(cell_neighbors[celli]) < 2:
                    boundary_node[node] = True

                #print "node=", node, "cell=", cell, "neighbors=", cell_neighbors[celli] 
                celli += 1

            adjacency_list.append(cell_neighbors)



    def get_nodes_from_elements(self, elem_values):
        if (len(elem_values) != self.num_elems):
            print "Error: number of elements in grid does not match number of elements passed"

        node_values = numpy.zeros(self.num_nodes)

        for i in range(self.num_nodes):
            num_neighbors = len(self.inverse_mapping[i])
            value = 0.0

            #add value of neighboring element
            for neighbor in self.inverse_mapping[i]:
                value += elem_values[neighbor]

            #compute the averaged value
            value /= num_neighbors

            #store result
            node_values[i] = value

        return node_values


    def compute_node_fracarea(self, i, elem_values, elem_area):
            num_neighbors = len(self.inverse_mapping[i])
            value = 0.0
            area = 0.0

            #add value of neighboring element
            for neighbor in self.inverse_mapping[i]:
                #value += elem_values[neighbor] 
                fracarea = elem_area[neighbor]/3.0
                area += fracarea
                value += elem_values[neighbor] * fracarea
            #compute the (area) averaged value
            #value /= num_neighbors
            value /= area

            return value

    def get_nodes_from_elements_fracarea(self, elem_values, elem_area):
        if (len(elem_values) != self.num_elems):
            print "Error: number of elements in grid does not match number of elements passed"

        node_values = numpy.zeros(self.num_nodes)

        for i in range(self.num_nodes):
            node_values[i] = self.compute_node_fracarea(i, elem_values, elem_area)

        return node_values




    """
      This function implements an algorithm for redistributing values from elements to nodes. 
      It is adapted from the gathering stage for computing subcell momenta
      the original algorithm is described in "A subcell remapping method on staggered polygonal
      grids for arbitrary-Lagrangian-Eulerian methods" by R. Loubere and M. Shashkov
      Journal of Computational Physics, 2005, Volume 209, Pages 105--138
    """
    def get_nodes_from_elements_gather(self, elem_values, elem_area):
        if (len(elem_values) != self.num_elems):
            print "Error: number of elements in grid does not match number of elements passed"

        node_values = numpy.zeros(self.num_nodes)

        for i in range(self.num_nodes):

            #if node is on boundary we for the moment defer to fracarea method
            num_neighbors = len(self.inverse_mapping[i])
            value = 0.0
            area = 0.0

            celli = 0 #index within adjacency_list
            #add value of neighboring element
            for cell in self.inverse_mapping[i]:
                #value += elem_values[neighbor] 
                subcellarea = elem_area[cell]/3.0
                area += subcellarea
                cell_neighbors = self.adjacency_list[i][celli]

                if len(cell_neighbors) == 2:
                    subcellval = 2.0 * elem_values[cell] 
                    subcellval -= 0.5 * elem_values[cell_neighbors[0]]
                    subcellval -= 0.5 * elem_values[cell_neighbors[1]]
                if len(cell_neighbors) == 1:
                    subcellval = 1.5 * elem_values[cell] 
                    subcellval -= 0.5 * elem_values[cell_neighbors[0]]
                if len(cell_neighbors) == 0:
                    subcellval = 1.0 * elem_values[cell] 

                value += subcellval * subcellarea
                celli += 1

            #compute the (area) averaged value
            #value /= num_neighbors
            value /= area

            #store result
            node_values[i] = value

        return node_values









    def get_elements_from_nodes(self, node_values):
        if (len(node_values) != self.num_nodes):
            print "Error: number of nodes in grid does not match number of nodes passed"
        
        elem_values = numpy.zeros(self.num_elems)

        for i in range(self.num_elems):
            neighbors = self.triangles[i]
            value = ( node_values[neighbors[0]] + node_values[neighbors[1]] + node_values[neighbors[2]] ) / 3.0
            elem_values[i] = value

        return elem_values










    def area_sum(self, elem_values, elem_area):
        sum = numpy.zeros(1, dtype=numpy.double)
        for i in range(self.num_elems):
            sum += elem_values[i] * elem_area[i]
        return sum


    def print_area_sums(self, original_values, elem_values):
        sum0 = c.area_sum(original_values, self.elem_area)
        ssum0 = c.area_sum(original_values, self.elem_area_spherical)
        sum1 = c.area_sum(elem_values, self.elem_area)
        ssum1 = c.area_sum(elem_values, self.elem_area_spherical)
        print "sum = ", sum1, "spherical area sum= ", ssum1
        terror = (original_values.sum() - elem_values.sum())/original_values.sum()
        error = (sum0-sum1)/sum0
        serror = (ssum0-ssum1)/ssum0
        print "total value error=%e" % terror , "error based on euclidian area=%e" % error, ", error based on spherical area=%e" % serror
        


if __name__ == "__main__":
    filename = "grids/nc_inundation_v6c.grd"
    c = adcirc_grid_converter(filename, "spherical")


    print 'size=', c.size




    from matplotlib import pyplot

    ntest = 1

#    elem_values = numpy.zeros(c.num_elems) + 1
    elem_values = numpy.array(numpy.random.random(c.num_elems), dtype=numpy.double)
    original_values = elem_values[:] #hoping this creates a copy

    numpy.set_printoptions(precision=35)
    
    c.print_area_sums(original_values, elem_values)

    print "Computing ",ntest," back and forth transformations using number of neighbors weights"

    for i in range(ntest):
        node_values = c.get_nodes_from_elements(elem_values)
        elem_values = c.get_elements_from_nodes(node_values)

    c.print_area_sums(original_values, elem_values)

    #restore original values
    elem_values = original_values[:] #assuming this creates a copy

    print "Computing ",ntest," back and forth transformations using area weights (Euclidian)"

    for i in range(ntest):
        node_values = c.get_nodes_from_elements_fracarea(elem_values, c.elem_area)
        elem_values = c.get_elements_from_nodes(node_values)

    c.print_area_sums(original_values, elem_values)

    #restore original values
    elem_values = original_values[:] #assuming this creates a copy

    print "Computing ",ntest," back and forth transformations using area weights (Spherical)"

    for i in range(ntest):
        node_values = c.get_nodes_from_elements_fracarea(elem_values, c.elem_area_spherical)
        elem_values = c.get_elements_from_nodes(node_values)

    c.print_area_sums(original_values, elem_values)

    #restore original values
    spheric_elem_values = elem_values[:]
    spheric_node_values = node_values[:]
    elem_values = original_values[:] #assuming this creates a copy

    print "Computing ",ntest," back and forth transformations using area weights (Gather)"

    for i in range(ntest):
        node_values = c.get_nodes_from_elements_gather(elem_values, c.elem_area_spherical)
        elem_values = c.get_elements_from_nodes(node_values)

    c.print_area_sums(original_values, elem_values)

    #now plot both to see if theres any visual differences
    cmin = original_values.min()
    cmax = original_values.max()

    f, sp = pyplot.subplots(nrows=2, ncols=3, sharex=True, sharey=True)
    sp = sp.flatten()
    for im in sp:
        im.set_adjustable('box-forced')

    kwargs = dict(shading='flat', edgecolor='k', cmap=pyplot.cm.jet, vmin=cmin, vmax=cmax)

    sp[0].tripcolor(c.lon, c.lat, c.triangles, original_values, **kwargs)
    sp[1].tripcolor(c.lon, c.lat, c.triangles, spheric_elem_values, **kwargs)
    sp[2].tripcolor(c.lon, c.lat, c.triangles, elem_values, **kwargs)


    sp[3].tripcolor(c.lon, c.lat, c.triangles, (elem_values - spheric_elem_values)/original_values, **kwargs)

    kwargs['shading'] = 'gouraud'

    sp[4].tripcolor(c.lon, c.lat, c.triangles, spheric_node_values, **kwargs)
    sp[5].tripcolor(c.lon, c.lat, c.triangles, node_values, **kwargs)




    pyplot.tight_layout()
    pyplot.show()


 
    raw_input()
