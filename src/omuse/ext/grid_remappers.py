from amuse.ext.grid_remappers import *

from amuse.datamodel.staggeredgrid import StaggeredGrid

class conservative_spherical_remapper(object):

    def __init__(self, source, target, axes_names=['lon', 'lat'], cdo_remapper=None):
        """ This class maps a source grid to a target grid using second-
            order conservative remapping by calling the re-implementation of
            SCRIP within CDO. The source grid should be a structured grid
            the target grid can be of any type. This class is able to deal
            with staggered grids for both source and target grid.
            Instantiating this class may take a while, as the remapping
            weights are being computed.
        """
        self.src_staggered = False
        self.source = source
        self.src_elements = source
        if type(source) is StaggeredGrid:
            self.src_staggered = True
            self.src_elements = source.elements
        if not type(self.src_elements) is StructuredGrid:
            raise Exception("Source grid should be of type StructuredGrid")

        self.tgt_staggered = False
        self.target = target
        self.tgt_elements = target
        if type(target) is StaggeredGrid:
            self.tgt_staggered = True
            self.tgt_elements = target.elements
    
        self._axes_names=list(axes_names)

        try:
            from omuse.community.cdo.interface import CDORemapper
        except:
            raise Exception("conservative spherical remapper requires omuse.community.cdo.interface")  

        if cdo_remapper == None:
            self.cdo_remapper = CDORemapper(channel="sockets", redirection="none")
            self.cdo_remapper.parameters.src_grid = self.src_elements
            self.cdo_remapper.parameters.dst_grid = self.tgt_elements

            #force start of the computation of remapping weights
            self.cdo_remapper.commit_parameters()
        else:
            self.cdo_remapper = cdo_remapper


    def _get_grid_copies_and_channel(self, source, target, attributes):
        source_copy=source.empty_copy()
        channel1=source.new_channel_to(source_copy)
        target_copy=target.empty_copy()
        channel2=target.new_channel_to(target_copy)
        channel3=target_copy.new_channel_to(target)

        channel1.copy_attributes(attributes)
        channel2.copy_attributes(self._axes_names)

        return source_copy, target_copy, channel3

    def forward_mapping(self, attributes):

        element_attributes = attributes
        node_attributes = []

        #if the grid is staggered split the list of attributes into node an element attributes
        if self.src_staggered:
            el_attr = self.source.elements.all_attributes()
            no_attr = self.source.nodes.all_attributes()
            element_attributes = set(el_attr).intersection(set(attributes))
            node_attributes = set(no_attr).intersection(set(attributes)).difference(element_attributes)

        self._forward_mapping_elements_to_elements(self.src_elements, self.tgt_elements, element_attributes)
        if len(node_attributes) > 0:
            self._forward_mapping_nodes_to_nodes(self.source, self.target, node_attributes)


    def _forward_mapping_elements_to_elements(self, source, target, attributes):

        #create in-memory copies of the grids and a channel to the target in-code grid
        source_copy, target_copy, channel3 = self._get_grid_copies_and_channel(source, target, attributes)

        #indices for interacting with CDORemapper
        index_i_src = list(range(source.size))         
        index_i_dst = list(range(target.size))         
       
        for attribute in attributes:
            #obtain source values and unit
            values=to_quantity( getattr(source_copy, attribute) )
            unit=values.unit
            values=numpy.array(values.number)

            #do the remapping
            self.cdo_remapper.set_src_grid_values(index_i_src, values.ravel('F'))
            self.cdo_remapper.perform_remap()
            result = self.cdo_remapper.get_dst_grid_values(index_i_dst).reshape(target.shape, order='F')

            #store result in copy target grid
            setattr(target_copy, attribute, (result if unit is units.none else (result | unit)))

        #push in-memory copy target grid to in-code storage grid
        channel3.copy_attributes(attributes)    


    def _forward_mapping_nodes_to_nodes(self, source, target, attributes):

        #create in-memory copies of the grids and a channel to the target in-code grid
        source_copy, target_copy, channel3 = self._get_grid_copies_and_channel(source.nodes, target.nodes, attributes)

        #indices for interacting with CDORemapper
        index_i_src = list(range(source.elements.size))
        index_i_dst = list(range(target.elements.size))
       
        for attribute in attributes:
            #obtain source values and unit
            values=to_quantity( getattr(source_copy, attribute) ) 
            unit=values.unit
            values=values.number
            if len(values.shape) > 1:
                values = numpy.swapaxes(values, 0, 1)

            #remap to elements within source grid
            values = source.map_nodes_to_elements(values)

            #do the remapping
            self.cdo_remapper.set_src_grid_values(index_i_src, values.flatten())
            self.cdo_remapper.perform_remap()
            result = self.cdo_remapper.get_dst_grid_values(index_i_dst)

            #remap to nodes within target grid
            result = target.map_elements_to_nodes(result)

            #store result in copy target grid
            setattr(target_copy, attribute, (result if unit is units.none else (result | unit)))

        #push in-memory copy target grid to in-code storage target grid
        channel3.copy_attributes(attributes)
