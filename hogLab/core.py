# -*- coding: utf-8 -*-
"""
Created on Tue Jan  9 13:06:21 2018

@author: Charlie and Rachel

Class implementing the Landlab version of Rachel's blocky hillslope model.
"""
from __future__ import division
import numpy as np
import sys

class HogLab(object):
    
    def __init__(self):
        pass
    
    def instantiate_retreat_rate_tracking_array(self, nr, num_timesteps):
    
        """
        Creates a 2-D array for saving erosion rate values,
        -Each column is a row number; bottom to top is left to right
        -Each row is a timestep; early to late is top to bottom
        """
        self.retreat_rate_tracking_array = np.zeros((num_timesteps, nr-2), 
                                              dtype = np.float64)
                                              
    
    def setup_resistant_layer(self, 
                              thickness,
                              block_height,
                              nrocks,
                              dip,
                              release_relief_threshold,
                              starting_point,
                              initial_elevation,
                              top_of_hard_layer,
                              bottom_of_hard_layer,
                              x_coordinate_array,
                              bedrock_elevation, 
                              node_spacing, 
                              core_nodes,
                              channel_nodes,
                              dx):
        """
        Sets up resistant layer, calculates number of blocks in each release 
        event. Sets up the 3-D dip line of the layer.
        
        INPUTS:
        -------
        thickness: dip-perpendicular thickness of hard layer [m]
        block_height:
        block_width:
        dip: [degrees]
        release_relief_threshold:
        starting_point:
        initial_elevation:
        top_of_hard_layer: mg.at_node['zb_hard_top']
        bottom_of_hard_layer: mg.at_node['zb_hard_bottom']
        x_coordinate_array: mg.node_x
        bedrock_elevation: mg.at_node['bedrock__elevation']
        """
        #setup hillslope vs. channel nodes
        self.hillslope_nodes = np.copy(core_nodes)
        self.blocks_to_channel_size = np.array([], dtype = np.float64)
        self.blocks_to_channel_location = np.array([], dtype = np.int64)
        for fish in range(len(self.hillslope_nodes)):
            if self.hillslope_nodes[fish] in channel_nodes:
                self.hillslope_nodes[fish] = -9999
        hillslope_nodes_delete = np.where(self.hillslope_nodes == -9999)[0]
        self.hillslope_nodes = np.delete(self.hillslope_nodes, hillslope_nodes_delete)
            
        
        
        #nrocks = int((np.power(node_spacing, 2) * block_height) / (block_height * np.power(block_width, 2)))
        slope_threshold = release_relief_threshold / node_spacing
        
        #convert slope to radians
        slope_rad = np.radians(dip)
        
        block_width = dx
        
        #define resistant layer location
        thickness = thickness * block_height * np.cos(slope_rad) * block_width
        #starting_point = ncols*dx/2 #x values where layer begins
        intercept = initial_elevation - (np.tan(slope_rad) * starting_point)
        top_of_hard_layer[:] = slope_rad * x_coordinate_array + intercept
        bottom_of_hard_layer[:] = top_of_hard_layer[:] - \
            (thickness / np.cos(slope_rad))
        #hard_layer = np.where((bedrock_elevation <= top_of_hard_layer) & \
        #    (bedrock_elevation >= bottom_of_hard_layer))[0]

        return (slope_threshold, slope_rad) #hard_layer
    
        
    def setup_horizontal_layer(self, 
                              thickness,
                              block_height,
                              nrocks,
                              release_relief_threshold,
                              starting_point,
                              initial_elevation,
                              top_of_hard_layer,
                              bottom_of_hard_layer,
                              x_coordinate_array,
                              bedrock_elevation, 
                              node_spacing, 
                              core_nodes,
                              channel_nodes):
        """
        Sets up horizontal resistant layer, calculates number of blocks in each release 
        event. 
        
        INPUTS:
        -------
        thickness: dip-perpendicular thickness of hard layer [m]
        block_height:
        block_width:
        release_relief_threshold:
        starting_point:
        initial_elevation:
        top_of_hard_layer: mg.at_node['zb_hard_top']
        bottom_of_hard_layer: mg.at_node['zb_hard_bottom']
        x_coordinate_array: mg.node_x
        bedrock_elevation: mg.at_node['bedrock__elevation']
        """
        #setup hillslope vs. channel nodes
        self.hillslope_nodes = np.copy(core_nodes)
        self.blocks_to_channel_size = np.array([], dtype = np.float64)
        self.blocks_to_channel_location = np.array([], dtype = np.int64)
        for fish in range(len(self.hillslope_nodes)):
            if self.hillslope_nodes[fish] in channel_nodes:
                self.hillslope_nodes[fish] = -9999
        hillslope_nodes_delete = np.where(self.hillslope_nodes == -9999)[0]
        self.hillslope_nodes = np.delete(self.hillslope_nodes, hillslope_nodes_delete)
            
        
        slope_rad = 0
        #nrocks = int((np.power(node_spacing, 2) * block_height) / (block_height * np.power(block_width, 2)))
        slope_threshold = release_relief_threshold / node_spacing
        
        
        #define resistant layer location
        top_of_hard_layer[:] = initial_elevation
        bottom_of_hard_layer[:] = initial_elevation - thickness


        return (slope_threshold,slope_rad) #hard_layer
   
        
        
        
    def setup_block_tracking_arrays(self):
        """
        set up block location array and block size array.
        """
        self.block_locations = np.zeros(shape = 0, dtype = np.int64)
        self.block_sizes = np.zeros(shape = 0)
        self.hillslope_nodes = np.zeros(shape = 0, dtype = np.int64)
        self.hard_layer = np.zeros(shape = 0, dtype = np.int64)
        self.blocks_to_channel_location = np.zeros(shape = 0, dtype = np.int64)
        self.blocks_to_channel_size = np.zeros(shape = 0, dtype = np.float64)
        self.cum_channel_blocks = 0
        self.cum_release_blocks = 0
                
    def find_release_location(self,
                              bedrock_elevation,
                              top_of_hard_layer,
                              bottom_of_hard_layer,
                              slope,
                              slope_threshold,
                              channel_nodes):
        """
        re-evaluate location of hard layer and find locations of blocks
        that will fall.
        
        INPUTS
        ------
        top_of_hard_layer: mg.at_node['zb_hard_top'] [m]
        bottom_of_hard_layer: mg.at_node['zb_hard_bottom'] [m]
        bedrock_elevation: mg.at_node['bedrock__elevation'] [m]
        slope: mg.at_node['topographic__steepest_slope'] [-]
        """
        self.hard_layer = np.where((bedrock_elevation <= top_of_hard_layer) & \
            (bedrock_elevation >= bottom_of_hard_layer))[0]
        
        if self.hard_layer.size:
            for r in range(len(self.hard_layer)):
                if self.hard_layer[r] in channel_nodes:
                    self.hard_layer[r] = -9999
                    
            
            hard_layer_delete = np.where(self.hard_layer == -9999)[0]
            self.hard_layer = np.delete(self.hard_layer, hard_layer_delete)
            
        
        
        if self.hard_layer.size:
            steepest_slope_hard = slope[self.hard_layer]       
            drop_locations = self.hard_layer[steepest_slope_hard >= slope_threshold]
            
#        drop_receive_locations = flow_receiver_node[drop_locations]
#        sorted_drop_receive_locations = np.sort(drop_receive_locations)
#        drop_receive_locations_diff = np.diff(sorted_drop_receive_locations)
#        
#        duplicate_locations = np.where(drop_receive_locations_diff == 0)[0]
        
        
        
        
        return (drop_locations)
        
    def move_blocks(self,
                    slope,
                    block_motion_threshold,
                    node_spacing,
                    soil_depth,
                    underlying_soil_depth,
                    bedrock_elevation,
                    topo_elevation,
                    flow_receiver_node,
                    core_nodes,
                    channel_nodes):
        """
        move blocks based on a user-defined relief threshold.        
        
        INPUTS
        ------
        flow_receiver_node: mg.at_node['flow__receiver_node']
        """
        if self.block_locations.size:
            print 'blocks to move: ', self.block_locations.size
            steepest_slope_block_all = slope[self.block_locations]
            
            block_threshold_all = block_motion_threshold * self.block_sizes
            block_threshold_slope_all = block_threshold_all / node_spacing
            
            move_locations = self.block_locations[steepest_slope_block_all >= block_threshold_slope_all]
               
            if move_locations.size:
                for l in (range(len(self.block_locations))):
                    
                    block_threshold = block_motion_threshold * self.block_sizes[l]
                    block_threshold_slope = block_threshold / node_spacing
                    steepest_slope_block = slope[self.block_locations[l]]
                    
                    if steepest_slope_block >= block_threshold_slope:
                        
                        soil_depth[self.block_locations[l]] = underlying_soil_depth[self.block_locations[l]]
                        bedrock_elevation[self.block_locations[l]] -= self.block_sizes[l] + underlying_soil_depth[self.block_locations[l]]
                        topo_elevation[self.block_locations[l]] = soil_depth[self.block_locations[l]] + bedrock_elevation[self.block_locations[l]]
                    
                        newlocation = flow_receiver_node[self.block_locations[l]]
                        
                        if newlocation in self.hillslope_nodes:
                            self.block_locations[l] = newlocation
                        
                            underlying_soil_depth[newlocation] = soil_depth[newlocation]
                            bedrock_elevation[newlocation] += soil_depth[newlocation] + self.block_sizes[l]
                            soil_depth[newlocation] = 0
                    
                    
                            topo_elevation[newlocation] = soil_depth[newlocation] + bedrock_elevation[newlocation]
                        
                        else:    
                            if newlocation in channel_nodes:
                                new_blocks_size = self.block_sizes[l]
                                #new_blocks_location = self.block_locations[l]
                                
                                self.blocks_to_channel_size = np.append(self.blocks_to_channel_size, 
                                                                   new_blocks_size)
                                self.blocks_to_channel_location = np.append(self.blocks_to_channel_location, 
                                                                       newlocation)
                            if newlocation not in channel_nodes:
                                sys.exit('move_blocks put block not in channel')
                                #print 'MOVE ' 
                                #print self.blocks_to_channel_size
                                #print self.blocks_to_channel_location
                                #print '------------------------------'
                                #print new_blocks_size
                                
                                self.cum_channel_blocks += 1
                                
                                self.block_locations[l] = -9999
                                self.block_sizes[l] = -9999
                            else:
                                self.block_locations[l] = -9999
                                self.block_sizes[l] = -9999
                            
                            underlying_soil_depth[newlocation] = 0

                move_locations_delete = np.where(self.block_locations == -9999)[0]
                self.block_locations = np.delete(self.block_locations, move_locations_delete)

                self.block_sizes = np.delete(self.block_sizes, move_locations_delete)
                
    def release_blocks(self,
                       topo_elevation,
                       bedrock_elevation,
                       slope_threshold,
                       core_nodes,
                       underlying_soil_depth,
                       soil_depth,
                       slope_rad,
                       drop_locations,
                       block_height,
                       flow_receiver_node,
                       channel_nodes,
                       slope, release_relief_threshold):
        """
        release blocks from the hard layer onto the hillslope.        
        
        INPUTS
        ------
        """        
        if drop_locations.size:
            
            #print drop_locations
            
            for n in range(len(drop_locations)):
                
                temp_relief = topo_elevation[drop_locations[n]] - topo_elevation[flow_receiver_node[drop_locations[n]]]
                
                if temp_relief >= release_relief_threshold:

                    self.cum_release_blocks += 1

                    #remove blocks from drop locations
                    topo_elevation[drop_locations[n]] -= block_height / np.cos(slope_rad)
                    bedrock_elevation[drop_locations[n]] -= block_height / np.cos(slope_rad)

                    #retreat layer edge by 1
                    #hardlayer[steepest_slope_hard >= slope_threshold] = np.empty(len(hardlayer[steepest_slope_hard >= slope_threshold]))


                    self.hard_layer[np.where(self.hard_layer == drop_locations[n])[0]] = -9999
                    
                    block_locations_new = flow_receiver_node[drop_locations[n]] #sets block locations in next downhill steepest cell

                    if block_locations_new in channel_nodes:
                        new_blocks_size = block_height
                        #new_blocks_location = block_locations_new[n]
                        self.blocks_to_channel_size = np.append(self.blocks_to_channel_size, 
                                                       new_blocks_size)
                        self.blocks_to_channel_location = np.append(self.blocks_to_channel_location, 
                                                           block_locations_new)
                        if block_locations_new not in channel_nodes:
                            sys.exit('release_blocks put block not in channel')
                        #print 'release '
                        #print self.blocks_to_channel_size
                        #print self.blocks_to_channel_location
                        #print '------------------------------'
                        
                        self.cum_channel_blocks += 1
                        
                    if block_locations_new in self.hillslope_nodes:
                        
                        topo_elevation[block_locations_new] += block_height
                        bedrock_elevation[block_locations_new] = topo_elevation[block_locations_new]

                        #Save current soil depth
                        underlying_soil_depth[block_locations_new] += soil_depth[block_locations_new]

                        #reassign soil depth to be zero where blocks land
                        soil_depth[block_locations_new] = 0
                        
                        self.block_locations = np.append(self.block_locations, block_locations_new)
                        block_sizes_new = block_height
                        self.block_sizes = np.append(self.block_sizes,block_sizes_new)
                        
            

                        

#                    for turtle in range(len(self.hard_layer)):
#            
#                        if self.hard_layer[turtle] in drop_locations:
#                            self.hard_layer[turtle] = -9999
                
            hard_layer_delete = np.where(self.hard_layer == -9999)[0]
            self.hard_layer = np.delete(self.hard_layer, hard_layer_delete)
            
            
            
#
#            block_locations_new = flow_receiver_node[drop_locations] #sets block locations in next downhill steepest cell
#            
#            for n in range(len(block_locations_new)):
#                
#
#              
#                
#                if block_locations_new[n] in channel_nodes:
#                    new_blocks_size = block_height
#                    new_blocks_location = block_locations_new[n]
#                    self.blocks_to_channel_size = np.append(self.blocks_to_channel_size, 
#                                                       new_blocks_size)
#                    self.blocks_to_channel_location = np.append(self.blocks_to_channel_location, 
#                                                           new_blocks_location)
#                    
##                    print self.blocks_to_channel_size
##                    print self.blocks_to_channel_location
##                    print new_blocks_size
#                    
#                    self.cum_channel_blocks += 1
#                    
#                                                           
#                if block_locations_new[n] in self.hillslope_nodes:
#                    pass
#                else:
#                    block_locations_new[n] = -9999
#            
#            block_locations_new_delete = np.where(block_locations_new == -9999)
#            block_locations_new = np.delete(block_locations_new, block_locations_new_delete)
#            
#                    
#            topo_elevation[block_locations_new] += block_height
#            bedrock_elevation[block_locations_new] = topo_elevation[block_locations_new]
#
#            #Save current soil depth
#            underlying_soil_depth[block_locations_new] += soil_depth[block_locations_new]
#
#            #reassign soil depth to be zero where blocks land
#            soil_depth[block_locations_new] = 0
#
#            self.block_locations = np.concatenate((self.block_locations,block_locations_new),axis = 0)
#            block_sizes_new = np.zeros(len(block_locations_new))
#            block_sizes_new += block_height
#            self.block_sizes = np.concatenate((self.block_sizes,block_sizes_new),axis = 0)



    def set_soil_production_rate_on_hard_layer(self,
                                           soil_production_rate):
            """
            INPUTS
            ------
            soil_production_rate: mg.at_node['soil_production__rate'] [m/yr]
            """
            soil_production_rate[self.hard_layer] = 0
    
    def weather_blocks(self,
                       soil_production_rate,
                       max_block_weathering_rate,
                       underlying_soil_depth,
                       soil_depth,
                       timestep):
        """
        weather blocks and delete them when they're gone.
        
        INPUTS
        ------
        
        """
        #Take control of block weathering rate
        if self.block_locations.size:
            soil_production_rate[self.block_locations.astype(int)] = max_block_weathering_rate ##%% change to include soil depth on top of blocks?
            self.block_sizes -= max_block_weathering_rate * timestep
            blocks_delete = np.where(self.block_sizes <= 0)[0]
            self.block_sizes = np.delete(self.block_sizes, blocks_delete)
            soil_depth[self.block_locations[blocks_delete]] = \
                underlying_soil_depth[self.block_locations[blocks_delete]] + \
                soil_depth[self.block_locations[blocks_delete]]
            underlying_soil_depth[self.block_locations[blocks_delete]] = 0
            self.block_locations = np.delete(self.block_locations, blocks_delete)
            
    def calculate_layer_retreat_rate(self, old, new, dt, dx, row):
        """
        Calculates channel incision rate by differencing channel elevations
        between timesteps
        """
        
        retreat_rate = (old - new)*dx / dt
        retreat_rate[retreat_rate<0] = 0
        
        self.retreat_rate_tracking_array[row, :] = retreat_rate
        
    def save_out_retreat_rate_array(self, directory, name_prefix):
        """
        At the end of the model run, save the matrix containing layer retreat rates from every node and timestep to a .npy binary.
        """
        np.save(directory + '/' + name_prefix + '_retreat_rate_record.npy', self.retreat_rate_tracking_array)