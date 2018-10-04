#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Wed May 23 11:48:50 2018

@author: Charlie
"""
import numpy as np
cimport numpy as np
import cython

DTYPE_FLOAT = np.double
ctypedef np.double_t DTYPE_FLOAT_t

DTYPE_INT = np.int
ctypedef np.int_t DTYPE_INT_t

@cython.boundscheck(False)
@cython.wraparound(False) 
def erode_bed(np.ndarray[DTYPE_INT_t, ndim=1] stack, 
              np.ndarray[DTYPE_FLOAT_t, ndim=1] z, 
              DTYPE_FLOAT_t bedrock__critical_stress, 
              np.ndarray[DTYPE_FLOAT_t, ndim=1] cover_fraction_array,
              np.ndarray[DTYPE_FLOAT_t, ndim=1] drag_stress_array,
              np.ndarray[DTYPE_FLOAT_t, ndim=1] flow_depth_array,
              DTYPE_FLOAT_t bedrock__erodibility,
              DTYPE_INT_t timestep,
              np.ndarray[DTYPE_FLOAT_t, ndim=1] downslope_link_length_array,
              np.ndarray[DTYPE_FLOAT_t, ndim=1] adjusted_shear_stress_array,
              DTYPE_INT_t baselevel_node):
    """
    erode_bed calculates vertical bedrock lowering given adjusted
    shear stress. It does so using the Fastscape algorithm of Braun and 
    Willett (2017), iterating up the "stack" from bottom to top.
    
    INPUTS
    ------
    stack: bottom-to-top list of nodes in drainage order
    br_elev: pointer to bedrock__elevation field
    bedrock__critical_stress: critical shear stress for bedrock erosion [Pa]
    cover_fraction_array: array of cover fractions
    drag_stress_array: array of dimensionless drag stresses [-]
    flow_depth_array: array of flow depths [m]
    bedrock__erodibility: bedrock "k" 
    timestep: model timestep [yr]
    downslope_link_length_array: array of lengths of downstream links
    baselevel_node: node number of baselevel node
    """
    #define internal variables
    cdef unsigned int node
    cdef float g
    cdef unsigned int dens_water
    cdef np.ndarray[DTYPE_FLOAT_t, ndim=1] excess_tau
    cdef np.ndarray[DTYPE_INT_t, ndim=1] stack_fsc
    cdef float f_open
    
    g = 9.81 #m/s2
    dens_water = 1000 #kg/m3
    excess_tau = adjusted_shear_stress_array - bedrock__critical_stress
    stack_fsc = np.insert(stack, 0, baselevel_node)        
    for node in range(1, len(stack_fsc)):
        if excess_tau[stack_fsc[node]] <= 0:
            pass #new elev is the same as old b/c no erosion
        else:
            f_open = 1 - cover_fraction_array[stack_fsc[node]]  
            #f_open = f_open.clip(min=0)
            z[stack_fsc[node]] = (z[stack_fsc[node]] + \
                (z[stack_fsc[node - 1]] * f_open * bedrock__erodibility * \
                dens_water * g * flow_depth_array[stack_fsc[node]] * timestep / \
                (downslope_link_length_array[stack_fsc[node]] * (1 + drag_stress_array[stack_fsc[node]]))) + (timestep * f_open * \
                 bedrock__erodibility * bedrock__critical_stress)) / (1 + (f_open * bedrock__erodibility * \
                 dens_water * g * flow_depth_array[stack_fsc[node]] * timestep / \
                 (downslope_link_length_array[stack_fsc[node]] * (1 + drag_stress_array[stack_fsc[node]]))))

@cython.boundscheck(False)
@cython.wraparound(False)            
def move_blocks(np.ndarray[np.int64_t, ndim=1] block_locations,
                np.ndarray[DTYPE_FLOAT_t, ndim=1] block_sizes,
                np.ndarray[DTYPE_FLOAT_t, ndim=1] slope,
                DTYPE_INT_t block_motion_threshold,
                DTYPE_FLOAT_t node_spacing,
                np.ndarray[DTYPE_FLOAT_t, ndim=1] soil_depth,
                np.ndarray[DTYPE_FLOAT_t, ndim=1] underlying_soil_depth,
                np.ndarray[DTYPE_FLOAT_t, ndim=1] bedrock_elevation,
                np.ndarray[DTYPE_FLOAT_t, ndim=1] topo_elevation,
                np.ndarray[np.int64_t, ndim=1] flow_receiver_node,
                np.ndarray[np.int64_t, ndim=1] core_nodes,
                np.ndarray[np.int64_t, ndim=1] channel_nodes,
                np.ndarray[np.int64_t, ndim=1] hillslope_nodes,
                np.ndarray[np.float64_t, ndim=1] blocks_to_channel_size,
                np.ndarray[np.int64_t, ndim=1] blocks_to_channel_location,
                DTYPE_INT_t cum_channel_blocks):
    """
    move blocks based on a user-defined relief threshold.        
    
    INPUTS
    ------
    flow_receiver_node: mg.at_node['flow__receiver_node']
    """
    #define internal variables
    cdef np.ndarray[DTYPE_FLOAT_t, ndim=1] steepest_slope_block_all
    #cdef np.ndarray[DTYPE_FLOAT_t, ndim=1] block_threshold_all
    cdef np.ndarray[DTYPE_FLOAT_t, ndim=1] block_threshold_slope_all
    cdef np.ndarray[DTYPE_INT_t, ndim=1] move_locations
    cdef DTYPE_FLOAT_t block_threshold
    cdef DTYPE_FLOAT_t block_threshold_slope
    cdef DTYPE_FLOAT_t steepest_slope_block
    cdef DTYPE_INT_t newlocation
    cdef np.ndarray[DTYPE_INT_t, ndim=1] move_locations_delete
    
    if block_locations.size:
        steepest_slope_block_all = slope[block_locations]
        #block_threshold_all = block_motion_threshold * block_sizes
        block_threshold_slope_all = block_motion_threshold * block_sizes / node_spacing
        
        move_locations = block_locations[steepest_slope_block_all >= block_threshold_slope_all]
           
        if move_locations.size:
            
            for l in (range(len(block_locations))):
                
                block_threshold = block_motion_threshold * block_sizes[l]
                block_threshold_slope = block_threshold / node_spacing
                steepest_slope_block = slope[block_locations[l]]
                
                if steepest_slope_block >= block_threshold_slope:
                    
                    soil_depth[block_locations[l]] = underlying_soil_depth[block_locations[l]]
                    bedrock_elevation[block_locations[l]] -= block_sizes[l] + underlying_soil_depth[block_locations[l]]
                    topo_elevation[block_locations[l]] = soil_depth[block_locations[l]] + bedrock_elevation[block_locations[l]]
                
                    newlocation = flow_receiver_node[block_locations[l]]
                    
                    if newlocation in hillslope_nodes:
                        block_locations[l] = newlocation
                    
                        underlying_soil_depth[newlocation] = soil_depth[newlocation]
                        bedrock_elevation[newlocation] += soil_depth[newlocation] + block_sizes[l]
                        soil_depth[newlocation] = 0
                
                
                        topo_elevation[newlocation] = soil_depth[newlocation] + bedrock_elevation[newlocation]
                    
                    else:    
                        if newlocation in channel_nodes:
                            new_blocks_size = block_sizes[l]
                            #new_blocks_location = self.block_locations[l]
                            
                            blocks_to_channel_size = np.append(blocks_to_channel_size, 
                                                               new_blocks_size)
                            blocks_to_channel_location = np.append(blocks_to_channel_location, 
                                                                   newlocation)
                            #if newlocation not in channel_nodes:
                            #    sys.exit('move_blocks put block not in channel')
                            #print 'MOVE ' 
                            #print self.blocks_to_channel_size
                            #print self.blocks_to_channel_location
                            #print '------------------------------'
                            #print new_blocks_size
                            
                            cum_channel_blocks += 1
                            
                            block_locations[l] = -9999
                            block_sizes[l] = -9999
                        else:
                            block_locations[l] = -9999
                            block_sizes[l] = -9999
                        
                        underlying_soil_depth[newlocation] = 0

            move_locations_delete = np.where(block_locations == -9999)[0]
            block_locations = np.delete(block_locations, move_locations_delete)
            block_sizes = np.delete(block_sizes, move_locations_delete)
    return block_locations, block_sizes, cum_channel_blocks, blocks_to_channel_size, blocks_to_channel_location