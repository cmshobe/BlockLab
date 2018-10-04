# -*- coding: utf-8 -*-
"""
Created on Thu Oct 26 13:59:12 2017

@author: Charlie

Class implementing the Landlab version of BRaKE
"""

#Packages:
from __future__ import division
import numpy as np
import sys
import os

class BrakeLandlab(object):
    
    def __init__(self):
        pass
    
    def instantiate_e_rate_tracking_array(self, num_chan_nodes, num_timesteps):
        """
        Creates a 2-D array for saving erosion rate values,
        -Each column is a channel node; bottom to top is left to right
        -Each row is a timestep; early to late is top to bottom
        """
        self.e_rate_tracking_array = np.zeros((num_timesteps, num_chan_nodes), 
                                              dtype = np.float64)
    
    def roughness_bisection(self, q, s, d, g):
        """
        roughness_bisection calculates flow depth given discharge, slope, etc
        
        INPUTS
        ------
        q: water discharge per unit width [m2/s]
        s: topographic steepest slope at node
        d: roughness length scale [m]
        g: acceleration due to gravity [m/s2]
        """
        if s <= 0:
            print '0 or - slope, asserting no erosion'
            h = 0
        else:
            a1 = 6.5
            a2 = 2.5
            coef_1 = (g * s * np.power(a1, 2)) / (np.power(q, 2) * np.power(d, 2))
            coef_2 = 0.0
            coef_3 = -1 / np.power(d, 5 / 3)
            coef_4 = -np.power(a1 / a2, 2)
            coef_array = np.array([coef_1, coef_2, coef_3, coef_4])
            roots = np.roots(coef_array)
            is_real = np.isreal(roots)
            root = np.real(roots[is_real])[0]
            h = np.power(root, 3 / 5)
            if np.isnan(np.sum(h)):
                print h
                sys.exit("NAN FOUND IN ROOT'S H-- KILLING MODEL")
        return h
        
    def calc_shear_stress_with_roughness(self, pre_adjusted_tau,  
                                         is_block_in_cell, drag_cube, 
                                         flow_depth, dx, roughness_height, 
                                         cell,
                                         drag_stress_array,
                                         chan_width):
        """
        calc_shear_stress_with_roughness calculates adjusted shear stress
        given the size and number of blocks in a cell.
        
        INPUTS
        ------
        pre_adjusted_tau: the depth-slope product [Pa]
        tracking_mat: the block tracking matrix
        is_block_in_cell: boolean array telling whether blocks are in this node
        drag_cube: drag coefficient on a cube [-]
        flow_depth: water depth [m]
        dx: node spacing [m]
        roughness_height: bed roughness length scale [m]
        slicing index: index for selecting extant blocks in tracking matrix
        cell: number of node being worked on
        """
        #blocks_above_flow = (is_block_in_cell) #& (tracking_mat[0:slicing_index, 2] >= flow_depth)        
        a1 = 6.5 #ferguson 2007
        a2 = 2.5 #ferguson 2007        
        if flow_depth == 0:
            sigma_d_blocks = 0
            adjusted_shear_stress = 0
        else:
            if np.count_nonzero(is_block_in_cell) == 0:
                sigma_d_blocks = 0
            else:
                #blocks_above_flow = (is_block_in_cell) #& (tracking_mat[0:slicing_index, 2] >= flow_depth)
                beta = (a1 * (flow_depth / roughness_height)) / np.power(np.power(flow_depth / roughness_height, 5 / 3) + np.power(a1 / a2, 2), 1/2)
                avg_diam_blocks = np.average(self.tracking_mat[0:self.for_slicing, 1][is_block_in_cell]) #!
                submerged_block = (is_block_in_cell) & (self.tracking_mat[0:self.for_slicing, 1] < flow_depth)#!
                emergent_block = (is_block_in_cell) & (self.tracking_mat[0:self.for_slicing, 1] >= flow_depth)#!
                self.tracking_mat[0:self.for_slicing, 3][submerged_block] = self.tracking_mat[0:self.for_slicing, 1][submerged_block]#!
                self.tracking_mat[0:self.for_slicing, 3][emergent_block] = flow_depth           #! 
                avg_submerged_height_blocks = np.average(self.tracking_mat[0:self.for_slicing, 3][is_block_in_cell]) #!
                #avg_spacing_blocks =  dx / np.count_nonzero(is_block_in_cell)
                avg_spacing_blocks_squared = (dx * chan_width) / np.count_nonzero(is_block_in_cell)
                #sigma_d_blocks = (1 / 2) * drag_cube * np.power(beta, 2) *(avg_submerged_height_blocks * avg_diam_blocks / np.power(avg_spacing_blocks, 2))
                sigma_d_blocks = (1 / 2) * drag_cube * np.power(beta, 2) * (avg_submerged_height_blocks * avg_diam_blocks / avg_spacing_blocks_squared)
            drag_stress_array[cell] = sigma_d_blocks
            adjusted_shear_stress = pre_adjusted_tau / (1 + sigma_d_blocks)
        return adjusted_shear_stress
        
    def track_new_blocks(self, size_and_position_array,
                         channel_nodes):
        """
        track_new_blocks adds blocks delivered from the hillslopes to
        the tracking matrix. Blocks start in the middle of whichever node 
        they're delivered to.
        
        INPUTS
        ------
        1) size_and_position_array: 
            -column 0: node to which block is delivered
            -column 1: block side length
        
        2) tracking_mat
        """
        num_new_blocks = len(size_and_position_array[:, 0])
        
        for new in range(0, num_new_blocks):
            self.for_slicing = 0
            #self.blocks += 1
            try:
                next_entry = max(max(np.where(np.isfinite(self.tracking_mat[:, 0])))) + 1 #adding one for next open entry #!
            except ValueError:
                next_entry = 0
            self.for_slicing = next_entry + 1 #b/c when you slice it takes the one before the end                
            if self.for_slicing >= self.tracking_mat.shape[0]:
                addition = np.zeros((10000, 4), dtype = np.float64) #!
                addition[:, :] = np.nan
                self.tracking_mat = np.concatenate((self.tracking_mat, addition))
                print 'ADDED ELEMENTS TO TRACKING MATRIX'
            #else:
            #    pass
            #self.tracking_mat[next_entry, 0] = 0 #piece is a block
            self.tracking_mat[next_entry, 0] = size_and_position_array[new, 0] #node! dropped in at center of cell
            self.tracking_mat[next_entry, 1] = size_and_position_array[new, 1] #side length
            self.tracking_mat[next_entry, 2] = np.power(self.tracking_mat[next_entry, 1], 3) #calculate volume
            if self.tracking_mat[next_entry, 0] not in channel_nodes:
                sys.exit('track_new_blocks not in channel')
            if new == 0:
                print 'in-channel'
                print self.tracking_mat[next_entry, 1]
                print self.tracking_mat[next_entry, 0]                             
                print '-----------'   
    def calc_flow_depth_and_velocity(self,
                                     node, 
                                     slope, 
                                     is_block_in_cell, 
                                     z0, 
                                     q,
                                     flow_depth_array,
                                     flow_velocity_array,
                                     raw_shear_stress_array,
                                     adjusted_shear_stress_array,
                                     dx,
                                     drag_cube,
                                     drag_stress_array,
                                     chan_width):
        """
        calc_flow_depth_and_velocity calculates flow depth, mean velocity,
        raw shear stress, and block-adjusted shear stress.
        
        calls roughness_bisection and calc_shear_stress_with_roughness
        
        INPUTS
        ------
        node: the number of the node currently being operated on
        slope: the topographic steepest slope at that node
        is_block_in_cell: boolean array of whether a given block is in the cell
        z0: bed roughness length scale [m]
        q: water discharge per unit width [m2/s]
        flow_depth_array: water depth array at all nodes [m]
        flow_velocity_array: water velocity array at all nodes [m/s]
        raw_shear_stress_array: array of raw shear stress at every node
        adjusted_shear_stress_array: array of adjusted shear stress at every node
        dx: grid node spacing
        drag_cube: drag coefficient for a cube
        for_slicing: index for slicing tracking matrix
        """
        g = 9.81 #m/s2
        dens_water = 1000 #kg/m3
        h = self.roughness_bisection(q, slope, z0, g)  
        if np.isnan(np.sum(h)):
            print h
            sys.exit("NAN FOUND IN FLOW DEPTH-- KILLING MODEL")
        flow_depth_array[node] = h
        if h == 0:
            tau_initial = 0
            v = 0
        else:
            v = q / h
        flow_velocity_array[node] = v
        tau_initial = dens_water * g * h * slope #shear stress at each node(Pa)
        raw_shear_stress_array[node] = tau_initial  
        if np.isnan(np.sum(raw_shear_stress_array)):
            print raw_shear_stress_array
            sys.exit("NAN FOUND IN UNCORRECTED TAU-- KILLING MODEL")                     
        tau = self.calc_shear_stress_with_roughness(tau_initial,  
                                                    is_block_in_cell, 
                                                    drag_cube, 
                                                    h, dx, z0, 
                                                    node,
                                                    drag_stress_array,
                                                    chan_width)
        adjusted_shear_stress_array[node] = tau
        if np.isnan(np.sum(adjusted_shear_stress_array)):
            print adjusted_shear_stress_array
            sys.exit("NAN FOUND IN CORRECTED TAU-- KILLING MODEL")
        #print 'FLOW DEPTH ', h
        #return (h, v, tau)
        
    def erode_bed(self, stack, z, 
                  bedrock__critical_stress, 
                  cover_fraction_array,
                  drag_stress_array,
                  flow_depth_array,
                  bedrock__erodibility,
                  timestep,
                  downslope_link_length_array,
                  adjusted_shear_stress_array,
                  baselevel_node):
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
        g = 9.81 #m/s2
        dens_water = 1000 #kg/m3
        excess_tau = adjusted_shear_stress_array - bedrock__critical_stress
        #TAKEN OUT FOR LANDLAB IMPLEMENTATION self.surface_elev_array[-1] -= (baselevel_drop * self.timestep) #adjust baselevel node     
        stack_fsc = np.insert(stack, 0, baselevel_node)        
        for node in range(1, len(stack_fsc)):
            if excess_tau[stack_fsc[node]] <= 0:
                pass #new elev is the same as old b/c no erosion
            else:
                f_open = 1 - cover_fraction_array[stack_fsc[node]]  
                f_open = f_open.clip(min=0)
                #print 'z before: ', z[stack_fsc[node]]
                z[stack_fsc[node]] = (z[stack_fsc[node]] + \
                    (z[stack_fsc[node - 1]] * f_open * bedrock__erodibility * \
                    dens_water * g * flow_depth_array[stack_fsc[node]] * timestep / \
                    (downslope_link_length_array[stack_fsc[node]] * (1 + drag_stress_array[stack_fsc[node]]))) + (timestep * f_open * \
                     bedrock__erodibility * bedrock__critical_stress)) / (1 + (f_open * bedrock__erodibility * \
                     dens_water * g * flow_depth_array[stack_fsc[node]] * timestep / \
                     (downslope_link_length_array[stack_fsc[node]] * (1 + drag_stress_array[stack_fsc[node]]))))
                #print 'z after: ', z[stack_fsc[node]]
        if np.isnan(np.sum(z)):
            print z
            sys.exit("NAN FOUND IN ELEV ARRAY-- KILLING MODEL")

        #euler forward
        #excess_adjusted = adjusted_shear_stress_array[stack] - bedrock__critical_stress
        #excess_adjusted = excess_adjusted.clip(min=0)
        #f_open = 1 - cover_fraction_array[stack]
        #z[stack] -= f_open * bedrock__erodibility * excess_adjusted * timestep
        
        
            
    def calc_force_balance(self, is_block_in_cell, h, v,
                           drag_cube, timestep, slope, receiver, channel_nodes,
                           baselevel_node, current_node):
        """
        calc_force_balance calculates block motion by a force balance following
        Lamb et al (2015).
        
        INPUTS
        ------
        is_block_in_cell: boolean array indicating which blocks are in the cell
        h: mean flow depth in the node being worked on [m]
        v: mean flow velocity in the node being worked on [m/s]
        for_slicing: index for slicing block tracking matrix
        drag_cube: cube drag coefficient [-]
        timestep: model timestep [yrs]
        slope: steepest descent slope from the node being worked on [-]
        receiver: node immediately downstream of the node being worked on
        """
        #instantiate constants and temporary arrays
        g = 9.81 #m/s2
        dens_water = 1000 #kg/m3
        dens_sediment = 2650 #kg/m3
        #lam_hop = (1. / 365 / 24 / 3600) #average rate of motion = 1 grain diameter per year   
        drag_force = np.zeros(len(self.tracking_mat[0:self.for_slicing, 0]), dtype = np.float64)
        weight_force = np.zeros(len(self.tracking_mat[0:self.for_slicing, 0]), dtype = np.float64)
        shear_stress_force = np.zeros(len(self.tracking_mat[0:self.for_slicing, 0]), dtype = np.float64)
        #hop_length = np.zeros((len(self.tracking_mat[0:for_slicing, 0])))
        
        #motion of blocks
        submerged_block = (is_block_in_cell) & (self.tracking_mat[0:self.for_slicing, 1] < h)
        emergent_block = (is_block_in_cell) & (self.tracking_mat[0:self.for_slicing, 1] >= h)
        self.tracking_mat[0:self.for_slicing, 3][submerged_block] = self.tracking_mat[0:self.for_slicing, 1][submerged_block]
        self.tracking_mat[0:self.for_slicing, 3][emergent_block] = h
        
        #FORCE BALANCE
        friction_angle = np.radians(20)
        drag_force[is_block_in_cell] = (1. / 2) * drag_cube * \
            self.tracking_mat[0:self.for_slicing, 3][is_block_in_cell] * \
            self.tracking_mat[0:self.for_slicing, 1][is_block_in_cell] * \
            dens_water * np.power(v, 2)
        weight_force[is_block_in_cell] = g * ((dens_sediment - \
            dens_water) * (self.tracking_mat[0:self.for_slicing, 3][is_block_in_cell] * \
            np.power(self.tracking_mat[0:self.for_slicing, 1][is_block_in_cell], 2)) + \
            dens_sediment * (self.tracking_mat[0:self.for_slicing, 1][is_block_in_cell] - \
            self.tracking_mat[0:self.for_slicing, 3][is_block_in_cell])\
            * np.power(self.tracking_mat[0:self.for_slicing, 1][is_block_in_cell], 2))
        shear_stress_force[emergent_block] = 0         
        shear_stress_force[submerged_block] = (dens_water * g * \
            slope * (h - self.tracking_mat[0:self.for_slicing, 3][submerged_block])) * \
            np.power(self.tracking_mat[0:self.for_slicing, 1][submerged_block], 2)      
        lift_force = 0.85 * (drag_force + shear_stress_force)
        slope_rad = np.arctan(slope) 
        total_motion_force = weight_force * np.sin(slope_rad) + shear_stress_force + drag_force
        total_resist_force = (weight_force * np.cos(slope_rad) - lift_force) * np.tan(friction_angle)

        is_moving_block = (total_motion_force > total_resist_force) & (is_block_in_cell)
        #if sum(is_moving_block) > 0:
            #print '# MOVING: ', sum(is_moving_block)
            #print is_block_in_cell
        #print 'is_moving: ', is_moving_block
        #print 'motion force: ', total_motion_force
        #print 'resist force: ', total_resist_force
        #print 'is_block: ', is_block_in_cell
        
        #not_moving_block = (total_motion_force <= total_resist_force) & (is_block_in_cell)
        current_chan_array = np.where(channel_nodes == current_node)[0][0]
        if current_chan_array == 0:
            next_chan_node = baselevel_node
        else:
            next_chan_node = channel_nodes[current_chan_array - 1]
        self.tracking_mat[0:self.for_slicing, 0][is_moving_block] = next_chan_node
        #if (sum(is_moving_block) > 0) & (current_node == 198):
        #    print '\a'
        #    print channel_nodes
        #    print current_node
        #    print current_chan_array
        #    print next_chan_node
        #    print receiver
        #    print slope
        #    print is_moving_block
        #hop_length[is_moving_block] = np.random.poisson(lam_hop * timestep, \
        #    sum(1 for x in is_moving_block if x)) * self.tracking_mat[0:for_slicing, 1][is_moving_block]
        #hop_length[not_moving_block] = 0
        #self.tracking_mat[0:for_slicing, 4] = self.tracking_mat[0:for_slicing, 4] + hop_length #adds hop length to downstream distance #!
        #hop_length[:] = 0
        
    def erode_blocks(self, is_block_in_cell, block__erodibility,
                     block__critical_stress, timestep,
                     raw_shear_stress, adjusted_shear_stress,
                     block_wx_rate):
        """
        erode_blocks removes mass and volume from blocks in proportion to the
        shear stress exerted on them.
        
        INPUTS
        ------
        is_block_in_cell: boolean array describing whether blocks are in cell
        block__erodibility: erodibility constant "k" of blocks
        block__critical_stress: critical shear stress to erode blocks [Pa]
        for_slicing: slicing index for black tracking matrix
        timestep: model timestep [yr]
        raw_shear_stress: depth-slope product at node [Pa]
        adjusted_shear_stress: remaining bed shear stress at node [Pa]
        """
        
        excess_block_tau = (raw_shear_stress - adjusted_shear_stress - \
            block__critical_stress)
        #print 'excess tau on block: ', excess_block_tau
        #print 'is_in_cell: ', is_block_in_cell
        if excess_block_tau < 0:
            excess_block_tau = 0
        else:
            pass
        #abrade blocks in proportion into shear stress
        self.tracking_mat[0:self.for_slicing, 1][is_block_in_cell] -= block__erodibility \
            * excess_block_tau * timestep
        #weather blocks at same rate as on hillslope
        self.tracking_mat[0:self.for_slicing, 1][is_block_in_cell] -= block_wx_rate \
            * timestep
        
        self.tracking_mat[0:self.for_slicing, 1][is_block_in_cell] = self.tracking_mat[0:self.for_slicing, 1][is_block_in_cell].clip(min = 0) #no negative blocks
        self.tracking_mat[0:self.for_slicing, 2][is_block_in_cell] = np.power(self.tracking_mat[0:self.for_slicing, 1][is_block_in_cell], 3)  
            
    def calc_cover_frac(self, is_block_in_cell, x, node_spacing, chan_width,
                        cover_fraction_array):
        """
        calc_cover_frac calculates the cover fraction at a given node based
        on information about the number and size of blocks in the cell.
        
        MODIFIED TO ACCOUNT FOR INDEPENDENT DX AND CHAN WIDTH
        
        INPUTS
        ------
        is_block_in_cell: boolean array telling about block presence/absence
        x: number of node being worked on
        node_spacing: grid property of spacing between nodes [m]
        """
        block_cover = sum(np.power(self.tracking_mat[0:self.for_slicing, 1][is_block_in_cell], 2)) #area covered by blocks    #! 
        #cover_frac = block_cover / np.power(node_spacing, 2) #how much of cell is COVERED by big bits
        cover_frac = block_cover / (node_spacing * chan_width)
        cover_fraction_array[x] = 1 - np.exp(-cover_frac)
        #get number of blocks in cell for plotting
        #self.blocks_in_cells[x] = sum(is_block_in_cell)
        
    def delete_eroded_or_gone_blocks(self, baselevel_node, for_slicing):
        """
        delete_eroded_or_gone_blocks is a utility to keep the tracking matrix
        clean by deleting blocks that have eroded down to nothing.
        """
        #now find and delete any tracking rows where blocks have degraded to 0.
        for_abrasion_deletion = np.where(self.tracking_mat[0:for_slicing, 1] == 0)[0] #for_deletion is an array of row indices
        self.tracking_mat = np.delete(self.tracking_mat, for_abrasion_deletion, axis = 0) #deletes relevant rows
        
        #now find and delete blocks in the baselevel node (may help with memory?)
        for_transport_deletion = np.where(self.tracking_mat[0:for_slicing, 0] == baselevel_node)[0]
        self.tracking_mat = np.delete(self.tracking_mat, for_transport_deletion, axis = 0) #deletes relevant rows
        
        #I think that transport deletion isn't necessary because blocks will end up in
        #a boundary node that is being manually lowered anyway.        
        
        #now find and delete any tracking rows where blocks have left domain
        #for_transport_deletion = np.where(self.tracking_mat[0:for_slicing, 0] > self.downstream_edges_of_cell[-1])[0] #!
        #self.tracking_mat = np.delete(self.tracking_mat, for_transport_deletion, axis = 0)
    
#    def instantiate_common_constants(self):  
#        np.random.seed(50) #TURN THIS LINE ON TO ELIMINATE INHERENT VARIABILITY
#        self.dens_water = 1000 #kg/m^3
#        self.dens_sediment = 2650 #kg/m^3
#        self.g = 9.81 #acceleration due to gravity (m/s^2)
#        self.drag_cube = .8 #from Carling 1998
#        self.coeff_static_fric = .6 #from Carling 1998
        
#    def set_shear_stress_params(self, bed_k, block_k, tau_c, tau_c_block):
#        self.ke_br_bed = bed_k #bedrock is not easy to erode obviously
#        self.ke_br_block = block_k
#        self.a_br = 1 #bedrock
#        self.tau_c_br = tau_c #pascals. bedrock detachment critical shear stress
#        self.tau_c_block = tau_c_block
        
#    def set_roughness_depth_calc_params(self, z_0):
#        #self.tolerance = .001 #1 mm
#        self.z0 = z_0 #roughness height NEED TO MAKE THIS DEPENDENT ON SED SIZE DISTRIBUTION
#        self.upper_test_lim = 10000 #upper q limit to test e.g. maximum possible flow height
#        self.lower_test_lim = 0 #lowest possible flow height
        
    def instantiate_tracking_matrix(self):    
        """
        Sets up block tracking matrix as a big matrix of NaNs.
        -Column 0: node at the center of the cell in which the block lives
        -Column 1: cube side length
        -Column 2: cube volume (side length ^3)
        -Column 3: to be used for emergence/submergence over flow
        ##-Column 4: along-channel distance within node (range is 0:(dx*sqrt(2)))
        
        Also instantiates cover fraction array
        """
        #if number_of_pieces > 0:
        #    if gamma == 0:
        #        self.tracking_mat = np.zeros((number_of_pieces, 4), dtype = np.float64)#!
        #    else:
        #        self.tracking_mat = np.zeros((number_of_pieces * 1000, 4), dtype = np.float64) #able to track 6 attributes #!
        #else:
        #    self.tracking_mat = np.zeros((100000, 4), dtype=np.float64) #!
        self.tracking_mat = np.zeros((10000, 4), dtype=np.float64)
        self.tracking_mat[:, :] = np.nan #THIS IS SO THAT I CAN DISTINGUISH UN-ENTERED VALUES FROM ACTUAL ZEROS
        self.for_slicing = 0
        #self.tracking_mat[0:number_of_pieces, 0] = starting_dist #starting distance #!
        #self.tracking_mat[boulders:boulders + cubes + 1, 1] = self.cube_side_length #!
        #self.tracking_mat[boulders:boulders + cubes + 1, 2] = np.power(self.tracking_mat[boulders:boulders + cubes + 1, 1], 3) #!
        #self.tracking_mat[boulders:boulders + cubes + 1, 3] = 0#np.nan #TO BE USED FOR SUBMERGENCE/EMERGENCE #!
        #self.cover_frac_array = np.zeros(self.n_cells, dtype = np.float64)
        #self.tracking_mat[0:number_of_pieces, 0] = np.random.random_integers(0, 99 * self.dx, sum(~np.isnan(self.tracking_mat[:, 0]))) #!
        
    def calculate_channel_incision_rate(self, old, new, dt, row):
        """
        Calculates channel incision rate by differencing channel elevations
        between timesteps
        """
        incision_rate = (old - new) / 1000
        self.e_rate_tracking_array[row, :] = incision_rate
       
        
    def save_out_channel_incision_rate_array(self, directory, name_prefix):
        """
        At the end of the model run, save the matrix containing channel
        incision rates from every node and timestep to a .npy binary.
        """
        np.save(directory + '/' + name_prefix + '_e_rate_record.npy', self.e_rate_tracking_array)