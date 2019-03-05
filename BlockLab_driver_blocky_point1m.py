# -*- coding: utf-8 -*-
"""
Driver script for BlockLab blocky landscape evolution model. This driver calls 
functions from Landlab (Hobley et al., 2017), BRaKE (Shobe et al., 2016; 2018),
and HogLab (Glade et al., 2017; 2018).

Authors: Rachel Glade and Charles Shobe, University of Colorado, Boulder

Language: Python 2.7 with Numpy and Cython

Development timeline:
    -Coupling between Landlab, BRaKE, and HogLab begin January 2018
    -Model coupling and testing completed June 2018
    -Archived on GitHub October 2018
    -Revised and re-released on GitHub March 2019
    
"""
#import all required model functionality
from __future__ import division
import os
import numpy as np
import brakeLandlab as bl
import hogLab as hl
from landlab import RasterModelGrid
from landlab.components import FlowDirectorD8
from landlab.components import DepthDependentDiffuser, ExponentialWeatherer
np.set_printoptions(threshold = np.nan)
from landlab.io.netcdf import write_netcdf
from cfuncs import erode_bed, move_blocks
from landlab import CLOSED_BOUNDARY

######INITIALIZE: MODEL SETUP AND PARAMETERS##################################
#name prefix for model run data
name_prefix = input("Enter run name, use quotes: ")

#directory to save model output
directory = 'runs/' + name_prefix
if not os.path.exists(directory):
    os.makedirs(directory)

#define run params
time_to_run = 500000 #years
timestep = 2 #years
plot_every = 1000 #years
nplots = time_to_run/plot_every

#setup grid
dx = 5.0
channel_width = 10 #m constant
nr = 200
nc = 400
save_topo_every = 1000 #save topo every ____ years
save_topo_flag = True #True to save netcdf every save_topo_every years.
save_profiles_every = 10000 #save out channel profile arrays every _____ years
save_profiles_flag = True #True to save info every dt, written out every save_profiles_every years
mg = RasterModelGrid((nr, nc), dx)

#define fluvial model parameters
bl_drop = 0.0001 #m/yr, baselevel lowering rate
z0 = 0.1 #m; bed roughness length scale
bedrock__erodibility = 0.000001
bedrock__critical_stress = 0. #critical shear stress to erode bedrock [Pa]
block__erodibility = 0.00000005
block__critical_stress = 0 #critical shear stress to erode blocks [Pa]
drag_cube = 0.8 #drag coefficient on a cube [-]
constant_volumetric_Q = 30. #m3/s
constant_q = constant_volumetric_Q / channel_width

#define hillslope parameters
linear_diffusivity = 0.2
soil_transport_decay_depth = 0.1
max_soil_production_rate = 1.e-3
soil_production_decay_depth = 0.2

#define hillslope nodes for which to save erosion rate time series
downstream_hillslope_node = 1 #node number along channel (bottom to top)
middle_hillslope_node = 60
upstream_hillslope_node = 120

#define parameters governing resistant rock layer
thickness = .1 #number of blocks thick
block_height = thickness
nrocks = 12
dip = 20
max_block_weathering_rate = 5e-5 
block_channel_weathering_rate = 0
release_relief_threshold = thickness + 1 #cliff release relief threshold
block_motion_threshold = 2 #must be multiple of block_height

bedrock_elevation = mg.add_zeros('node', 'bedrock__elevation')

z = mg.add_zeros('node', 'topographic__elevation')
bedrock_elevation[:] = z

#set boundary conditions
mg.set_closed_boundaries_at_grid_edges(bottom_is_closed = True,
                                       left_is_closed = True,
                                       right_is_closed = True,
                                       top_is_closed = True)

mg.set_watershed_boundary_condition_outlet_id(mg.nodes_at_right_edge[0] - int((nc/2) + 0.5), 
    mg['node']['topographic__elevation'], -9999) #open middle node

#instantiate hoglab
hog = hl.HogLab()

#grid fields
soil_production_rate = mg.add_zeros('node','soil_production__rate')
top_of_hard_layer = mg.add_zeros('node','zb_hard_top')
bottom_of_hard_layer = mg.add_zeros('node','zb_hard_bottom')
slope = mg.add_zeros('node','topographic__steepest_slope')

#add the fields required to do blocks
flow_depth_array = mg.add_zeros('node', 'flow_depth')
flow_velocity_array = mg.add_zeros('node', 'flow_velocity')
raw_shear_stress_array = mg.add_zeros('node', 'raw_shear_stress')
adjusted_shear_stress_array = mg.add_zeros('node', 'adjusted_shear_stress')
drag_stress_array = mg.add_zeros('node', 'drag_stress')
cover_fraction_array = mg.add_zeros('node', 'cover_fraction')
block_sizes = mg.add_zeros('node','block_sizes')

#fields for hillslopes
soil_depth = mg.add_zeros('node', 'soil__depth')
underlying_soil_depth = mg.add_zeros('node', 'underlying_soil_depth')
    
initial_elevation = 100 #initial elevation
mg.at_node['bedrock__elevation'][:] = initial_elevation

#give channel a slope towards baselevel #operate on BR elev b/c summation is after
#define channel nodes
channel_nodes = mg.nodes_at_right_edge - int((nc/2) + 0.5)
baselevel_node = mg.nodes_at_right_edge[0] - int((nc/2) + 0.5) #this is passed into brake.erode_bed

mg.at_node['bedrock__elevation'][channel_nodes] = initial_elevation - 0.001 * (max(mg.node_y[channel_nodes]) - mg.node_y[channel_nodes])

mg.at_node['bedrock__elevation'][baselevel_node] = min(mg.at_node['bedrock__elevation'][channel_nodes]) - 0.0001

channel_nodes = channel_nodes[1:-1]

#sum bedrock and soil to get topo elev
mg.at_node['topographic__elevation'][:] = mg.at_node['bedrock__elevation'] + mg.at_node['soil__depth'] 

#set second boundary layer to keep hillslope from going into baselevel node
second_boundary_layer = mg.nodes_at_bottom_edge + nc
boundary_remove = np.where(second_boundary_layer == channel_nodes[0])
second_boundary_layer = np.delete(second_boundary_layer,boundary_remove)
mg.status_at_node[second_boundary_layer] = CLOSED_BOUNDARY

#inititalize resistant layer
starting_point = max(mg.node_x) - (2 * dx)

#route flow
fr = FlowDirectorD8(mg,'topographic__elevation')

#instantiate BRaKE, the collection of functions governing river evolution
blocks = bl.BrakeLandlab() 

#instantiate hillslope diffusion model component, from Landlab
hillslopeflux = DepthDependentDiffuser(mg, linear_diffusivity = linear_diffusivity, 
                                       soil_transport_decay_depth = soil_transport_decay_depth)

#instantiate exponential weathering component, from Landlab
weatheringrate = ExponentialWeatherer(mg, soil_production__maximum_rate = max_soil_production_rate, 
                                      soil_production__decay_depth = soil_production_decay_depth)

#set boundary elevations to enable easier data visualization in Paraview
mg.at_node['topographic__elevation'][mg.boundary_nodes] -= (bl_drop * time_to_run)
mg.at_node['bedrock__elevation'][mg.boundary_nodes] -= (bl_drop * time_to_run)
mg.at_node['bedrock__elevation'][baselevel_node] = min(mg.at_node['bedrock__elevation'][channel_nodes]) - 0.0001
mg.at_node['topographic__elevation'][baselevel_node] = mg.at_node['bedrock__elevation'][baselevel_node]

#instantiate fluvial block tracking matrix
blocks.instantiate_tracking_matrix()

#instantiate channel elev and cover frac arrays for saving
#these are arrays that will get saved every xxxx years to conserve memory
channel_node_elevs_for_saving = np.zeros((int(save_profiles_every / timestep), int(len(channel_nodes) + 1)))
channel_node_cover_fracs_for_saving = np.zeros((int(save_profiles_every / timestep), int(len(channel_nodes) + 1)))

#instantiate arrays for tracking erosion rates at 3 points on the hillslope
downstream_hillslope = channel_nodes[downstream_hillslope_node] - 1
middle_hillslope = channel_nodes[middle_hillslope_node] - 1
upstream_hillslope = channel_nodes[upstream_hillslope_node] - 1
downstream_channel_adjacent_node = np.zeros(int(time_to_run/timestep))
middle_channel_adjacent_node = np.zeros(int(time_to_run/timestep))
upstream_channel_adjacent_node = np.zeros(int(time_to_run/timestep))

hog.setup_block_tracking_arrays() 

slope_threshold, slope_rad = hog.setup_horizontal_layer(thickness,
                                                           block_height,
                                                           nrocks,
                                                           release_relief_threshold,
                                                           starting_point,
                                                           initial_elevation,
                                                           top_of_hard_layer,
                                                           bottom_of_hard_layer,
                                                           mg.node_x,
                                                           bedrock_elevation,
                                                           mg.node_spacing,
                                                           mg.core_nodes,
                                                           channel_nodes)

#route flow 
fr.run_one_step()

###########RUN: TIME LOOP#####################################################
#time loop
elapsed_time = 0
loop_iter = 0
profile_save_counter = 0
hillslope_save_counter = 0

while elapsed_time < time_to_run:
    if elapsed_time == 0:
        #parameter dictionary
        param_dict = {'dx' : dx,
         'nr' : nr,
         'nc' : nc,
         'thickness' : thickness,
         'block_height' : block_height,
         'dip' : dip,
         'max_block_weathering_rate' : max_block_weathering_rate,
         'release_relief_threshold' : release_relief_threshold,
         'block_motion_threshold' : block_motion_threshold,
         'time_to_run' : time_to_run,
         'timestep' : timestep,
         'bl_drop' : bl_drop,
         'z0' : z0,
         'bedrock__erodibility' : bedrock__erodibility,
         'bedrock__critical_stress' : bedrock__critical_stress,
         'block__erodibility' : block__erodibility,
         'block__critical_stress' : block__critical_stress,
         'drag_cube' : drag_cube,
         'constant_q' : constant_q,
        }
        
        #save param dictionary to npy binary 
        np.save(directory + '/' + name_prefix + '_param_dict.npy', param_dict) #load with "params = np.load('param_dict.py').item()"
    
    drop_locations = hog.find_release_location(bedrock_elevation,
                                               top_of_hard_layer,
                                               bottom_of_hard_layer,
                                               slope,
                                               slope_threshold,
                                               channel_nodes)    
    hog.release_blocks(z,
                    bedrock_elevation,
                    slope_threshold,
                    mg.core_nodes,
                    underlying_soil_depth,
                    soil_depth,
                    slope_rad,
                    drop_locations,
                    block_height,
                    mg.at_node['flow__receiver_node'], channel_nodes, slope, release_relief_threshold)
                    
    
    hog.block_locations, hog.block_sizes, hog.cum_channel_blocks, hog.blocks_to_channel_size, hog.blocks_to_channel_location = move_blocks(hog.block_locations,
                                             hog.block_sizes,
                                             mg.at_node['topographic__steepest_slope'], 
                                             block_motion_threshold,
                                             mg.node_spacing,
                                             soil_depth,
                                             underlying_soil_depth,
                                             bedrock_elevation,
                                             z,
                                             mg.at_node['flow__receiver_node'], 
                                             mg.core_nodes,
                                             channel_nodes,
                                             hog.hillslope_nodes,
                                             hog.blocks_to_channel_size,
                                             hog.blocks_to_channel_location,
                                             hog.cum_channel_blocks)
                                               
    #route flow 
    fr.run_one_step()
                                                                            
    one_block_volume = hog.blocks_to_channel_size * np.power(mg.node_spacing, 2) / nrocks    
    one_block_side_length = np.power(one_block_volume, 1. / 3)                                                                    
    side_length_times_nrocks = np.repeat(one_block_side_length, nrocks)
    block_location_times_nrocks = np.repeat(hog.blocks_to_channel_location, nrocks)                                                   
    delivery_size_and_position_array = np.column_stack((block_location_times_nrocks, side_length_times_nrocks))
    
    hog.blocks_to_channel_location = np.zeros(shape = 0, dtype = np.int64)
    hog.blocks_to_channel_size = np.zeros(shape = 0, dtype = np.float64)
    
    #can do this outside of spatial loop b/c it's just getting blocks into tracking matrix
    blocks.track_new_blocks(delivery_size_and_position_array, 
                                          channel_nodes) 
    
    stack = channel_nodes #get stack
    
    #inner loop iterates across all channel nodes 
    for x in stack: #np.flipud(stack):    
        if mg.status_at_node[x] == 0:          
            #define which blocks are in this node:
            is_block_in_cell = blocks.tracking_mat[0:blocks.for_slicing, 0] == x 
            
            #calculate flow depth and velocity from discharge (calls roughness solver)
            #this function also calculates basic and modified shear stresses
            blocks.calc_flow_depth_and_velocity(x,
                                                mg['node']['topographic__steepest_slope'][x],
                                                is_block_in_cell,
                                                z0, 
                                                constant_q,
                                                flow_depth_array,
                                                flow_velocity_array,
                                                raw_shear_stress_array,
                                                adjusted_shear_stress_array,
                                                mg.node_spacing,
                                                drag_cube,
                                                drag_stress_array,
                                                channel_width) 
        
            #calculate the cover fraction
            blocks.calc_cover_frac(is_block_in_cell, x, mg.node_spacing,
                                   channel_width,
                                   cover_fraction_array) 
        else:
            pass
        
    #now out of the 'for' loop to erode the channel bed 
    downslope_link_length_array = mg.length_of_d8[mg.at_node['flow__link_to_receiver_node']]
    erode_bed(stack, z, 
                     bedrock__critical_stress, 
                     cover_fraction_array,
                     drag_stress_array,
                     flow_depth_array,
                     bedrock__erodibility,
                     timestep,
                     downslope_link_length_array,
                     adjusted_shear_stress_array,
                     baselevel_node) 
                     
    #need bedrock erosion to equal new topographic elevation
    mg.at_node['bedrock__elevation'][channel_nodes] = mg.at_node['topographic__elevation'][channel_nodes]
        
    #now back into a second 'for' loop to do block transport and degradation
    for x in stack:
        if mg.status_at_node[x] == 0:
            receiver = mg.at_node['flow__receiver_node'][x] #receiver of current node
            #calculate force balance for block motion
            is_block_in_cell = blocks.tracking_mat[0:blocks.for_slicing, 0] == x 
            blocks.calc_force_balance(is_block_in_cell, 
                                      mg.at_node['flow_depth'][x], 
                                      mg.at_node['flow_velocity'][x], 
                                      drag_cube, timestep, 
                                      mg.at_node['topographic__steepest_slope'][x], 
                                      receiver,
                                      channel_nodes,
                                      baselevel_node, x) 
            
            #re-do is_block_in_cell after block motion
            is_block_in_cell = blocks.tracking_mat[0:blocks.for_slicing, 0] == x 
            #is_block_in_cell = is_block_in_cell[0]
            
            #erode blocks
            blocks.erode_blocks(is_block_in_cell, block__erodibility,
                                block__critical_stress, timestep,
                                raw_shear_stress_array[x], 
                                adjusted_shear_stress_array[x],
                                block_channel_weathering_rate)
            
            #delete blocks that have passed out of domain or eroded down to nothing
            blocks.delete_eroded_or_gone_blocks(baselevel_node, blocks.for_slicing) 
        else:
            pass
        
    #now, the hillslope evolution processes
    
    #weathering
    weatheringrate.calc_soil_prod_rate()
                                                                        
    hog.set_soil_production_rate_on_hard_layer(soil_production_rate)
    
    mg.at_node['soil_production__rate'][channel_nodes] = 0
    
    hog.weather_blocks(soil_production_rate,
                       max_block_weathering_rate,
                       underlying_soil_depth,
                       soil_depth,
                       timestep)    
    
    #flux
    hillslopeflux.soilflux(timestep)
    
    mg.at_node['topographic__elevation'][channel_nodes] = mg.at_node['bedrock__elevation'][channel_nodes]
    mg.at_node['soil__depth'][channel_nodes] = 0
    
    mg.at_node['bedrock__elevation'][baselevel_node] -= bl_drop * timestep
    
    #sum bedrock and soil to get topographic elevation
    mg.at_node['topographic__elevation'][:] = mg.at_node['bedrock__elevation'] + mg.at_node['soil__depth']

    #get block locations
    block_sizes[:] = 0
    block_sizes[hog.block_locations] = hog.block_sizes 
    
    #calc and save hillslope topography at 3 points every timestep
    downstream_channel_adjacent_node[hillslope_save_counter] = mg.at_node['topographic__elevation'][downstream_hillslope] - block_sizes[downstream_hillslope]
    middle_channel_adjacent_node[hillslope_save_counter] = mg.at_node['topographic__elevation'][middle_hillslope] - block_sizes[middle_hillslope]
    upstream_channel_adjacent_node[hillslope_save_counter] = mg.at_node['topographic__elevation'][upstream_hillslope] - block_sizes[upstream_hillslope]

    hillslope_save_counter += 1
    
    #populate channel node elevation array, which gets saved every _____ years
    temp_channel = mg.at_node['topographic__elevation'][channel_nodes]
    temp_channel = np.insert(temp_channel, 0, mg.at_node['topographic__elevation'][baselevel_node])
    channel_node_elevs_for_saving[profile_save_counter, :] = temp_channel
    
    #populate channel cover fraction array, which gets saved every _____ years
    temp_cover = mg.at_node['cover_fraction'][channel_nodes]
    temp_cover = np.insert(temp_cover, 0, mg.at_node['cover_fraction'][baselevel_node])
    channel_node_cover_fracs_for_saving[profile_save_counter, :] = temp_cover
    
    profile_save_counter += 1
    
    ##FINALIZE: SAVE DATA#####################################################
    
    elapsed_time += timestep
    
    if (save_profiles_flag == True) & (np.isclose(elapsed_time % save_profiles_every, 0.)):
        #save channel elevations over last save_profiles_every years
        np.save(directory + '/' + name_prefix + '_chan_elevs' + str(elapsed_time) + '.npy', channel_node_elevs_for_saving)
        channel_node_elevs_for_saving[:] = 0
        
        #save channel cover fractions over last save_profiles_every years
        np.save(directory + '/' + name_prefix + '_cover_fracs' + str(elapsed_time) + '.npy', channel_node_cover_fracs_for_saving)
        channel_node_cover_fracs_for_saving[:] = 0
        
        #reset profile_save_counter
        profile_save_counter = 0
        
    if np.isclose(elapsed_time % save_topo_every, 0.): 
        
        #save topography as netcdf files
        if save_topo_flag == True:
            block_sizes[:] = 0
            block_sizes[hog.block_locations] = hog.block_sizes
            temp_chan_sizes = np.zeros(len(channel_nodes))
            for cn in range(len(channel_nodes)):
                is_block_in_cell = blocks.tracking_mat[0:blocks.for_slicing, 0] == channel_nodes[cn] 
                mean_block_size = np.mean(blocks.tracking_mat[0:blocks.for_slicing, 1][is_block_in_cell])
                hillslope_size = (np.power(mean_block_size, 3) * nrocks) / np.power(dx, 2)
                temp_chan_sizes[cn] = hillslope_size

            block_sizes[channel_nodes] = temp_chan_sizes[:]
            write_netcdf(directory + '/' + name_prefix + str(elapsed_time) + '.nc', mg, format='NETCDF3_64BIT', names=('topographic__elevation','block_sizes'))
            block_sizes[channel_nodes] = 0
        
        print elapsed_time
    
    loop_iter += 1
    
#once model is finished, save out erosion rate time series at chosen hillslope nodes
np.save(directory + '/' + name_prefix + '_downstream_hillslope.npy', downstream_channel_adjacent_node)
np.save(directory + '/' + name_prefix + '_middle_hillslope.npy', middle_channel_adjacent_node)
np.save(directory + '/' + name_prefix + '_upstream_hillslope.npy', upstream_channel_adjacent_node)