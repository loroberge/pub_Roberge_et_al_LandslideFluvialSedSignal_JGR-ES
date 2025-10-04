# -*- coding: utf-8 -*-
"""
Created on Mon Oct 17 16:41:56 2022

@author: LaurentRoberge
"""

import sys
import yaml
import time
import warnings
import numpy as np
from pathlib import Path

from landlab.io.netcdf import read_netcdf
from landlab import load_params
from landlab.components import (PriorityFloodFlowRouter,
                                SpaceLargeScaleEroder,
                                ExponentialWeatherer,
                                DepthDependentDiffuser,
                                BedrockLandslider,
                                SteepnessFinder
                                )

from library import (create_run_names,
                     create_empty_arrays,
                     create_runtime_arrays,
                     calc_runtimes,
                     export_run_time,
                     export_grids,
                     export_grids_final,
                     collect_stats_grid,
                     collect_stats_ls,
                     calculate_ls_soil_fraction,
                     combine_and_export_timeseries
                     )

# %% Filter warnings
warnings.filterwarnings("ignore")

# %% Set directory paths
input_file = sys.argv[1]
#input_file = 'input_EQ.txt'

path_project = Path(__file__).resolve().parent
path_input = path_project.joinpath(input_file)

# %% Load input parameters
params = load_params(str(path_input))

params_fr = params['params_fr']
params_sp = params['params_sp']
params_ew = params['params_ew']
params_ddd = params['params_ddd']
params_ls = params['params_ls']

# %% Generate model run names for import and export
runtype = "EQ"
run_name_import, run_name_export = create_run_names(params, params_sp, params_ls, runtype)

path_import = path_project.joinpath('Output', run_name_import) 
path_export = path_project.joinpath('Output', run_name_export) 
if not path_export.exists():
    path_export.mkdir(parents = True)

# %% Set up: grid, time, and other parameters
nrows = params['nrows']
ncols = params['ncols']
n_core_nodes = (nrows-2)*(ncols-2)
dx = params['dx']
dt = params['dt']
channel_threshold_area = params['channel_threshold_area']
theta = params_sp['m_sp']/params_sp['n_sp']

progress_interval = params['progress_interval']
output_interval = params['output_interval']
total_t = params['total_t'] + dt
ndt = int(total_t // dt)

uplift_per_step = params['uplift_rate'] * dt

outlet_node = 0
node_next_to_outlet = ncols+1

# %% Load model grid from NetCDF file

# Import topographic and bedrock elevation fields from NetCDF file
mg = read_netcdf(str(path_import.joinpath(run_name_import + '_Grid_FINAL.nc')))

# Add soil depth field
_ = mg.add_zeros('soil__depth', at='node', units= ['m','m'])
mg.at_node['soil__depth'][:] = (mg.at_node['topographic__elevation']
                                - mg.at_node["bedrock__elevation"])

# Add fields for cumulative channel age, LS erosion, and LS deposition
_ = mg.add_zeros('channel__time_since_occupation', at='node', units= ['m','m'])
mg.at_node['channel__time_since_occupation'][:] = np.nan
_ = mg.add_zeros('landslide__erosion_cumulative', at='node', units= ['m','m'])
_ = mg.add_zeros('landslide__deposition_cumulative', at='node', units= ['m','m'])

# Define grid units and set boundary conditions
mg.axis_units = ('m', 'm')
mg.set_status_at_node_on_edges(right=4,
                                top=4,
                                left=4,
                                bottom=4)
mg.status_at_node[outlet_node] = mg.BC_NODE_IS_FIXED_VALUE
    
# %% Instantiate model components
fr = PriorityFloodFlowRouter(mg, **params_fr)
fr.run_one_step()

sp = SpaceLargeScaleEroder(mg, **params_sp)

ew = ExponentialWeatherer(mg, **params_ew)

ddd = DepthDependentDiffuser(mg, **params_ddd)

ls = BedrockLandslider(mg, **params_ls)

sf = SteepnessFinder(mg,
                     reference_concavity=theta,
                     min_drainage_area=channel_threshold_area
                     )

#%% Create empty arrays to fill during loop
ts_data, riv_data, vol_data = create_empty_arrays(ndt, n_core_nodes)
t_start, t_end, runtimes = create_empty_arrays(ndt)

# Used for landscape memory correlation coefficient
topo_ref = mg.at_node["topographic__elevation"][mg.core_nodes]
topo_ref_mean = np.mean(topo_ref)
topo_ref_val = topo_ref - topo_ref_mean
topo_ref_val_sq = np.sqrt(np.sum(np.square(topo_ref - topo_ref_mean)))

# Landslide magnitudes
landslides_size_all_steps = []

# %% Model Run

# Set initial timestamp (used to print progress updates)
start_time = time.time()

# Set elapsed model time to 0 years
elapsed_time = 0

for i in range(ndt):
    # Update time counter
    elapsed_time = i*dt 
    
    # Add uplift
    mg.at_node['bedrock__elevation'][mg.core_nodes] += uplift_per_step
    
    # Run ExponentialWeatherer
    t_start['ew'] = time.time()
    ew.calc_soil_prod_rate()
    soil_prod_per_step = mg.at_node['soil_production__rate'][mg.core_nodes] * dt
    t_end['ew'] = time.time()
    
    # Convert bedrock to soil using soil production rate
    mg.at_node['bedrock__elevation'][mg.core_nodes] -= soil_prod_per_step
    mg.at_node['soil__depth'][mg.core_nodes] += soil_prod_per_step
    
    # Update topographic elevation to match bedrock and soil depths
    mg.at_node['topographic__elevation'][:] = (mg.at_node["bedrock__elevation"]
                                               + mg.at_node["soil__depth"])
    
    # Topography prior to fluvial incision
    topo_pre_sp = mg.at_node['topographic__elevation'][mg.core_nodes]
    
    # Mass conservation check: total topographic volume before fluvial incision
    vol_data["vol_topo_pre_sp"] = np.sum(mg.at_node['topographic__elevation'][mg.core_nodes]) * dx * dx
    vol_data["vol_br_pre_sp"] = np.sum(mg.at_node['bedrock__elevation'][mg.core_nodes]) * dx * dx
    vol_data["vol_sed_pre_sp"] = np.sum(mg.at_node['soil__depth'][mg.core_nodes]) * dx * dx
    
    # River channel nodes prior to running flow routing/accumulation
    riv_data["Nodes_r"] = mg.at_node['drainage_area'][mg.core_nodes] >= channel_threshold_area
    
    # Run PriorityFloodFlowRouter, SPACE, SteepnessFinder
    t_start['fr'] = time.time()
    fr.run_one_step()
    t_end['fr'] = time.time()
    
    t_start['sp'] = time.time()
    sp.run_one_step(dt=dt)
    t_end['sp'] = time.time()
    
    t_start['sf'] = time.time()
    sf.calculate_steepnesses()
    t_end['sf'] = time.time()
        
    # River channel nodes after running flow routing/accumulation
    riv_data["Nodes_r2"] = mg.at_node['drainage_area'][mg.core_nodes] >= channel_threshold_area
    
    # Mass conservation check: total topographic volume after fluvial incision
    vol_data["vol_topo_post_sp"] = np.sum(mg.at_node['topographic__elevation'][mg.core_nodes]) * dx * dx
    vol_data["vol_br_post_sp"] = np.sum(mg.at_node['bedrock__elevation'][mg.core_nodes]) * dx * dx
    vol_data["vol_sed_post_sp"] = np.sum(mg.at_node['soil__depth'][mg.core_nodes]) * dx * dx
    
    # Run DepthDependentDiffuser, BedrockLandslider
    t_start['ddd'] = time.time()
    ddd.run_one_step(dt=dt)
    t_end['ddd'] = time.time()
    
    t_start['ls'] = time.time()
    ls.run_one_step(dt=dt)
    t_end['ls'] = time.time()
    
    # Topography after landsliding
    topo_post_ls = mg.at_node['topographic__elevation'][mg.core_nodes]
    
    # Mass conservation check: total topographic volume after landsliding
    vol_data["vol_topo_post_ls"] = np.sum(mg.at_node['topographic__elevation'][mg.core_nodes]) * dx * dx
    vol_data["vol_br_post_ls"] = np.sum(mg.at_node['bedrock__elevation'][mg.core_nodes]) * dx * dx
    vol_data["vol_sed_post_ls"] = np.sum(mg.at_node['soil__depth'][mg.core_nodes]) * dx * dx
    
    # Update cumulative channel age field
    mg.at_node['channel__time_since_occupation'][np.where(riv_data["Nodes_r"]==1)] = 0
    mg.at_node['channel__time_since_occupation'][np.where(riv_data["Nodes_r"]==0) and
                               ~np.isnan(mg.at_node['channel__time_since_occupation'][:])] += dt
    
    # Update cumulative LS erosion and deposition fields
    mg.at_node['landslide__erosion_cumulative'][mg.core_nodes] += mg.at_node['landslide__erosion'][mg.core_nodes]
    mg.at_node['landslide__deposition_cumulative'][mg.core_nodes] += mg.at_node['landslide__deposition'][mg.core_nodes]
    
    collect_stats_grid(ts_data, vol_data, riv_data, elapsed_time, i, dt, mg, 
                       uplift_per_step, node_next_to_outlet, n_core_nodes, 
                       channel_threshold_area, topo_pre_sp, topo_post_ls, 
                       topo_ref_val, topo_ref_val_sq
                       )
    collect_stats_ls(ts_data, vol_data, i, dt, dx, ls)
    calc_runtimes(t_start, t_end, runtimes, i)
    
    # Store landslide sizes from current time step into general ls_size list
    landslides_size_all_steps = np.append(landslides_size_all_steps, ls.landslides_size)
    
    # Print progress and output grids
    if np.mod(elapsed_time, progress_interval) == 0:
        print(elapsed_time/1000," ka",("-- %s seconds --" % round((time.time() - start_time),1)))
        print("%.1f percent complete" % (elapsed_time*100 / total_t))
        print("Max elev: ", round(max(mg.at_node['topographic__elevation'][mg.core_nodes]),1))
    if np.mod(elapsed_time, output_interval) == 0:
        export_grids(mg,
                    params['grid_fields_to_save'],
                    path_export,
                    run_name_export, 
                    elapsed_time
                    )
        mg.at_node['landslide__erosion_cumulative'][mg.core_nodes] *= 0.
        mg.at_node['landslide__deposition_cumulative'][mg.core_nodes] *= 0.
                    

# Export grid fields from final timestep
export_grids_final(mg,
                   params['grid_fields_to_save'],
                   path_export,
                   run_name_export, 
                   elapsed_time
                   )

# Calculate fraction of soil in landslides
calculate_ls_soil_fraction(ts_data)

# Combine all timeseries data into one dataframe and export as csv.
combine_and_export_timeseries(ts_data, path_export, run_name_export)

# Export record of all landslide sizes (for frequency-magnitude analysis)
landslides_size_all_steps *= (dx**2)
np.savetxt(str(path_export.joinpath(run_name_export + '_DataLSsizes.csv')), landslides_size_all_steps, delimiter=",") 

# Export record of parameters used in model run
with open(str(path_export.joinpath(run_name_export + '_Parameters.txt')), 'w') as file:
     file.write(yaml.dump(params, indent=4, sort_keys=True))

# Print and export total model run time
end_time = time.time()
export_run_time(start_time, end_time, path_export)

