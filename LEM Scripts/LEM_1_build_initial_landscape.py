# -*- coding: utf-8 -*-
"""
Created on Tue Oct 18 11:57:09 2022

@author: LaurentRoberge
"""

import sys
import yaml
import time
import warnings
import numpy as np
from pathlib import Path

from landlab import RasterModelGrid, load_params
from landlab.components import (PriorityFloodFlowRouter,
                                SpaceLargeScaleEroder,
                                ExponentialWeatherer,
                                DepthDependentDiffuser
                                )

from library import (create_run_names,
                     export_grids,
                     export_grids_final,
                     check_for_eq
                     )

# %% Filter warnings
warnings.filterwarnings('ignore')

# %% Set directory paths
input_file = sys.argv[1]
#input_file = 'input_INITIAL.txt'

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
runtype = "BASE"
run_name_import, run_name_export = create_run_names(params, params_sp, params_ls, runtype)

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

progress_interval = params['progress_interval']
output_interval = params['output_interval']
total_t = params['total_t'] + dt
ndt = int(total_t // dt)

T_eq_window = np.linspace(0,params['eq_interval'],params['eq_interval']//dt)
Zm_buffer = params['Zm_buffer']/params['eq_interval']
Zm_slope_lims = [0-Zm_buffer, 0+Zm_buffer]

uplift_per_step = params['uplift_rate'] * dt

outlet_node = 0
node_next_to_outlet = ncols+1

# %% Generate model grid from initial noise

mg = RasterModelGrid((nrows, ncols), dx)
mg.axis_units = ('m', 'm')
mg.set_status_at_node_on_edges(right=4,
                               top=4,
                               left=4,
                               bottom=4)
mg.status_at_node[outlet_node] = mg.BC_NODE_IS_FIXED_VALUE

# Soil depth
_ = mg.add_zeros('soil__depth', at='node', units= ['m','m'])
mg.at_node['soil__depth'] += params['soil_depth_initial']

# Bedrock elevation
np.random.seed(params['grid_seed'])
_ = mg.add_zeros('bedrock__elevation', at='node', units= ['m','m'])
mg.at_node['bedrock__elevation'] += np.random.rand(mg.number_of_nodes) / 10
mg.at_node['bedrock__elevation'][outlet_node] = 0

# Topographic elevation
_ = mg.add_zeros('topographic__elevation', at='node', units= ['m','m'])
mg.at_node['topographic__elevation'][:] += mg.at_node['bedrock__elevation']
mg.at_node['topographic__elevation'][:] += mg.at_node['soil__depth']
    
n_core_nodes = len(mg.core_nodes)
    
# %% Instantiate model components
fr = PriorityFloodFlowRouter(mg, **params_fr)
fr.run_one_step()

sp = SpaceLargeScaleEroder(mg, **params_sp)

ew = ExponentialWeatherer(mg, **params_ew)

ddd = DepthDependentDiffuser(mg, **params_ddd)

#%% Create empty arrays to fill during loop

T = np.zeros(ndt)       # Time 
Zm = np.zeros(ndt)  # Mean elevation timeseries to check for equilibrium

# %% Model Run

# Set elapsed model time to 0 years
elapsed_time = 0

# Set switch to end run when dynamic equilibrium is achieved
end_run = False

# Set initial timestamp (used to print progress updates)
start_time = time.time()

for i in range(ndt):
    # Update time counter
    elapsed_time = i*dt
    
    # Add uplift
    mg.at_node['bedrock__elevation'][mg.core_nodes] += uplift_per_step
    
    # Run ExponentialWeatherer
    ew.calc_soil_prod_rate()
    soil_prod_per_step = mg.at_node['soil_production__rate'][mg.core_nodes] * dt
    
    # Convert bedrock to soil using soil production rate
    mg.at_node['bedrock__elevation'][mg.core_nodes] -= soil_prod_per_step
    mg.at_node['soil__depth'][mg.core_nodes] += soil_prod_per_step
    
    # Update topographic elevation to match bedrock and soil depths
    mg.at_node['topographic__elevation'][:] = (mg.at_node["bedrock__elevation"]
                                               + mg.at_node["soil__depth"])
    
    vol_topo_pre_step = np.sum(mg.at_node['topographic__elevation'][mg.core_nodes]) * dx * dx

    # Run PriorityFloodFlowRouter, SPACE, DepthDependentDiffuser
    fr.run_one_step()
    sp.run_one_step(dt=dt)
    ddd.run_one_step(dt=dt)
    
    vol_topo_post_step = np.sum(mg.at_node['topographic__elevation'][mg.core_nodes]) * dx * dx
    
    T[i] = elapsed_time
    #Qs[i] = (vol_topo_pre_step - vol_topo_post_step) / dt
    Zm[i] = np.mean(mg.at_node['topographic__elevation'][mg.core_nodes])
    
    # Print progress and output grids
    if np.mod(elapsed_time, progress_interval) == 0:
        print(elapsed_time/1000," ka",("-- %s seconds --" % round((time.time() - start_time),1)))
        print("%.1f percent complete" % (elapsed_time*100 / total_t))
        print("Max elev: ", round(max(mg.at_node['topographic__elevation'][mg.core_nodes]),1))
    
        export_grids(mg,
                     params['grid_fields_to_save'],
                     path_export,
                     run_name_export, 
                     elapsed_time
                     )
    
    if np.mod(elapsed_time, output_interval) == 0 and elapsed_time > params['eq_interval']:
        end_run = check_for_eq(i, T_eq_window, Zm, Zm_slope_lims, end_run) #Qs, Qs_eq_lims, Qs_slope_lims,
        
        if end_run is True:

            break

export_grids_final(mg,
                   params['grid_fields_to_save'],
                   path_export,
                   run_name_export, 
                   elapsed_time
                   )

# Export record of parameters used in model run
with open(str(path_export.joinpath(run_name_export + '_Parameters.txt')), 'w') as file:
     file.write(yaml.dump(params, indent=4, sort_keys=True))

