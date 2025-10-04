# -*- coding: utf-8 -*-
"""
Created on Mon Oct 24 14:40:40 2022

@author: LaurentRoberge
"""

# calc_hillslopelength
# calc_hillchannelthreshold
# calc_hillchannelthreshold_main
# create_run_names
# export_grids
# export_grids_final
# check_for_eq
# export_run_time
# create_runtime_arrays
# create_empty_arrays
# calc_runtimes
# collect_stats_grid
# collect_stats_ls
# calculate_ls_soil_fraction
# combine_and_export_timeseries

from landlab.io.netcdf import write_netcdf

import numpy as np
import pandas as pd
from numpy.polynomial import Polynomial
from landlab.components import DrainageDensity
from SlopeBreaks import get_slopeBreaks

# Functions

def calc_hillslopelength(channel_threshold, base=10):
    
    return int(base * round(np.sqrt(channel_threshold)/base))


def calc_hillchannelthreshold(grid, base, custom_bins=False, custom_bin_edges=[]): 
    """ Calculate hillslope-channel threshold from Slope-Area relationship.
    Power trends are fitted to binned means of S-A and hinge point is returned.
    """
    
    # Define Slope and Drainage Area, log tranform them, remove NaN and inf.
    
    slope = grid.at_node["topographic__steepest_slope"][grid.core_nodes]
    da = grid.at_node["drainage_area"][grid.core_nodes]
    x1,y1 = np.log10(da),np.log10(slope)
    x1 = x1[~np.isnan(y1)]
    y1 = y1[~np.isnan(y1)]
    x1 = x1[~np.isinf(y1)]
    y1 = y1[~np.isinf(y1)]
    
    # Generate bins within which to average S-A values. Generate bin centres.
    if custom_bins: 
        bin_edges = custom_bin_edges
    else:
        bin_edges = np.linspace(1,10,40)
        
    binc_c = bin_edges[:-1]+ np.diff(bin_edges)/2
    
    # Generate array of slopes to fit within bins, calculate mean for each bin.
    slopes = np.zeros((len(bin_edges)-1))
    
    for bb in range(0, len(bin_edges)-1):
        slopes[bb] = np.mean(y1[np.logical_and(x1>bin_edges[bb], x1<bin_edges[bb+1])])
    
    binc_c = binc_c[~np.isnan(slopes)]
    slopes = slopes[~np.isnan(slopes)]
    
    # Use power law to fit data and find slope breaks (hinge points).
    my_pwlf, breaks, x_p, y_p = get_slopeBreaks(binc_c,slopes,plot=False, verbose=True)    
    
    if breaks is None:
        breaks = np.zeros(2)
    
    return  int(base * round((10**breaks[1])/base))

def calc_hillchannelthreshold_log(grid, custom_bins=False, custom_bin_edges=[]): 
    """ Calculate hillslope-channel threshold from Slope-Area relationship.
    Power trends are fitted to binned means of S-A and hinge point is returned.
    """
    
    # Define Slope and Drainage Area, log tranform them, remove NaN and inf.
    
    slope = grid.at_node["topographic__steepest_slope"][grid.core_nodes]
    da = grid.at_node["drainage_area"][grid.core_nodes]
    x1,y1 = np.log10(da),np.log10(slope)
    x1 = x1[~np.isnan(y1)]
    y1 = y1[~np.isnan(y1)]
    x1 = x1[~np.isinf(y1)]
    y1 = y1[~np.isinf(y1)]
    
    # Generate bins within which to average S-A values. Generate bin centres.
    if custom_bins: 
        bin_edges = custom_bin_edges
    else:
        bin_edges = np.linspace(1,10,40)
        
    binc_c = bin_edges[:-1]+ np.diff(bin_edges)/2
    
    # Generate array of slopes to fit within bins, calculate mean for each bin.
    slopes = np.zeros((len(bin_edges)-1))
    
    for bb in range(0, len(bin_edges)-1):
        slopes[bb] = np.mean(y1[np.logical_and(x1>bin_edges[bb], x1<bin_edges[bb+1])])
    
    binc_c = binc_c[~np.isnan(slopes)]
    slopes = slopes[~np.isnan(slopes)]
    
    # Use power law to fit data and find slope breaks (hinge points).
    my_pwlf, breaks, x_p, y_p = get_slopeBreaks(binc_c,slopes,plot=False, verbose=True)    
    
    if breaks is None:
        breaks = np.zeros(2)
    
    return  breaks[1]

def calc_hillchannelthreshold_main(grid, base, main_channel_ids): 
    """ Calculate hillslope-channel threshold from Slope-Area relationship.
    Power trends are fitted to binned means of S-A and hinge point is returned.
    """
    
    # Define Slope and Drainage Area, log tranform them, remove NaN and inf.
    
    slope = grid.at_node["topographic__steepest_slope"][main_channel_ids]
    da = grid.at_node["drainage_area"][main_channel_ids]

    x1,y1 = np.log10(da),np.log10(slope)
    x1 = x1[~np.isnan(y1)]
    y1 = y1[~np.isnan(y1)]
    x1 = x1[~np.isinf(y1)]
    y1 = y1[~np.isinf(y1)]
    
    # Generate bins within which to average S-A values. Generate bin centres.
    bin_edges = np.linspace(2.5,10,20)    
    binc_c = bin_edges[:-1]+ np.diff(bin_edges)/2
    
    # Generate array of slopes to fit within bins, calculate mean for each bin.
    slopes = np.zeros((len(bin_edges)-1))
    
    for bb in range(0, len(bin_edges)-1):
        slopes[bb] = np.mean(y1[np.logical_and(x1>bin_edges[bb], x1<bin_edges[bb+1])])
    
    binc_c = binc_c[~np.isnan(slopes)]
    slopes = slopes[~np.isnan(slopes)]
    
    # Use power law to fit data and find slope breaks (hinge points).
    my_pwlf, breaks, x_p, y_p = get_slopeBreaks(binc_c,slopes,plot=False, verbose=True)    
    
    if breaks is None:
        breaks = np.zeros(2)
    
    return  int(base * round((10**breaks[1])/base))



def create_run_names(params, params_sp, params_ls, runtype):

    DA_km2 = int(params['nrows']*params['ncols']*params['dx']*params['dx']/1000000)
    
    if params['nrows'] == params['ncols']:
        shape = '_Square'
    elif params['nrows'] == params['ncols']/4:
        shape = '_Long'
    else:
        raise Exception("'ncols' and 'nrows' are of incorrect dimensions.")
        
    if params_sp['v_s'] == 0.01:
        V_category = 'vs'
    elif params_sp['v_s'] == 0.1:
        V_category = 's'
    elif params_sp['v_s'] == 1:
        V_category = 'm'
    elif params_sp['v_s'] == 10:
        V_category = 'l'
    elif params_sp['v_s'] == 100:
        V_category = 'vl'
    else:
        raise Exception("Net effective settling velocity (v_s) must be 0.1, 1, or 10 m/yr.")
        
    if params_ls['landslides_return_time'] == 1e1:
        tLS_category = '1e1'
    elif params_ls['landslides_return_time'] == 1e3:
        tLS_category = '1e3'
    elif params_ls['landslides_return_time'] == 1e5:
        tLS_category = '1e5'
    elif params_ls['landslides_return_time'] == 1e7:
        tLS_category = '1e7'
    elif params_ls['landslides_return_time'] == 1e9:
        tLS_category = '1e9'
    else:
        raise Exception("Landslide return time (tLS) must be 10^ 1, 3, 5, 7, or 9 yrs.")
    
    if runtype == "BASE":
        
        run_name_import = "No run import necessary."

        run_name_export = ('BASE' + 
                            '_DA' + str(DA_km2) + shape +
                            '_V' + V_category
                            )
        
    elif runtype == "TRANSIENT":
        
        run_name_import = ('BASE' + 
                            '_DA' + str(DA_km2) + shape +
                            '_V' + V_category
                            )
        
        run_name_export = ('TRANSIENT' + 
                            '_DA' + str(DA_km2) + shape +
                            '_V' + V_category +
                            '_tLS' + tLS_category
                            )
        
    elif runtype == "EQ":
        
        run_name_import = ('TRANSIENT' + 
                            '_DA' + str(DA_km2) + shape +
                            '_V' + V_category +
                            '_tLS' + tLS_category
                            )
        
        run_name_export = ('EQ' + 
                            '_DA' + str(DA_km2) + shape +
                            '_V' + V_category +
                            '_tLS' + tLS_category
                            )
    
    return run_name_import, run_name_export

def export_grids(mg, grid_fields_to_save, path_export, run_name_export, elapsed_time):
    """Exports chosen fields from grid to NetCDF file."""
    
    fname = str(path_export.joinpath(run_name_export + '_Grid_' + str(int(elapsed_time/1000)) + 'ka.nc'))
    write_netcdf(fname, mg, names=grid_fields_to_save)

def export_grids_final(mg, grid_fields_to_save, path_export, run_name_export, elapsed_time):
    """Exports chosen fields from grid to NetCDF file."""
    
    fname = str(path_export.joinpath(run_name_export + '_Grid_FINAL.nc'))
    write_netcdf(fname, mg, names=grid_fields_to_save)
    
def check_for_eq(i, T_eq_window, Zm, Zm_slope_lims, end_run):
    """Checks slope of recent mean elevation to determine whether dynamic
    equilibrium has been achieved."""
    
    Zm_eq_window = Zm[i-len(T_eq_window) : i]
        
    p = Polynomial.fit(T_eq_window, Zm_eq_window, deg=1)
    Zm_intercept, Zm_slope = p.convert()
    
    if Zm_slope_lims[0] <= Zm_slope <= Zm_slope_lims[1]:
        end_run = True

    return end_run

def export_run_time(start_time, end_time, path_export):
    """Exports text file containing a string showing days, hrs, mins, secs."""
    
    run_time = end_time - start_time

    days = str(round(run_time/86400))
    days_remainder = np.remainder(run_time,86400)
    hrs = str(round(days_remainder/3600))
    hrs_remainder = np.remainder(run_time,3600)
    mins = str(round(hrs_remainder/60))
    mins_remainder = np.remainder(run_time,60)
    secs = str(round(mins_remainder,1))

    run_time_str = days + "d " + hrs + "h " + mins + "m " + secs + "s"
    
    print("Model run time: " + run_time_str)
    
    with open(str(path_export.joinpath('Runtime.txt')), 'w') as file:
         file.write(run_time_str)

# Functions for collecting and outputting data

def create_runtime_arrays(ndt):
    t_start = {"ew": np.zeros(ndt),
               "fr": np.zeros(ndt),
               "sp": np.zeros(ndt),
               "sf": np.zeros(ndt),
               "ddd": np.zeros(ndt),
               "ls": np.zeros(ndt),
               }
    t_end = {"ew": np.zeros(ndt),
               "fr": np.zeros(ndt),
               "sp": np.zeros(ndt),
               "sf": np.zeros(ndt),
               "ddd": np.zeros(ndt),
               "ls": np.zeros(ndt),
               }
    
    runtimes = {"ew": np.zeros(ndt),
               "fr": np.zeros(ndt),
               "sp": np.zeros(ndt),
               "sf": np.zeros(ndt),
               "ddd": np.zeros(ndt),
               "ls": np.zeros(ndt),
               }
    return t_start, t_end, runtimes

def create_empty_arrays(ndt, n_core_nodes):
    ts_data = {"Time": np.zeros(ndt),
               "Qs_outlet": np.zeros(ndt),
               "Qs_sp" : np.zeros(ndt),
               "Qs_ls" : np.zeros(ndt),
               "Qs_check" : np.zeros(ndt),
               "Wi" : np.zeros(ndt),           # Wiggle number (Kwang et al., 2021)
               "M_riv" : np.zeros(ndt),        # Mobility index (Campforts et al., 2022)
               "C_mem" : np.zeros(ndt),        # Degree of memory of initial landscape (Kwang et al., 2021)
               "Dd_simple" : np.zeros(ndt),    # Simple drainage density (length of "river" nodes to total catchment area)
               "Dd" : np.zeros(ndt),           # Drainage density calculated by distance-to-channel (Tucker et al., 2001)
               "Z_min" : np.zeros(ndt),
               "Z_max" : np.zeros(ndt),
               "Z_mean" : np.zeros(ndt),
               "Z_median" : np.zeros(ndt),
               "Soil_depth_in_channel_min" : np.zeros(ndt),
               "Soil_depth_in_channel_max" : np.zeros(ndt),
               "Soil_depth_in_channel_mean" : np.zeros(ndt),
               "Soil_depth_in_channel_median" : np.zeros(ndt),
               "Soil_depth_in_channel_total" : np.zeros(ndt),
               "Soil_depth_on_hillslopes_min" : np.zeros(ndt),
               "Soil_depth_on_hillslopes_max" : np.zeros(ndt),
               "Soil_depth_on_hillslopes_mean" : np.zeros(ndt),
               "Soil_depth_on_hillslopes_median" : np.zeros(ndt),
               "Soil_depth_on_hillslopes_total" : np.zeros(ndt),
               "Ksn_min" : np.zeros(ndt),
               "Ksn_max" : np.zeros(ndt),
               "Ksn_mean" : np.zeros(ndt),
               "Ksn_median" : np.zeros(ndt),
               "LS_number" : np.zeros(ndt),
               "LS_size_min" : np.zeros(ndt),
               "LS_size_max" : np.zeros(ndt),
               "LS_size_mean" : np.zeros(ndt),
               "LS_size_median" : np.zeros(ndt),
               "LS_size_total" : np.zeros(ndt),
               "LS_vol_min" : np.zeros(ndt),
               "LS_vol_max" : np.zeros(ndt),
               "LS_vol_mean" : np.zeros(ndt),
               "LS_vol_median" : np.zeros(ndt),
               "LS_vol_total" : np.zeros(ndt),
               "LS_vol_br_min" : np.zeros(ndt),
               "LS_vol_br_max" : np.zeros(ndt),
               "LS_vol_br_mean" : np.zeros(ndt),
               "LS_vol_br_median" : np.zeros(ndt),
               "LS_vol_br_total" : np.zeros(ndt),
               "LS_vol_sed_min" : np.zeros(ndt),
               "LS_vol_sed_max" : np.zeros(ndt),
               "LS_vol_sed_mean" : np.zeros(ndt),
               "LS_vol_sed_median" : np.zeros(ndt),
               "LS_vol_sed_total" : np.zeros(ndt),
               "LS_fraction_soil" : np.zeros(ndt),
               "Br_eroded_by_sp" : np.zeros(ndt),
               "Sed_eroded_by_sp" : np.zeros(ndt),
               "Br_Sed_check_sp" : np.zeros(ndt),
               "Br_eroded_by_ls" : np.zeros(ndt),
               "Sed_generated_by_ls" : np.zeros(ndt),
               "Br_Sed_check_ls" : np.zeros(ndt)
               }
    
    # Used for channel mobility/dynamism
    riv_data = {"Nodes_r" : np.zeros(n_core_nodes),    # nodes that are river (before fr)
                "Nodes_r2" : np.zeros(n_core_nodes),   # nodes that are rivers (after fr)
                "Nodes_r_nr" : np.zeros(n_core_nodes), # nodes that change from river to non-river
                "Nodes_nr_r" : np.zeros(n_core_nodes), # nodes that change from non-river to river
                }
    
    vol_data = {"vol_topo_pre_sp" : 0.,
                "vol_br_pre_sp" : 0.,
                "vol_sed_pre_sp" : 0.,
                "vol_topo_post_sp" : 0.,
                "vol_br_post_sp" : 0.,
                "vol_sed_post_sp" : 0.,
                "vol_topo_post_ls" : 0.,
                "vol_br_post_ls" : 0.,
                "vol_sed_post_ls" : 0.
            }
    
    return ts_data, riv_data, vol_data

def calc_runtimes(t_start, t_end, runtimes, i):
    
    runtimes['ew'][i] = t_end['ew'] - t_start['ew']
    runtimes['fr'][i] = t_end['fr'] - t_start['fr']
    runtimes['sp'][i] = t_end['sp'] - t_start['sp']
    runtimes['sf'][i] = t_end['sf'] - t_start['sf']
    runtimes['ddd'][i] = t_end['ddd'] - t_start['ddd']
    runtimes['ls'][i] = t_end['ls'] - t_start['ls']
    
    return runtimes

def collect_stats_grid(ts_data, vol_data, riv_data, elapsed_time, i, dt, dx, 
                       mg, uplift_per_step, node_next_to_outlet, n_core_nodes, 
                       channel_threshold_area, topo_pre_sp, topo_post_ls, 
                       topo_ref_val, topo_ref_val_sq
                       ):
    # Time (years)
    ts_data["Time"][i] = elapsed_time
    
    # Sediment flux at outlet
    ts_data["Qs_outlet"][i] = mg.at_node['sediment__flux'][node_next_to_outlet]
    
    # Mass conservation check of sediment flux at outlet
    ts_data["Qs_sp"][i] = (vol_data["vol_topo_pre_sp"] - vol_data["vol_topo_post_sp"]) / dt
    ts_data["Br_eroded_by_sp"][i] = (vol_data["vol_br_pre_sp"] - vol_data["vol_br_post_sp"]) / dt
    ts_data["Sed_eroded_by_sp"][i] = (vol_data["vol_sed_pre_sp"] - vol_data["vol_sed_post_sp"]) / dt
    ts_data["Br_Sed_check_sp"][i] = (ts_data["Br_eroded_by_sp"][i] + ts_data["Sed_eroded_by_sp"][i]) - ts_data["Qs_sp"][i]
    
    # River mobility index
    for n in range(n_core_nodes):
        riv_data["Nodes_r_nr"][n] = riv_data["Nodes_r"][n] == True and riv_data["Nodes_r2"][n] == False # nodes that changed from rivers to not rivers
        riv_data["Nodes_nr_r"][n] = riv_data["Nodes_r"][n] == False and riv_data["Nodes_r2"][n] == True # nodes that changed from not rivers to rivers
    
    nNodes_r = np.count_nonzero(riv_data["Nodes_r"])
    nNodes_r_nr = np.count_nonzero(riv_data["Nodes_r_nr"])
    nNodes_nr_r = np.count_nonzero(riv_data["Nodes_nr_r"])
    
    ts_data["M_riv"][i] = np.max([nNodes_r_nr,nNodes_nr_r]) / (dt * nNodes_r)
    
    # Wiggle number
    I = topo_pre_sp-topo_post_ls
    U_I = np.sum(abs(uplift_per_step - I))
    
    ts_data["Wi"][i] = U_I / (uplift_per_step*n_core_nodes)
        
    # Landscape memory correlation coefficient
    topo_mean = np.mean(topo_post_ls)
    topo_val = topo_post_ls - topo_mean
    topo_val_sq = np.sqrt(np.sum(np.square(topo_post_ls - topo_mean)))
    
    ts_data["C_mem"][i] = np.sum(topo_val*topo_ref_val) / (topo_val_sq*topo_ref_val_sq)
    
    # Drainage density (Ratio of channel length per area assumeing length of channel per river node is dx*sqrt(2))
    ts_data["Dd_simple"][i] = (nNodes_r*(dx*np.sqrt(2))) / (n_core_nodes*dx*dx)
    
    # Drainage density (using distance-to-channel method, Tucker et al., 2001)
    channel_mask = np.array(mg.at_node['drainage_area'] > channel_threshold_area, dtype=np.uint8)
    dd = DrainageDensity(mg, channel__mask=channel_mask)
    ts_data["Dd"][i] = dd.calculate_drainage_density()
    
    # Elevation statistics
    ts_data["Z_min"][i] = np.min(mg.at_node["topographic__elevation"][mg.core_nodes])
    ts_data["Z_max"][i] = np.max(mg.at_node["topographic__elevation"][mg.core_nodes])
    ts_data["Z_mean"][i] = np.mean(mg.at_node["topographic__elevation"][mg.core_nodes])
    ts_data["Z_median"][i] = np.median(mg.at_node["topographic__elevation"][mg.core_nodes])
    
    # Hillslope soil depth statistics
    hillslope_nodes = np.where(mg.at_node["drainage_area"][mg.core_nodes] < channel_threshold_area)
    channel_nodes = np.where(mg.at_node["drainage_area"][mg.core_nodes] >= channel_threshold_area)
    
    ts_data["Soil_depth_on_hillslopes_min"][i] = np.min(mg.at_node["soil__depth"][hillslope_nodes])
    ts_data["Soil_depth_on_hillslopes_max"][i] = np.max(mg.at_node["soil__depth"][hillslope_nodes])
    ts_data["Soil_depth_on_hillslopes_mean"][i] = np.mean(mg.at_node["soil__depth"][hillslope_nodes])
    ts_data["Soil_depth_on_hillslopes_median"][i] = np.median(mg.at_node["soil__depth"][hillslope_nodes])
    ts_data["Soil_depth_on_hillslopes_total"][i] = np.sum(mg.at_node["soil__depth"][hillslope_nodes])
    
    # In-channel soil depth statistics
    ts_data["Soil_depth_in_channel_min"][i] = np.min(mg.at_node["soil__depth"][channel_nodes])
    ts_data["Soil_depth_in_channel_max"][i] = np.max(mg.at_node["soil__depth"][channel_nodes])
    ts_data["Soil_depth_in_channel_mean"][i] = np.mean(mg.at_node["soil__depth"][channel_nodes])
    ts_data["Soil_depth_in_channel_median"][i] = np.median(mg.at_node["soil__depth"][channel_nodes])
    ts_data["Soil_depth_in_channel_total"][i] = np.sum(mg.at_node["soil__depth"][channel_nodes])

    # Ksn statistics
    ts_data["Ksn_min"][i] = np.min(mg.at_node["channel__steepness_index"][channel_nodes])
    ts_data["Ksn_max"][i] = np.max(mg.at_node["channel__steepness_index"][channel_nodes])
    ts_data["Ksn_mean"][i] = np.mean(mg.at_node["channel__steepness_index"][channel_nodes])
    ts_data["Ksn_median"][i] = np.median(mg.at_node["channel__steepness_index"][channel_nodes])

def collect_stats_ls(ts_data, vol_data, i, dt, dx, ls):
    
    ts_data["LS_number"][i] = np.size(ls.landslides_size)
    
    # Mass conservation check of Qs from landslides and for the whole timestep
    ts_data["Qs_ls"][i] = (vol_data["vol_topo_post_sp"] - vol_data["vol_topo_post_ls"]) / dt
    ts_data["Qs_check"][i] = ts_data["Qs_sp"][i] + ts_data["Qs_ls"][i]
    ts_data["Br_eroded_by_ls"][i] = (vol_data["vol_br_post_sp"] - vol_data["vol_br_post_ls"]) / dt
    ts_data["Sed_generated_by_ls"][i] = (vol_data["vol_sed_post_ls"] - vol_data["vol_sed_post_sp"]) / dt
    ts_data["Br_Sed_check_ls"][i] = ts_data["Br_eroded_by_ls"][i] - (ts_data["Sed_generated_by_ls"][i] + ts_data["Qs_ls"][i])
    
    if np.size(ls.landslides_size) > 0:
        ts_data["LS_size_min"][i] = np.min(ls.landslides_size) * dx**2
        ts_data["LS_size_max"][i] = np.max(ls.landslides_size) * dx**2
        ts_data["LS_size_mean"][i] = np.mean(ls.landslides_size) * dx**2
        ts_data["LS_size_median"][i] = np.median(ls.landslides_size) * dx**2
        ts_data["LS_size_total"][i] = np.sum(ls.landslides_size) * dx**2
    
        ts_data["LS_vol_min"][i] = np.min(ls.landslides_volume)
        ts_data["LS_vol_max"][i] = np.max(ls.landslides_volume)
        ts_data["LS_vol_mean"][i] = np.mean(ls.landslides_volume)
        ts_data["LS_vol_median"][i] = np.median(ls.landslides_volume)
        ts_data["LS_vol_total"][i] = np.sum(ls.landslides_volume)
    
        ts_data["LS_vol_br_min"][i] = np.min(ls.landslides_volume_bed)
        ts_data["LS_vol_br_max"][i] = np.max(ls.landslides_volume_bed)
        ts_data["LS_vol_br_mean"][i] = np.mean(ls.landslides_volume_bed)
        ts_data["LS_vol_br_median"][i] = np.median(ls.landslides_volume_bed)
        ts_data["LS_vol_br_total"][i] = np.sum(ls.landslides_volume_bed)
    
        ts_data["LS_vol_sed_min"][i] = np.min(ls.landslides_volume_sed)
        ts_data["LS_vol_sed_max"][i] = np.max(ls.landslides_volume_sed)
        ts_data["LS_vol_sed_mean"][i] = np.mean(ls.landslides_volume_sed)
        ts_data["LS_vol_sed_median"][i] = np.median(ls.landslides_volume_sed)
        ts_data['LS_vol_sed_total'][i] = np.sum(ls.landslides_volume_sed)
        

def calculate_ls_soil_fraction(ts_data):
    ts_data["LS_fraction_soil"] = ts_data["LS_vol_sed_total"] / ts_data["LS_vol_total"]
    ts_data["LS_fraction_soil"] = np.nan_to_num(ts_data["LS_fraction_soil"])


def combine_and_export_timeseries(ts_data, path_export, run_name_export):
    col_data = ts_data.values()
    col_names = ["Time (y)",
                 "Qs at outlet (m^3/y)",
                 "Qs from fluvial incision (m^3/y)",
                 "Qs from landslides (m^3/y)",
                 "Qs from mass balance (m^3/y)",
                 "EXPERIMENTAL Wiggle number (-)",
                 "EXPERIMENTAL River mobility index (fraction of channel changing per year)",
                 "EXPERIMENTAL Landscape memory correlation coefficient (-)", 
                 "EXPERIMENTAL Drainage density simple (m^-1)",
                 "EXPERIMENTAL Drainage density (m^-1)",
                 "Elevation min (m)",
                 "Elevation max (m)",
                 "Elevation mean (m)",
                 "Elevation median (m)",
                 "Soil depth in channel min (m)",
                 "Soil depth in channel max (m)",
                 "Soil depth in channel mean (m)",
                 "Soil depth in channel median (m)",
                 "Soil depth in channel total (m)",
                 "Soil depth on hillslopes min (m)",
                 "Soil depth on hillslopes max (m)",
                 "Soil depth on hillslopes mean (m)",
                 "Soil depth on hillslopes median (m)",
                 "Soil depth on hillslopes total (m)",
                 "Ksn min (m^0.9)",
                 "Ksn max (m^0.9)",
                 "Ksn mean (m^0.9)",
                 "Ksn median (m^0.9)",
                 "Number of landslides",
                 "Landslide size min (m^2)",
                 "Landslide size max (m^2)",
                 "Landslide size mean (m^2)",
                 "Landslide size median (m^2)",
                 "Landslide size total (m^2)",
                 "Landslide volume min (m^3)",
                 "Landslide volume max (m^3)",
                 "Landslide volume mean (m^3)",
                 "Landslide volume median (m^3)",
                 "Landslide volume total (m^3)",
                 "Landslide volume of bedrock min (m^3)",
                 "Landslide volume of bedrock max (m^3)",
                 "Landslide volume of bedrock mean (m^3)",
                 "Landslide volume of bedrock median (m^3)",
                 "Landslide volume of bedrock total (m^3)",
                 "Landslide volume of soil min (m^3)",
                 "Landslide volume of soil max (m^3)",
                 "Landslide volume of soil mean (m^3)",
                 "Landslide volume of soil median (m^3)",
                 "Landslide volume of soil total (m^3)",
                 "Landslide fraction of soil by volume (-)",
                 "Bedrock eroded by SP (m^3/y)",
                 "Sediment eroded by SP (m^3/y)",
                 "Mass balance check: SP bedrock and soil (m^3/y)",
                 "Bedrock eroded by LS (m^3/y)",
                 "Bedrock converted to sediment by LS (m^3/y)",
                 "Mass balance check: LS bedrock to soil (m^3/y)"
                 ]
    
    ts_data = pd.DataFrame.from_dict(dict(zip(col_names, col_data)))
    fname = str(path_export.joinpath(run_name_export + '_DataTimeseries.csv'))
    ts_data.to_csv(fname, na_rep='NaN', index=False)
