# -*- coding: utf-8 -*-
"""
Created on Tue Oct 29 17:06:45 2024

@author: LaurentRoberge
"""

"""
Plot effects of different catchment size values on landscape, Qs timeseries, 
and Qs power spectra for JGR:ES paper: Roberge et al., 2025.
"""

# %%
# Import math & data libraries
import numpy as np
import pandas as pd

# Import plotting libraries
import matplotlib.pyplot as plt
import matplotlib.image as image
import matplotlib.gridspec as gridspec

# %% Set parameters

figure_title = 'Catchment size'
scenarios = ['Small (4 km$^2$)', 'Moderate (25 km$^2$)', 'Large (100 km$^2$)']

project_path = r'C:\Users\laure\Documents\Roberge_et_al_2025_JGR-ES'
figures_path = r'\Figures'

image_path = project_path + figures_path + '\\3D Views'
image_file_name = ['DA_4-Vm-Km-1e3.png',
                   'DA_25-Vm-Km-1e3.png',
                   'DA_100-Vm-Km-1e3.png',
                   ]

timeseries_path = project_path + '\\Timeseries'
timeseries_file_name = 'Qs_timeseries_Square_catchments_cleaned.csv'

mtm_spectra_path = project_path + '\\Spectra'
mtm_norm_file_name = 'Qs_50m_Square_nw7_MTMnorm.csv'
mtm_mean_file_name = 'Qs_50m_Square_nw7_MTMmean.csv'
mtm_period_file_name = 'Qs_50m_Square_nw7_period.csv'

col_names = ['EQ_DA4_Square_Vm_tLS1e3',
             'EQ_DA25_Square_Vm_tLS1e3',
             'EQ_DA100_Square_Vm_tLS1e3',
             ]

TRW = [2700,3000,3400]
TWB = [20000,16000,16000]

# %% LOAD DATA
# Load images into data structure
images = []
for file_name in image_file_name:
    images.append(image.imread(image_path + '\\' + file_name))
    
# Load time series into dataframe
timeseries = pd.read_csv(timeseries_path + '\\' + timeseries_file_name)
timeseries['Time (y)'] -= 2e6

# Load MTM power spectra into dataframe
mtm_period = pd.read_csv(mtm_spectra_path + '\\' + mtm_period_file_name, header=None)
mtm_period.loc[0, 'Period'] = 0
mtm_norm = pd.read_csv(mtm_spectra_path + '\\' + mtm_norm_file_name)
mtm_mean = pd.read_csv(mtm_spectra_path + '\\' + mtm_mean_file_name)
mtm_norm.insert(0,'Period',mtm_period[0])
mtm_mean.insert(0,'Period',mtm_period[0])

# %% Plot Figure 8 of paper

start_time = 0
stop_time = 2e6
start_time_zoom = 4e6
stop_time_zoom = 4.02e6

xlims = (start_time,stop_time)
xlims_zoom = (start_time_zoom,stop_time_zoom)
ylims = (-5,15)
ylims_zoom = (-4,7)

spectra_xlims = (6e1,1.5e5)
spectra_ylims = (1e-2,1e1)

text_row_1 = ['A)','B)','C)']
text_row_2 = ['D)','E)','F)']
text_row_3 = ['G)','H)','I)']
text_row_4 = ['J)','K)','L)']

ticks_major_location = np.linspace(start_time,stop_time,5)
ticks_minor_location = np.arange(0.25e6,stop_time,5e5)
ticks_major_labels = ticks_major_location/1e6
ticks_minor_labels = None

ticks_major_location_zoom = np.linspace(start_time_zoom,stop_time_zoom,5)
ticks_minor_location_zoom = np.arange(start_time_zoom+2.5e3,stop_time_zoom,5e5)
ticks_major_labels_zoom = np.linspace(0,20,5,dtype=int)
ticks_minor_labels_zoom = None

# Plot figure with subplots
fig = plt.figure(1)
fig.suptitle(figure_title + '\n', fontsize=16,)
fig.set_size_inches(w=12,h=12)
# set up subplot grid
gridspec.GridSpec(4,3)

for scenario in range(len(scenarios)):
    # Top row: 3D images of landscape grid
    plt.subplot2grid((4,3), (0,scenario))
    # plt.locator_params(axis='x', nbins=5)
    # plt.locator_params(axis='y', nbins=5)
    plt.text(160,-160,text_row_1[scenario],fontsize=14)
    plt.title(f'{scenarios[scenario]}')
    plt.imshow(images[scenario])
    plt.xticks([])
    plt.yticks([])
    plt.box(on=False)
    
    # Middle row: timeseries
    plt.subplot2grid((4,3), (1,scenario))
    # plt.locator_params(axis='x', nbins=5)
    # plt.locator_params(axis='y', nbins=5)
    # plt.title(f'{scenarios[scenario]}')
    plt.text(start_time+6e4,12.5,text_row_2[scenario],fontsize=14)
    plt.plot(timeseries['Time (y)'],timeseries[col_names[scenario]])
    plt.xlim(xlims)
    plt.ylim(ylims)
    plt.xlabel('Time (My)', fontsize=12, color='black')
    if scenario == 0:
        plt.ylabel('$Q_s$ anomaly (m$^3$/y)', fontsize=12, color='black')
    else:
        ax=plt.gca()
        ax.set_yticklabels([])
    plt.xticks(ticks=ticks_major_location,labels=ticks_major_labels,minor=False)
    plt.xticks(ticks_minor_location,labels=ticks_minor_labels,minor=True)
    plt.axhline(y=0, lw=0.5, linestyle=(0,(10,7)), color='k')
    
    # Middle row: timeseries (zoomed in)
    plt.subplot2grid((4,3), (2,scenario))
    # plt.locator_params(axis='x', nbins=5)
    # plt.locator_params(axis='y', nbins=5)
    # plt.title(f'{scenarios[scenario]}')
    plt.text(start_time_zoom+700,5.5,text_row_3[scenario],fontsize=14)
    plt.plot(timeseries['Time (y)'],timeseries[col_names[scenario]])
    plt.xlim(xlims_zoom)
    plt.ylim(ylims_zoom)
    plt.xlabel('Time (ky)', fontsize=12, color='black')
    if scenario == 0:
        plt.ylabel('$Q_s$ anomaly (m$^3$/y)', fontsize=12, color='black')
    else:
        ax=plt.gca()
        ax.set_yticklabels([])
    plt.xticks(ticks=ticks_major_location_zoom,labels=ticks_major_labels_zoom,minor=False)
    plt.xticks(ticks_minor_location_zoom,labels=ticks_minor_labels_zoom,minor=True)
    #plt.axhline(y=0, lw=0.5, linestyle=(0,(10,7)), color='k')
    
    # Bottom row: power spectra
    plt.subplot2grid((4,3), (3,scenario))
    # plt.locator_params(axis='x', nbins=5)
    # plt.locator_params(axis='y', nbins=5)
    # plt.title(f'{scenarios[scenario]}')
    plt.text(80,4,text_row_4[scenario],fontsize=14)
    plt.xlabel('Period (y)', fontsize=12, color='black')
    plt.loglog(mtm_mean['Period'],mtm_mean[col_names[scenario]])
    if TRW[scenario]:
        plt.loglog((TRW[scenario],TRW[scenario]),spectra_ylims,'r--')
    if TWB[scenario]:
        plt.loglog((TWB[scenario],TWB[scenario]),spectra_ylims,'b--')
    if scenario == 0:
        plt.ylabel('Power', fontsize=12, color='black')
    else:
        ax=plt.gca()
        ax.set_yticklabels([])
    plt.xlim(spectra_xlims)
    plt.ylim(spectra_ylims)

# Fit subplots and save fig
fig.tight_layout()
# fig_name = 'fig_8.png'
# fig.savefig(fig_name)
plt.show()