--- INSTRUCTIONS ---
The LEM must be run as three successive scripts, each requiring input from the previous.
- LEM_1 builds a base landscape using the deterministic hillslope and fluvial processes,
  stopping when it has achieved static steady state.
- LEM_2 takes the initial landscape and adds landslides, causing a transient period as
  the landscape adjusts. It stops when the landscape achieves a dynamic steady state.
- LEM_3 is the final model scenario. It runs for 22 million years of dynamic equilibrium
  during which all metrics are calculated and periodic landscape outputs are saved.

IMPORTANT: Make sure that LEM_3 has the same input parameters as LEM_2. Otherwise, there
will be a transient period at the beginning of the dynamic equilibrium run, which will
cause problems for the timeseries analysis.

STEP 1: Place all files in the same folder:
   - Input files (input_INITIAL, input_TRANSIENT, input_EQ)
   - LEM files (LEM_1, LEM_2, LEM_3)
   - Function libraries (library.py, SlopeBreaks.py)
STEP 2: Run "LEM_1_build_initial_landscape.py" with "input INITIAL" as input.
STEP 3: Run "LEM_2_transient_with_landslides.py" with "input INITIAL" as input.
STEP 4: Update "input_EQ" parameter: "channel_threshold_area" to the correct value
        based on a reasonable threshold area from the end of the LEM_2 model run.
STEP 5: Run "LEM_3_equilibrium_with_landslides.py" with "input INITIAL" as input.

NOTE: Several outputs metrics are tagged as "EXPERIMENTAL" as they have not been 
tested and I cannot guarantee they are calculated correctly. They are:
- Wiggle number (-)
- River mobility index (fraction of channel changing per year)
- Landscape memory correlation coefficient
- Drainage density simple
- Drainage density

# --------------------------------
# Description of input parameters
# --------------------------------
# Grid, time, and other parameters

nrows                   : Number of rows
ncols                   : Number of columns
uplift_rate             : Uplift rate (m/yr)
channel_threshold_area  : Minimum drainage area to define a channel (m^2)

dx                      : Node spacing (m)
dt                      : Timestep (yrs)
total_t                 : Total run time (yrs)
grid_seed               : Random seed used to generate micro-topography
progress_interval       : Interval at which to print progress (yrs)
output_interval         : Interval at which to output grids (yrs)
soil_depth_initial      : Initial soil depth (m)

# ------------------------------------
params_sp:   # Parameters for SpaceLargeScaleEroder

  v_s         : Net effective settling velocity (m/yr) 

  m_sp        : Slope exponent (-)
  n_sp        : Area exponent (-)
  K_sed       : Sediment erodibility (m^-1)
  K_br        : Bedrock erodibility (m^-1)
  F_f         : Fraction of fine fluvial sediment (-)
  phi         : Porosity of the bed sediment (-)
  H_star      : Bedrock roughness length scale (m)
  v_s_lake    : Net effective settling velocity in lakes (m/yr)
  sp_crit_sed : Entrainment threshold for sediment (E/(yr m^2))
  sp_crit_br  : Erosion threshold for bedrock (E/(yr m^2))

# ------------------------------------
params_ls:   # Parameters for BedrockLandslider

  angle_int_frict              : Tangent of the angle of internal friction (m/m)
  landslides_return_time       : Mean recurrence interval between landslide trigger events (y)

  min_deposition_slope         : Tangent of the minimal deposition angle (m/m)
  threshold_slope              : Threshold slope for landslide deposition (m/m)

  cohesion_eff                 : Cohesion (kPa)
  fraction_fines_LS            : Fraction of fine hillslope sediment (-)
  phi                          : Porosity of sediment (-)
  max_pixelsize_landslide      : Maximum size for landslides (# pixels)
  seed                         : Seed to set stochastic model
  verbose_landslides           : Print number of simulated landslides per timestep (boolean)
  landslides_on_boundary_nodes : Allow/disallow landslides to act on boundary nodes (boolean)

# ------------------------------------
params_fr:   # Parameters for PriorityFloodFlowRouter

  flow_metric             : Flow routing method: 'D8','Rho8','Quinn','Freeman','Holmgren','Dinf'
  update_flow_depressions : Activate/deactivate depression handler (boolean)
  depression_handler      : Method for dealing with depressions: 'fill or 'breach'
  accumulate_flow_hill    : Activate/deactivate flow accumulations for hillslope flow direction component (boolean)
  separate_hill_flow      : Activate/deactivate both single and multiple flow direction/accumulation (boolean)
  update_hill_depressions : Activate/deactivate hillslopw depression handler (boolean)
  hill_flow_metric        : Flow routing method for hillslopes: 'D4','D8','Rho4','Rho8','Quinn','Freeman','Holmgren','Dinf'

# ------------------------------------
params_ew:   # Parameters for ExponentialWeatherer

  soil_production__maximum_rate  : Maximum soil production rate (m)
  soil_production__decay_depth   : Soil production decay depth (m)

# ------------------------------------
params_ddd:   # Parameters for DepthDependentDiffuser

  linear_diffusivity          : Diffusion coefficient (m^2/yr)
  soil_transport_decay_depth  : Characteristic soil depth for soil creep (m)