#ECH2O configuration file v2

# Please, check Appendix A of Documentation
# for units of parameters and variables
# (http://ech2o-iso.readthedocs.io/en/latest/Keywords.html)

#
#Folder section
#

#Maps_Folder = ./Input_Maps_50m
#Clim_Maps_Folder = ./Input_Climate
#Output_Folder = ./Results_50m_spin21Y_map5

## Water tracking (isotopes, age...)
#Tracking = 1
#TrackingConfig = ./Input_Configs/configTrck.ini

#
# Options section
#

MapTypes = csf
Species_State_Variable_Input_Method = tables # maps or tables

Map_Output_Format = pcraster # pcraster or netcdf
Ts_Output_Format = pcraster # pcraster or netcdf

#== Boolean switches
Reinfiltration = 1
Channel = 0
Channel_Infiltration = 0
# Are there lateral outputs (surface and saturated subsurface) from the domain ?
Lateral_Exit = 1
# Irrigation inputs ? Requires a separate forcings (+tracking if applicable)
Irrigation_Input = 0

# Limiting vertically routed at the the layers bottom to the 'reachable height' storage (Kv*dt)
Layer_thickness_for_conductivity = 0
# Experimental feature (may make the model unstable if activated): 
Root_Uptake_Redistribution = 0


# -- Toggle switches

# Hydraulic conductivity and porosity profiles --> 3 options each: 
# 0=constant, 1=exponential profile, 2=defined for each layer
Hydraulic_Conductivity_profile = 2
Porosity_profile = 2
# Pedotransfer parameter profiles (residual moisture, air-entry pressure
# and Brooks & Corey lambda), maps read for each layer if activated: 
Pedotransfer_Parameters_profile = 1
# Plant-level water uptake mode: 
# 0 = based root profile (kroot translated into layer fraction)
# 1 = using easing function based on water potential and depth
Uptake_Profile_opt = 1

# Vegetation dynamics (via allocation): 3 modes
# 0 -> desactivated, no dynamic allocation and constant LAI to initial value
# 1 -> fully activated
# 2 -> partially activated, except that LAI is prescribed via an input file
Vegetation_dynamics = 2
# Used only if Vegetation_dynamics = 2. Files names for each species is
# name below + "_"+ species number (starting at 0) + ".bin"
TimeSeries_LAI = LAI

# Aerodynamic resistance choices: 
# 0 = Penman Monteith option 
# 1 = Thom and Oliver 1977 
Aerodyn_resist_opt = 0 

# Soil resistance to vapor diffusion choices: 
# 0 = No resistance
# 1 = Passerat de Silans et al. 1989
# 2 = Sellers et al. 1992
# 3 = Sakaguchi and Zeng 2009 (CLM 3.5)
Soil_resistance_opt = 3 

#
# Time variables section
#
Simul_start = 0 # always 0
Simul_end = 646531200 # (seconds)
Simul_tstep = 86400 # seconds
Clim_input_tstep = 86400
Report_interval = 86400 #1 days
ReportMap_interval = 504921600 # gap between 2003-12-31 and 2023-12-24 (5844 days)
ReportMap_starttime = 86400 # from the start
#
# Climate input information
# Maps in this section to be contained in folder pointed by Clim_Maps_Folder
#
ClimateZones = ClimZones_P5.map
Isohyet_map = isohyet_P5.map
# ADDED LATER -- Isohyet_map = isohyet.map  # Precipitation multiplier map
Precipitation = Precip.bin  # Precip rate in meters/second
AirTemperature = Tavg.bin  # Average air temperature in degC
MaxAirTemp = Tmax.bin  # Maximum air temperature in degC
MinAirTemp = Tmin.bin  # Minimum air temperature in degC
RelativeHumidity = RH.bin  # air relative humidity in kPa/kPa
WindSpeed = Windspeed.bin  # Wind speed in meters/second
IncomingLongWave = Ldown.bin  # Downwelling longwave radiation in W/sq.meter
IncomingShortWave = Sdown.bin  # Solar radiation in W/sq.meter

#
# Spatial input information
# Maps below this line to be contained in folder pointed by Maps_Folder
#
# Site geometry : Digital Elevation Model and slope:
DEM = DEM.map
Slope = slope.map 
#
# Local drainage direction
local_drain_direc = ldd.map
#   
# Subsurface domain depth 
Soil_depth = soildepth.map 
Depth_soil_layer_1 = soildepth.L1.map 
Depth_soil_layer_2 = soildepth.L2.map 
#   
# Hydrologic Initial Conditions  
# Forest Initial states are included as maps or tables
#   
Streamflow = Init_Streamflow.map 
snow_water_equivalent = Init_SWE.map 
Soil_moisture_1 = Init_SWC.L1.map 
Soil_moisture_2 = Init_SWC.L2.map 
Soil_moisture_3 = Init_SWC.L3.map  
Soil_temperature = Init_SoilTemp.map 

# Channel parameters  
channel_width = chanwidth.map
channel_gw_transfer_param = chanparam.map # Back compatibility
channel_gw_transfer_param_Layer1 = chanparam.map
channel_gw_transfer_param_Layer2 = chanparam.map
channel_gw_transfer_param_Layer3 = chanparam.map
mannings_n = chanmanningn.map

# Subsurface parameters  
Horiz_Hydraulic_Conductivity = Khsat.map 
# if Ksat profile = 1
Horiz_Hydraulic_Conductivity_Profile_Coeff = kKhsat.map
# if Ksat profile = 2
Horiz_Hydraulic_Conductivity_Layer2 = Khsat.L2.map 
Horiz_Hydraulic_Conductivity_Layer3 = Khsat.L3.map

# if Ksat profile > 0
Vert_Horz_Anis_ratio = KvKh.map 
Vert_Horz_Anis_ratio_Layer2 = KvKh.L2.map 
Vert_Horz_Anis_ratio_Layer3 = KvKh.L3.map 
Soil_bedrock_leakance = leakance.map 

# porosity: full profile / top-of-profile / L1 value (resp.) 
# if poros profile = 0 / 1 / 2 (resp.)
Porosity = poros.map 
# if poros profile = 1
Porosity_Profile_Coeff = kporos.map
# if poros profile = 2
Porosity_Layer2 = poros.L2.map
Porosity_Layer3 = poros.L3.map
# if Porosity_profile > 0
Residual_soil_moisture_Layer2 = theta_r.L2.map
Residual_soil_moisture_Layer3 = theta_r.L3.map

# Water tension at field capacity (in meters)
Field_Capacity_Tension_Layer1 = psi_fc.L1.map
Field_Capacity_Tension_Layer2 = psi_fc.L2.map
Field_Capacity_Tension_Layer3 = psi_fc.L3.map

Air_entry_pressure = psi_ae.map 
Brooks_Corey_lambda = BClambda.map
Residual_soil_moisture = theta_r.map
# if pedotransfer param profile = 1
Air_entry_pressure_Layer2 = psi_ae.L2.map # if pedotransfer param profile = 1
Air_entry_pressure_Layer3 = psi_ae.L3.map
Brooks_Corey_lambda_Layer2 = BClambda.L2.map
Brooks_Corey_lambda_Layer3 = BClambda.L3.map

Veget_water_use_param1 = Wc.map 
Veget_water_use_param2 = Wp.map 

# Soil/surface energy balance
Terrain_Random_Roughness = randrough.map 
Albedo = albedo.soil.map
Surface_emissivity = emissivity.map
Dry_Soil_Heat_Capacity = soilheatcap.map
Dry_Soil_Therm_Cond = soilthermalK.map
Damping_depth = dampdepth.map
Temp_at_damp_depth = temp_damp.map
Snow_rain_temp_threshold = SnowRainTemp.map
Snow_Melt_Coeff = snowmeltCoeff.map

#   
#Forest Parameters and initial states 
#
ForestPatches = patches.map
Number_of_Species = 1
#Species_Parameters = SpeciesParams_new.tab  # Added by python script

#Tables below are only needed if Species_State_Variable_Input_Method = tables 
Species_Proportion_Table = SpecsProp.tab 
Species_StemDensity_Table = SpecsStemDens.tab 
Species_LAI_Table = SpecsLAI.tab 
Species_AGE_Table = SpecsAge.tab 
Species_BasalArea_Table = SpeciesBasalArea.tab 
Species_Height_table = SpeciesHeight.tab 
Species_RootMass_table = SpecsRootDensity.tab 
# 2. if Species_State_Variable_Input_Method = maps, then initial conditions
#    are read from maps with pre-established names: 
#    - p_i.map for species fraction cover
#    - ntr_i.map for stem density (tree m-2)
#    - lai_i.map for LAI (not needed of Vegetation_dynamics = 2)
#    - age_i.map for vegetation age (years)
#    - bas_i.map for stem basal area (m2)
#    - hgt_i.map for vegetation height (m)
#    - root_i.map for root density (g m-2)


# --------------------------------------------------------------
# == Report flags for output files
# (files will be written in be folder pointed by Output_Folder)

# Report catchment-scale vegetation state
Report_Vegetation_Summary = 1
#   

# Output file names, if netCDF is used
outMaps_Climate = Maps_Climate.nc
outMaps_Hydro_Storage = Maps_HydroStorage.nc
outMaps_Hydro_Fluxes = Maps_HydroFluxes.nc
outMaps_Cumulated_Hydro_Fluxes = Maps_HydroFluxesCumulated.nc
outMaps_Energy_Balance = Maps_Energy.nc
outMaps_Vegetation_State = Maps_Vegetation.nc
outTimeSeries_base = TimeSeries.nc
#   

# - Climate
Report_Long_Rad_Down = 0 
Report_Short_Rad_Down = 0 
Report_Precip = 0
Report_Rel_Humidity = 0
Report_Wind_Speed = 0
Report_AvgAir_Temperature = 0
Report_MinAir_Temperature = 0
Report_MaxAir_Temperature = 0
Report_Irrig = 0

# - Hydro Storage
Report_Canopy_Water_Stor_sum = 0 
Report_Canopy_Water_Stor = 0 
Report_SWE = 0 
Report_Saturation_Area = 0 
Report_Ponding = 0 
Report_Channel_Storage = 0 
Report_Pore_Relative_Humidity_L1 = 0
Report_Soil_Water_Content_Average = 0 
Report_Soil_Water_Content_Up = 0 
Report_Soil_Water_Content_L1 = 0
Report_Soil_Water_Content_L2 = 0 
Report_Soil_Water_Content_L3 = 0 
Report_WaterTableDepth = 0 
Report_WaterTableDepth_perched = 0 
Report_WaterTableDepth_L1 = 0 
Report_WaterTableDepth_L2 = 0 
Report_WaterTableDepth_L3 = 0 
Report_Soil_Sat_Deficit = 0 
Report_Ground_Water = 0 
# Maps of time-constant variables
Report_Field_Capacity_L1 = 0 
Report_Field_Capacity_L2 = 0 
Report_Field_Capacity_L3 = 0 

# - Hydro Fluxes
Report_Streamflow = 0 
Report_Total_ET = 0 
Report_Transpiration_sum = 0 
Report_Transpiration_Layer1 = 0 
Report_Transpiration_Layer2 = 0 
Report_Transpiration_Layer3 = 0 
Report_Einterception_sum = 0 
Report_Esoil_sum = 0 
Report_species_ET = 0 
Report_Transpiration = 0 
Report_Einterception = 0 
Report_Esoil = 0 
Report_Leakage_Out_of_System = 0

# Internal fluxes
Report_Throughfall = 0
Report_Snowmelt = 0 
Report_Infiltration = 0
Report_Infilt_Cap = 0 
Report_Percolation_to_Layer2 = 0
Report_Percolation_to_Layer3 = 0
Report_Groundwater_Recharge = 0
Report_GWseepage_Layer1_to_Channel = 0 
Report_GWseepage_Layer2_to_Channel = 0 
Report_GWseepage_Layer3_to_Channel = 0 
Report_Surface_to_Channel = 0 
Report_Return_Flow_Surface = 0
Report_Return_Flow_to_Layer1 = 0
Report_Return_Flow_to_Layer2 = 0
Report_Overland_Inflow = 0
Report_Stream_Inflow = 0
Report_Groundwater_Inflow = 0
Report_Groundwater_Inflow_Layer1 = 0 
Report_Groundwater_Inflow_Layer2 = 0 
Report_Groundwater_Inflow_Layer3 = 0 
Report_Overland_Outflow = 0
Report_Stream_Outflow = 0
Report_Groundwater_Outflow = 0
Report_Groundwater_Outflow_Layer1 = 0
Report_Groundwater_Outflow_Layer2 = 0
Report_Groundwater_Outflow_Layer3 = 0

# - Hydro Fluxes (cumulated)
Report_Transpiration_acc = 0
Report_Transpiration_Layer1_acc = 0
Report_Transpiration_Layer2_acc = 0
Report_Transpiration_Layer3_acc = 0
Report_E_Interception_acc = 0
Report_Soil_Evaporation_acc = 0
Report_Leakage_Out_of_System_acc = 0
Report_GW_to_Channel_acc = 0
Report_GWseepage_Layer1_to_Channel_acc = 0 
Report_GWseepage_Layer2_to_Channel_acc = 0 
Report_GWseepage_Layer3_to_Channel_acc = 0 
Report_Surface_to_Channel_acc = 0 
Report_Infiltration_acc = 0
Report_Percolation_to_Layer2_acc = 0
Report_Percolation_to_Layer3_acc = 0
Report_Groundwater_Recharge_acc = 0
Report_Return_Flow_Surface_acc = 0
Report_Return_Flow_to_Layer1_acc = 0
Report_Return_Flow_to_Layer2_acc = 0
Report_Overland_Inflow_acc = 0
Report_Stream_Inflow_acc = 0
Report_Groundwater_Inflow_acc = 0
Report_Groundwater_Inflow_Layer1_acc = 0
Report_Groundwater_Inflow_Layer2_acc = 0
Report_Groundwater_Inflow_Layer3_acc = 0
Report_Overland_Outflow_acc = 0
Report_Stream_Outflow_acc = 0
Report_Groundwater_Outflow_acc = 0
Report_Groundwater_Outflow_Layer1_acc = 0
Report_Groundwater_Outflow_Layer2_acc = 0
Report_Groundwater_Outflow_Layer3_acc = 0

# - Energy balance
Report_Surface_Net_Rad = 0 
Report_Vegetation_Net_Rad = 0 
Report_Total_Net_Rad = 0 
Report_Surface_Latent_Heat = 0 
Report_Vegetation_Latent_Heat = 0 
Report_Total_Latent_Heat = 0 
Report_Surface_Sens_Heat = 0 
Report_Vegetation_Sens_Heat = 0 
Report_Total_Sens_Heat = 0 
Report_Grnd_Heat = 0 
Report_Snow_Heat = 0 
Report_Soil_Temperature = 0 
Report_Skin_Temperature = 0 
Report_Canopy_Temp = 0 
Report_Canopy_NetR = 0 
Report_Canopy_LE_E = 0 
Report_Canopy_LE_T = 0 
Report_Canopy_Sens_Heat = 0 

# - Vegetation state
Report_GPP_sum = 0 
Report_NPP_sum = 0 
Report_Root_Allocation_sum = 0 
Report_Stem_Allocation_sum = 0 
Report_Leaf_Allocation_sum = 0 
Report_Root_Biomass_sum = 0 
Report_Stem_Biomass_sum = 0 
Report_Leaf_Biomass_sum = 0 
Report_Root_Decay_sum = 0 
Report_Leaf_Fall_sum = 0 
Report_Veget_frac = 0 
Report_Stem_Density = 0 
Report_RootFracL1_species = 0 
Report_RootFracL2_species = 0 
Report_Leaf_Area_Index = 0 
Report_Stand_Age = 0 
Report_Canopy_Conductance = 0 
Report_GPP = 0 
Report_NPP = 0 
Report_Basal_Area = 0 
Report_Tree_Height = 0 
Report_Root_Mass = 0 
# Maps of time-constant variables
Report_RootZone_in_L1 = 0 
Report_RootZone_in_L2 = 0 
Report_RootZone_in_L3 = 0 

# ---------------------------------------------------- 
# Report time section (at locations set by TS_mask map) 
#   

TS_mask = Tsmask.map 
#
Ts_OutletDischarge = 0
Ts_Long_Rad_Down = 0 
Ts_Short_Rad_Down = 0 
Ts_Precip = 0 
Ts_Rel_Humidity = 0 
Ts_Wind_Speed = 0 
Ts_AvgAir_Temperature = 0 
Ts_MinAir_Temperature = 0 
Ts_MaxAir_Temperature = 0
Ts_Irrig = 0

Ts_SWE = 0 
Ts_Snowmelt = 0
Ts_Infilt_Cap = 0 
Ts_Streamflow = 1 
Ts_Surface_Water = 0
Ts_Ponding = 0
Ts_Channel_Storage = 1
Ts_Pore_Relative_Humidity_L1 = 0
Ts_Soil_Water_Content_Average = 1 
Ts_Soil_Water_Content_Up = 0 
Ts_Soil_Water_Content_L1 = 1 
Ts_Soil_Water_Content_L2 = 1 
Ts_Soil_Water_Content_L3 = 1 
Ts_WaterTableDepth = 1
Ts_WaterTableDepth_perched = 1
Ts_WaterTableDepth_L1 = 1 
Ts_WaterTableDepth_L2 = 1 
Ts_WaterTableDepth_L3 = 1 
Ts_Field_Capacity_L1 = 0 
Ts_Field_Capacity_L2 = 0 
Ts_Field_Capacity_L3 = 0 
Ts_Soil_Sat_Deficit = 0 
Ts_Ground_Water = 0 
Ts_Surface_Net_Rad = 0
Ts_Vegetation_Net_Rad = 0
Ts_Total_Net_Rad = 0
Ts_Surface_Latent_Heat = 0 
Ts_Vegetation_Latent_Heat = 0
Ts_Total_Latent_Heat = 0
Ts_Surface_Sens_Heat = 0
Ts_Vegetation_Sens_Heat = 0
Ts_Total_Sens_Heat = 0
Ts_Grnd_Heat = 0 
Ts_Snow_Heat = 0 
Ts_Soil_Temperature = 0 
Ts_Skin_Temperature = 0

Ts_Total_ET = 1
Ts_Transpiration_sum = 1
Ts_Transpiration_Layer1 = 0
Ts_Transpiration_Layer2 = 0
Ts_Transpiration_Layer3 = 0 
Ts_Einterception_sum = 0
Ts_Esoil_sum = 1
Ts_Leakage_Out_of_System = 1

Ts_Canopy_Water_Stor_sum = 0 
Ts_GPP_sum = 0
Ts_NPP_sum = 0 
Ts_Root_Allocation_sum = 0 
Ts_Stem_Allocation_sum = 0 
Ts_Leaf_Allocation_sum = 0 
Ts_Root_Biomass_sum = 0
Ts_Stem_Biomass_sum = 0 
Ts_Leaf_Biomass_sum = 0 
Ts_Root_Decay_sum = 0
Ts_Leaf_Fall_sum = 0

Ts_Veget_frac = 0 
Ts_Stem_Density = 0 
Ts_RootFracL1_species = 0
Ts_RootFracL2_species = 0
Ts_Leaf_Area_Index = 1
Ts_Stand_Age = 0 
Ts_Canopy_Conductance = 0
Ts_GPP = 0 
Ts_NPP = 0 
Ts_Basal_Area = 0 
Ts_Tree_Height = 0 
Ts_Root_Mass = 0 
Ts_Canopy_Temp = 0
Ts_Canopy_NetR = 0
Ts_Canopy_LE_E = 0 
Ts_Canopy_LE_T = 0 
Ts_Canopy_Sens_Heat = 0
Ts_Canopy_Water_Stor = 0 
Ts_species_ET = 0
Ts_Transpiration = 0
Ts_Einterception = 0
Ts_Esoil = 0

# Internal fluxes
Ts_GW_to_Channel = 0
Ts_GWseepage_Layer1_to_Channel = 0 
Ts_GWseepage_Layer2_to_Channel = 0 
Ts_GWseepage_Layer3_to_Channel = 0 
Ts_Surface_to_Channel = 0 
Ts_Throughfall = 0
Ts_Infiltration = 0
Ts_Return_Flow_Surface = 0
Ts_Percolation_to_Layer2 = 0
Ts_Return_Flow_to_Layer1 = 0
Ts_Percolation_to_Layer3 = 0
Ts_Groundwater_Recharge = 0
Ts_Return_Flow_to_Layer2 = 0
Ts_Overland_Inflow = 0
Ts_Stream_Inflow = 0
Ts_Groundwater_Inflow = 0
Ts_Groundwater_Inflow_Layer1 = 0 
Ts_Groundwater_Inflow_Layer2 = 0 
Ts_Groundwater_Inflow_Layer3 = 0 
Ts_Overland_Outflow = 0
Ts_Stream_Outflow = 0
Ts_Groundwater_Outflow = 0
Ts_Groundwater_Outflow_Layer1 = 0
Ts_Groundwater_Outflow_Layer2 = 0
Ts_Groundwater_Outflow_Layer3 = 0
