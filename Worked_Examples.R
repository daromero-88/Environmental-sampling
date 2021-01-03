# WORKED EXAMPLES
### This script contains a worked example on how to use all the functions created to design
### a survey of potential sites inside a region of interest. The selection is based on points
### or suitability values from a Species Distribution Model (SDM) and can be focused in either
### covering all the range of suitability values, if the goal is to improve the estimation
### of the SDM, or covering the regions of high suitability, if the goal is to select sites
### where the likelihood of detecting the species is high.

# Read R code ----------------
# These three scripts contain the functions needed for the analyses
source("E-space-functions.R")
source("Hutchinson-functions.R")
source("Post-track-functions.R")

# Load packages -------------
library(gatepoints)
library(raster)
library(rgdal)
library(maptools)
library(plyr)
library(sp)

# Set global parameters
# worldwide shapefile 
data("wrld_simpl", package = "maptools")

# Read data for all examples----------------------

#' environemntal_variables: The set of environmental variables selected in the modeling process
all_envs = stack(list.files('./environmental_variables', full.names = T))
# number of layers
length(all_envs@layers)

#Read data for examples 1 and 2: MELIOIDOSIS SAMPLING--------------------

#' categorized_models: models categorized using four thresholds as described in the main text
#' Minimum training presence, E= 5%, E= 10%, and E= 20% 
cat_models = stack(list.files('./categorized_models', full.names = T))
# number of layers
length(cat_models@layers)

#' uncertainty_models: Interquartile range of bootstraped models using final model selected parameters
uncert_models = stack(list.files ('./uncertainty_models', full.names = T))
# number of layers
length(uncert_models@layers)

#Read data for example 3: POINTS TO ENVIRONMENTS AND SAMPLING----------------

ldb = read.csv ('./env_samp_points/spp_ex1.csv')
amr = readOGR ('./env_samp_points/shape_amrs/Americas.shp')
amr = buffer(amr, width = 0.2, dissolve = T)


# WORKED EXAMPLE 1, FULLY COMMENTED: Ceara, Brazil ---------------------

# Read shapefiles and select region of interest

# World shapefiles from Natural Earth: 
provs = readOGR('./shapefiles2/ne_10m_admin_1_states_provinces.shp')

# Select specific region: CEARA
#   Provinces of Brazil
BR = provs@data[which(provs@data$admin =='Brazil'),] 
#   Look for the specifi province
BR$name
#  Select specific province by subsetting the data frame
cear = subset (provs, name == 'Ceará') #if problems with UNICODE, change name accordingly
# Create a buffer around the region of interest. Warnings are related with the planar projection. 
# this is to guarantee that all pixels in the raster layers are captured
cear_buf = buffer(cear, width = 0.2, dissolve = T)
# you can used different buffer sizes
cear_buf2 = buffer(cear, width = 0.5, dissolve = T) 
# Plot regions to verify our selection
#x11() #if needed for plotting 
plot(cear)
plot(cear_buf,add=T)
plot(cear_buf2,add=T)

# Crop and mask rasters -----------------

# Environmental variables
cear_envs = crop(all_envs, cear_buf)
cear_envs = mask(cear_envs, cear_buf)
# plot one of the layers to verify that the process was done correctly
plot (cear_envs[[1]])

# Categorized rasters
# Note: cear_mods contains 3 raster files corresponding to categorized models 
# using the no extrapolation (NE), extrapolation (E), and extrapolation and 
# clamping (EC) features from Maxent to transfer models to a given area. 
# We will work with the third raster which corresponds to NE. 
cear_mods = crop(cat_models, cear_buf)
cear_mods= mask(cear_mods, cear_buf)
# plot one of the layers to verify that the process was done correctly
plot (cear_mods[[3]]) 

# Uncertainty maps
# The uncertainty plots are available for all the transfered rasters including 
# NE, E, and EC. We will use the NE layer throughout this example. 
cear_unc = crop(uncert_models, cear_buf)
cear_unc = mask(cear_unc, cear_buf)
# plot one of the layers to verify that the process was done correctly
plot(cear_unc[[3]])

# Using functions -------------------------------------

# Select color ramp to be used in the visualizations
col <- c("blueviolet", "springgreen3")

# FUNCTION 1: e_space

# Displaying the environmental variables in E-space
# and the region covered by the layers in G-space
f1_cear = e_space(stck = cear_envs,pflag = T)

# FUNCTION 2: e_space_cat

# Displaying the environmental variables in E-space using different colors for different
# suitability categories and the sites that cover the region of interest also labeled
# according to its suitability category. This function recovers the database needed for 
# the Hutchinson functions. 
f2_cear = e_space_cat(stck = cear_envs, ctgr = cear_mods[[3]], pflag = T, col.use = col)

# FUNCTION 3: e_space_cat_back
# The output plots are similar to the ones obtained with the function e_space_cat()
# except that in this case there is a new category of points called Background.
# Background points are inside the region of interest but they have no suitability
# values assigned by the model.
f3_cear = e_space_cat_back(stck = cear_envs, ctgr = cear_mods[[3]],
                           bck = cear_envs, pflag = T, col.use = col)

# FUNCTION 4: hutchinson

# Option 1: from E-space to G-space
# Select the columns that contain the enviromental variables to be plotted
# and use them as the 'calls' argument of this function
names (f2_cear) 
# temperature & humidity
cear_tmp_hum = hutchinson_cat(EtoG=T, data=f2_cear, calls=c(6,4), plyg=cear, ntr=3, col.use=col) #if NA, na.omit will be neccesary, try to avoid the presence of NA in the database
# temperature & soil
cear_tmp_soil = hutchinson_cat(EtoG=T, data=f2_cear, calls=c(6,5), plyg=cear, ntr=3, col.use=col)
# humidity & soil
cear_hum_soil = hutchinson_cat(EtoG=T, data=f2_cear, calls=c(4,5), plyg=cear, ntr=3, col.use=col)

# Option 2: from G-space to E-space
# temperature & humidity
cear2_tmp_hum = hutchinson_cat(EtoG=F, data=f2_cear, calls=c(6,4), plyg=cear, ntr=3, col.use=col) #if NA, na.omit will be neccesary, try to avoid the presence of NA in the database
# temperature & soil
cear2_tmp_soil = hutchinson_cat(EtoG=F, data=f2_cear, calls=c(6,5), plyg=cear, ntr=3, col.use=col)
# humidity & soil
cear2_hum_soil = hutchinson_cat(EtoG=F, data=f2_cear, calls=c(4,5), plyg=cear, ntr=3, col.use=col)

# The sampling exercise has the goal of maximizing the selection of different suitability
# categories with the selection of different transects in either E-space or G-space.
# Therefore, once the transects are selected, a final step is needed to check the levels
# of uncertainty in the selected sites:

# FUNCTION 5: post_track_cat()

# Combine the dataframes if your are willing to include more than two environmental variables
cear_sampling = rbind(cear_tmp_hum, cear_tmp_soil, cear_hum_soil)
dim(cear_sampling)

# Select uncertainty layer and apply function 
uncer_check = post_track_cat(cear_sampling, cear_unc[[3]], cear, col.use=col)

#' Because different environmental tracks were selected using different environmental 
#' variables, some information might be repeated, however, the post_track_cat function
#' eliminates duplicates
dim(uncer_check)

# Save resulting tables for further analyses and visualizations
write.csv(f2_cear, './ceara_df1.csv', row.names = F)
write.csv(uncer_check, './ceara_res.csv', row.names = F)

# WORKED EXAMPLE 2: Texas, US ---------------------

# Read shapefiles and select region of interest

# World shapefiles from Natural Earth: 
provs = readOGR('./shapefiles2/ne_10m_admin_1_states_provinces.shp')

# Select specific region: TEXAS
#   STATES OF THE US
US = provs@data[which(provs@data$admin =='United States of America'),] 
#   Look for the specifi province
US$name
#  Select specific province by subsetting the data frame
tex1 = subset (provs, name == 'Texas')

#buffer
tex1_buf = buffer(tex1, width = 0.2, dissolve = T)

#plotting
plot(tex1)
plot(tex1_buf,add=T)


# Preparing data

# Environmental variables
tex1_envs = crop(all_envs, tex1_buf)
tex1_envs = mask(tex1_envs, tex1_buf)
# plot one of the layers to verify that the process was done correctly
plot (tex1_envs[[1]])

# Categorized rasters
tex1_mods = crop(cat_models, tex1_buf)
tex1_mods= mask(tex1_mods, tex1_buf)
# plot one of the layers to verify that the process was done correctly
plot (tex1_mods[[3]]) 

# Uncertainty maps
tex1_unc = crop(uncert_models, tex1_buf)
tex1_unc = mask(tex1_unc, tex1_buf)
# plot one of the layers to verify that the process was done correctly
plot(tex1_unc[[3]])

# Using functions 

# Select color ramp to be used in the visualizations
col <- c("blueviolet", "springgreen3")

# FUNCTION 1: e_space

f1_tex1 = e_space(stck = tex1_envs,pflag = T)

# FUNCTION 2: e_space_cat

f2_tex1 = e_space_cat(stck = tex1_envs, ctgr = tex1_mods[[3]], pflag = T, col.use = col)

# FUNCTION 3: e_space_cat_back

f3_tex1 = e_space_cat_back(stck = tex1_envs, ctgr = tex1_mods[[3]],
                           bck = tex1_envs, pflag = T, col.use = col)

# FUNCTION 4: hutchinson

# Option 1: from E-space to G-space
# Select the columns that contain the enviromental variables to be plotted
# and use them as the 'calls' argument of this function
names (f2_tex1) 
# temperature & humidity
tex1_tmp_hum = hutchinson_cat(EtoG=T, data=f2_tex1, calls=c(6,4), plyg=tex1, ntr=3, col.use=col) #if NA, na.omit will be neccesary, try to avoid the presence of NA in the database
# temperature & soil
tex1_tmp_soil = hutchinson_cat(EtoG=T, data=f2_tex1, calls=c(6,5), plyg=tex1, ntr=3, col.use=col)
# humidity & soil
tex1_hum_soil = hutchinson_cat(EtoG=T, data=f2_tex1, calls=c(4,5), plyg=tex1, ntr=3, col.use=col)

# Option 1: from G-space to E-space
# temperature & humidity
tex12_tmp_hum = hutchinson_cat(EtoG=F, data=f2_tex1, calls=c(6,4), plyg=tex1, ntr=3, col.use=col) #if NA, na.omit will be neccesary, try to avoid the presence of NA in the database
# temperature & soil
tex12_tmp_soil = hutchinson_cat(EtoG=F, data=f2_tex1, calls=c(6,5), plyg=tex1, ntr=3, col.use=col)
# humidity & soil
tex12_hum_soil = hutchinson_cat(EtoG=F, data=f2_tex1, calls=c(4,5), plyg=tex1, ntr=3, col.use=col)

# FUNCTION 5: post_track_cat()

# Combine the dataframes if your are willing to include more than two environmental variables
tex1_sampling = rbind(tex1_tmp_hum, tex1_tmp_soil, tex1_hum_soil)
dim(tex1_sampling)

# Select uncertainty layer and apply function
uncer_check = post_track_cat(tex1_sampling, tex1_unc[[3]], tex1, col.use=col)

#post_track_cat function eliminates duplicates
dim(uncer_check)

# Save resulting tables for further analyses and visualizations
write.csv(f2_tex1, './tex1_df1.csv', row.names = F)
write.csv(uncer_check, './tex1_res.csv', row.names = F)


#WORKED EXAMPLE 3: points to environmental sampling ------------------

#creatinng main data for Hutchinson-based sampling: 

#environmental variables: corpping word variables to the Americas
wrd_merra2 = crop (all_envs, amr)
wrd_merra2 = mask (wrd_merra2, amr)

#creating main dataframe: visualizing points in E and G space

rr = e_space(ldb[,2:3], wrd_merra2, pflag = T)

#examining only a stack of raster in E and G space 

#new extent for examples: 
ee = extent(-50, 50, -50, 50)
ee2 = extent (-60, -30, -70,-10)

#raster stack cropped with the created extent 
ex_stck = crop (wrd_merra2, ee) 
ex_stck2 = crop (wrd_merra2, ee2)

#check the environmental space of the raster stack: 
rr1 = e_space(stck = ex_stck, pflag = T) #in the absence of coordinates, please define the raster as shown here: 'stck= ex_stck' 

#examining points in G space and in E overlapped with a background:
#here we are plotting the points with values of the americas (wrd_merra2) in a selected background (ex_stck, the extent cropped above)
rr2 = e_space_back(ldb[,2:3], wrd_merra2, ex_stck, pflag = T)

#In the absence of points, we can compare one stack, against other 
#here we are plotting one stack (ex_stack2) with a background (ex_stack) 
#without points, directly using the raster stacks. Notice that you have to define each argument as in the example. 
rr3 = e_space_back(stck= ex_stck2, bck = ex_stck, pflag = T)

#Hutchinson sampling: From E to G
#sampling scheme: 
qq1 = hutchinson(EtoG = T, na.omit(rr), c(3,4), amr, 2) #if NA, na.omit will be neccesary, try to avoid the presence of NA in the database
qq2 = hutchinson(EtoG = T, na.omit(rr), c(3,5), amr, 2) 
qq3 = hutchinson(EtoG = T, na.omit(rr), c(4,5), amr, 2)

#Hutchinson sampling: from G to E
qq4 = hutchinson(EtoG = F, na.omit(rr), c(3,4), amr, 2) #if NA, na.omit will be neccesary, try to avoid the presence of NA in the database
qq5 = hutchinson(EtoG = F, na.omit(rr), c(3,5), amr, 2) 
qq6 = hutchinson(EtoG = F, na.omit(rr), c(4,5), amr, 2)

#defining database with transects from different environmental dimensions (E to G): 
q_df = rbind (qq1, qq2, qq3)
dim(q_df)


#trimming database for final coordinates and uncertainty/raster evaluation: 
ss = post_track(q_df, wrd_merra2[[1]], amr)
dim(ss)

#Manually trimming the database: 
#eliminating duplicates
q_df$dup = paste(q_df[,1], q_df[,2], sep= '_')
q_df2 = q_df[!duplicated(q_df$dup),]
dim(q_df2)
q_df2$dup = NULL


#Manually plotting individual results:
plot (amr)
points(q_df2[,1:2], pch = q_df2[,dim(q_df2)[2]], col= c(q_df2[,dim(q_df2)[2]]))
nm_lg = paste('Transect', unique(q_df2[,dim(q_df2)[2]]))
legend('bottomleft', legend = nm_lg, pch = unique(q_df2[,dim(q_df2)[2]]), 
       col = c(unique(q_df2[,dim(q_df2)[2]])))

#
# Daniel Romero-Alvarez & Laura Jimenez, 2020