# Environmental Sampling
 
The Hutchinson’s duality allows linking a point in geography to a point in an n-dimensional environmental space where ecological niches are defined. We took advantage of this relationship to design future surveys of a species of interest based on the suitibility values from an Ecological Niche Model (ENM) or Species Distribution Model (SDM). Sampling in different environmental gradients has the advantage of obtaining information that was excluded for the model calibration phase due to the lack of occurrence points. 

Potential detections of the species of interest from different suitability categories, especially those that are different from the maximum value of suitability, will improve the fitting of ecological niche models by adding occurrences.

The examples included here focuses on the species Burkholderia pseudomallei (the Gram-negative bacteria known to cause the infectious disease called melioidosis) in the regions Ceará, Brazil and Texas, United States (this last one as presented in Romero-Alvarez et al., 2020). These examples illustrate how the proposed methodology of environmental sampling can be replicated for other states, countries, or regions. 

Finally, the functions available also allow the exploration of environmental spaces defined by decimal geographical coordinates and the development of sampling excercises using transects in this space and how they translate to geogprahies. 
 
In order to reproduce the examples, please follow these steps:

(1) Download or clone this repository in your computer and un-zip the folders that contain the datasets: 'shapefiles2.zip', 'categorized_models.zip', 'environmental_variables.zip', 'uncertainty_models.zip', and 'env_samp_points.zip'. Notice that you will need to replace these files in order to reproduce the example with another species, model, points or region of study. You will need to have folders that contain: the polygons of the regions of interest, categorized ENMs/SDMs outputs, environmental layers used in the ENMs/SDMs, rasters representing uncertainty or other predictors to be explored, and coordinates. 

(2) Open the R project called Environmental-sampling and from there, open the R script called Worked_Examples.R. This code will guide you to reproduce one analysis for the region of Ceará and another for Texas following the results of the paper accompanying these functions (see below). The code also allow to use the functions using only coordinates as it is represented in the third available example. 

(3) Run all the lines under each example and read the comments to understand what each line of code is doing. Detailed comments are included in the source codes of the functions, see the scripts called E-space-functions.R, Hutchinson-functions.R, and Post-track-functions.R for further detailed information.

We also recommend you to check out the paper in which we first presented and used this methodology:
- Romero-Alvarez et al. ( ) Data-oriented sampling of Burkholderia pseudomallei environmental distribution via ecological niche modeling: an example in Texas, United States. ADD DOI.

How to cite this code:
