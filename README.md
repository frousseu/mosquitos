**Code for mapping mosquito abundance data in southern Quebec**

The code is organized in several files as follows:

## *predata.R*

This file is for preparing some of the data that will be loaded in the following *data.R* file. Specifically, it is used for downloading the Daymet data and then for calculating tmean from tmin and tmax. A new set of tmean rasters is produced from this file. Once tmean is available, anomalies (anom) are extracted using daily averages built from GAM models over the whole period.

## *data.R*

This file is for downloading and preparing the data that is used in the analysis. Its loads the mosquito abundance data, the Daymet weather data and the land cover LULC data and it also contains some visualizations of the data. It also contains the building of the grid that will be used for mapping predictions. Hence, some choices have to be made in this section in order to determine what to map. This is done here instead of in the following *model_parameters.R* file in order to do the Daymet and LULC data extraction once. 

## *model_parameters.R*

This file is used for setting the parameters of the analysis. This includes building the model sets, checking VIFs, subsetting the data, choosing the species, transforming and scaling variables, building the mesh, the spde, the priors, the A matrix, the stacks, the index and the data that will be used for prediction graphs.

## *model_selection.R*

This file is for looping over the different models to get the DIC values. It only uses the stack that is made for estimation, but not the parts that are meant for mapping predictions. The file is meant to be ran on an external server for greater computation power.

## *model_output.R*

This file runs the chosen model with the complete stack (estimation, mapping and graphical predictions) and samples from the posterior to get values that will be used in predictions and model checking.

## *results_display.R*

This file is used for model checking, validation and for producing all graphical or textual results.
