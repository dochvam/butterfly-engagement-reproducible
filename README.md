# Logistical and preference bias in participatory science butterfly data

This repository contains all code needed to replicate the analyses performed in "Logistical and preference bias in participatory science butterfly data" by Goldstein, Stoudt, Lewthwaite, Shirey, Mendoza, and Guzman.

## Components

The main code to be executed is contained in the "code" folder and should be run in numeric order (00 through 03) with two helper scripts. 

All but one of the data inputs can be found in the "data" folder. The iNaturalist data is required but not published within this directory. It may be obtained from GBIF at [this link](https://doi.org/10.15468/dl.rhmxtn) and should be saved into the "data" folder. 

Note that, at the request of eButterfly, we have removed location information for eButterfly observations of sensitive species in this repository. Please contact eButterfly directly for the uncensored data.

## Other notes

The code in the "01" file takes a long time to run (all of the GAMs for each species). You can skip ahead to "02" to run meta-analysis code on the species-level results, which we provide in "output" as "estimated_indices_fromGAMs.csv".




