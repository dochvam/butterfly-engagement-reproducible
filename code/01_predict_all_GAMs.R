library(mgcv)
library(tidyverse)
library(taxadb)
library(parallel)
library(tryCatchLog)
library(futile.logger)
library(R.utils)
library(nimble)


source("code/helper_fns.R")

# We create some directories. lpmtx will store some GAM outputs that we won't
# actually use but needed for reproducing the surfaces. ind_results is where
# we'll store the results of interest. Warnings is a folder where warnings and
# errors will be written (if they occur--they shouldn't if unless we're timing
# out). If an individual species throws an error, the error will be saved to a
# file and computation will continue
if (!dir.exists("intermediate/warnings")) dir.create("intermediate/warnings")
if (!dir.exists("intermediate/ind_results")) dir.create("intermediate/ind_results")
if (!dir.exists("intermediate/lpmtx")) dir.create("intermediate/lpmtx")

timeout <- 1e7 # Stop after a while for each dataset (Helps catch super slow species)

# Grab inputs
if (file.exists("intermediate/gridcts_east.csv")) {
  hex_data <- read_csv("intermediate/gridcts_east.csv") 
} else {
  stop("You need to first run the data creation file, data_code/createDataGrid.R")
}


output_file <- "intermediate/GAM_diffs.csv"

target_species_df <- read_csv("intermediate/target_species.csv")
target_species_df <- sample_n(target_species_df, size = nrow(target_species_df),
                              replace = FALSE)


# 1

for (k in 1:nrow(target_species_df)) {
  this_species <- target_species_df$species[k]
  this_region <- target_species_df$region[k]
  
  this_name_clean <- gsub(" ", "_", this_species)
  
  # Set up warning printing and output filepath
  # flog.appender(appender.file(paste0("intermediate/warnings/", this_name_clean, "_warning.log")))
  this_output_file <- paste0("intermediate/ind_results/res_", this_name_clean, 
                             "_", this_region, ".csv")
  
  # Only run if we haven't done this species yet
  if (!file.exists(this_output_file)) {
    cat("===== Processing:", this_species, "=====\n")
    
    # We do write a file and use the Rscript command to make sure that we're
    # running each model in a fresh R session. The idea here is to encourage
    # the garbage collector to do its job; before using this method we were
    # experiencing some hanging/memory overflow issues
    cat(paste0(c('
     library(mgcv)
     library(tidyverse)
     library(tryCatchLog)
     library(futile.logger)
     library(R.utils)

     hex_data <- read_csv("intermediate/gridcts_', this_region, '.csv")

     source("code/helper_fns.R")
     tryLog(process_spec_diff("',
                 this_species, '", "', this_region,
                 '", hex_data = hex_data, verbose = TRUE, timeout =',
                 timeout,
                 ', maxK = 15))'), collapse = ""), file = "tempfile.R")
    
    # sink("captured.txt")
    captured <- system("Rscript tempfile.R")
    sink()
  }
}


# Write the results (not hugely necessary)
result_df <- do.call(rbind, lapply(list.files("intermediate/ind_results", 
                                              full.names = T, pattern = "res_"), 
                                   read_csv))
write_csv(result_df, output_file)

predicted_df <- do.call(rbind, lapply(list.files("intermediate/ind_results", 
                                                 full.names = T, pattern = "predicted_"), 
                                      function(x) {
                                        temp <- read_csv(x)
                                        chunks <- str_split(x, "_")[[1]]
                                        temp$species <- chunks[4]
                                        temp$region <- substr(chunks[5], 1, 4)
                                        temp
                                      }))
write_csv(predicted_df, "intermediate/GAM_predictions.csv")

