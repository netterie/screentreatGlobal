#############################################################################
## File Name: wrapper.R

## File Purpose: Run screening-with-treatment model, saving input code and data so 
## results can be reproduced, and saving metadata about the model run
## Author: Jeanette Birnbaum
## Date: 6/12/2013
## Edited on:6/21/2013

## Additional Comments: 
## Process is as follows:
## 1. Edit [area folder, e.g. diagnostics]/code/user_options_master.r
## 2. Edit the "model run info" section of this file, 
##    model_run_wrapper.r, to describe the model chosen in #1
## 3. Open a new instance of R from here, and source this file. (This
##    is because the initial working directory is assumed to be
##    [area folder, e.g. diagnostics]/code)
#############################################################################

# Empty the workspace
rm(list=ls())

# EDIT THIS SECTION:
############################################################
# model run information
############################################################

# What is the name of the user options file associated with
# this wrapper file? Location is presumed to be in [area]/code,
# e.g. diagnostics/code
user_options_file = 'user_options.R'
input_data_file = 'input_tza_barrier-nobarrier.csv'

# Should a copy of the input dataset be stored in the model folder?
# It will be named input_data.csv, so describe it below in 
# model_description
copy_data = TRUE

# Establish model folders and describe the model
model_type = 'screentreatGlobal'
model_version = 'tza_3_chemo1cages50-69'
model_description = c(country='Tanzania',
                      data_stage='Catherine lit search approx',
                      data_surv='CT for early, Galukande 2015 for advanced',
                      data_treat='EBCTCG trials',
                      data_other='ER dist from Catherine',
                      design='two trials: withoutER, withER',
                      stages='early and advanced',
                      stageshift='none',
                      prognostic_factors='ER',
                      treatments='None, Tamoxifen, Chemo',
                      efficacies='from EBCTCG',
                      life_table='IHME',
                      other='like _3c, but ages 50-69')

# Do you want to actually run the model now, or 
# do everything in this file except run the model?
run = FALSE


# THE FOLLOWING SECTIONS SHOULD NOT NEED EDITING:
using_wrapper=TRUE
############################################################
# directory setup
############################################################

setwd('~')
if (grepl('jbirnbau', getwd())) rootdir <- getwd()
if (grepl('jeanette', getwd())) rootdir <- file.path(getwd(), 'Documents', 'jbirnbau')
if (grepl('Elan', getwd())) rootdir <- getwd()

base_path <- file.path(rootdir, model_type, 'examples')

############################################################
# create directories if this is a new model run
############################################################

# Look for directory to evaluate whether we need to create folders
create_dirs = tryCatch({dir.create(file.path(base_path, model_version))},
                       warning=function(w) FALSE)

if (create_dirs) {
    dir.create(file.path(base_path, model_version, 'input'))
    dir.create(file.path(base_path, model_version, 'output'))
}

############################################################
# record model description in a master file as well 
# as within this run's folder
############################################################

# Define a little function to do version tracking

track_version = function(version_info, filename) {
    if (file.exists(filename)) {
        tmp = read.csv(filename, header=TRUE)
        if (!version_info$version%in%tmp$version) {
            write.table(version_info, 
                        filename, 
                        sep=",",
                        append=TRUE, 
                        row.names=FALSE,
                        col.names=FALSE)
        }
    } else write.table(version_info, filename, sep=",", row.names=FALSE)
}

# Compile information
vinfo = data.frame(version=model_version,
                   t(data.frame(model_description)),
                   row.names=NULL)

# Save info in master version file
track_version(vinfo, 
              file.path(rootdir, model_type, 
                        'version_guide.csv'))

# Save within this run's folder
write.csv(t(vinfo), file.path(base_path,  model_version, 'model_info.csv'))

############################################################
# store code within this run's input folder
############################################################

file.copy(file.path(base_path, input_data_file),
          file.path(base_path, model_version, 'input', 'input.csv'),
          overwrite=TRUE)
file.copy(file.path(rootdir, model_type, 'code', 'wrapper.R'),
          file.path(base_path, model_version, 'input', 'wrapper.R'),
          overwrite=TRUE)
file.copy(file.path(rootdir, model_type, 'code', 'run_file.R'),
          file.path(base_path, model_version, 'input', 'run_file.R'),
          overwrite=TRUE)
file.copy(file.path(rootdir, model_type, 'code', user_options_file),
          file.path(base_path, model_version, 'input', 'user_options.R'),
          overwrite=TRUE)

############################################################
# run the model
############################################################

setwd(file.path(rootdir, model_type, 'code'))
if (run) source(user_options_file)

