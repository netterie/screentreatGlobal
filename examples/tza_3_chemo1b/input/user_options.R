############################################################################
## File Name: user_options.R

## File Purpose: Specify inputs for run_file.R
## Author: Jeanette Birnbaum
## Date: 10/14/2014
## Edited on: 

## Additional Comments: 
#############################################################################

############################################################
# Establish model version if this file is not being called 
# by a wrapper
############################################################
# TODO
if (!'using_wrapper'%in%ls()) {
    warning('Empyting the workspace')
    rm(list=ls())
    
    model_version <- 'tza_3_chemo1b'

    setwd('~')
    if (grepl('jbirnbau', getwd())) rootdir <- getwd()
    if (grepl('jeanette', getwd())) rootdir <- file.path(getwd(), 'Documents', 'jbirnbau')
    base_path <- file.path(rootdir, 'screentreatGlobal/examples')
}

############################################################
# Simulation features
############################################################
country = 'tza'
nsim = 100
times = c(5,10)
pop_size = 100000
study_year = 2013 # Approx year of incidence/life table data
inc_source = 'globocan'
standard_pop = TRUE

############################################################
# Input data files
############################################################

treat_file = file.path(base_path, model_version, 'input', 'input.csv')
incidence_file = file.path(rootdir, 'screentreatGlobal/data', 
                           paste0(country, '_incidence.csv'))
library_file = file.path(rootdir, 'screentreatGlobal/code/screentreat_library.R')
life_table_file = file.path(rootdir, 'screentreatGlobal/data', 
                           paste0(country, '_lifetable.csv'))
age_file = file.path(rootdir, 'screentreatGlobal/data', 
                     paste0(country, '_age.csv'))
if (standard_pop) age_file = file.path(rootdir, 'screentreatGlobal/data', 
                     'std_age.csv')


############################################################
# Population features
############################################################

# 8/19/16 note: using a function format_age to get single-year
# ages from 5-yr age groups
if ('age_file'%in%ls()) {
    source(library_file)
    ages <- format_age(age_file, minAge=0, maxAge=80)
    pop_chars = 
        list(age=ages,
             male=data.frame(male=c(0), prop=c(1)))
} else {
    pop_chars = 
        list(age=data.frame(age=c(40), prop=c(1)),
             male=data.frame(male=c(0), prop=c(1)))
}

# Is age in the data age at clinical incidence? 
# If not, provide incidence table
age_is_ageclin = FALSE
if (!age_is_ageclin) {
    inc_table = read.csv(incidence_file, header=FALSE, 
                         stringsAsFactors=FALSE)
}

# Denominator for reporting results
denom <- 100000

############################################################
# Screening, treatment and cancer mortality
############################################################

# Stage shift
HR_advanced = 0.6/0.85

# Within stage effects
instage_screen_benefit_early = 1
instage_screen_benefit_advanced = 1

# Add lead time? Default is undefined or FALSE
# If true, add mean lead time in years
lead_time = FALSE
if (lead_time) lt_mean = (40/12)

# Treatment HRs and distributions by subgroup-stage
treat_chars = read.csv(treat_file, header=TRUE, 
                        stringsAsFactors=FALSE)

# Survival distribuion: exponential or weibull?
surv_distr = 'exponential'

# Baseline mortality rates and population proportions by
# subgroup-stages. Subgroup stages specified here must
# match those given in the scrtrt_file
control_notreat = data.frame(stage=c(rep('Early',2),
                                     rep('Advanced',2)),
                             subgroup=rep(c('ER+',
                                            'ER-'),2),
                             mortrate=c(rep(.0446,2),rep(0.21, 2)),
                             prop=c(0.045, 0.105, 0.255, 0.595)
                                # Early ER+, Early ER-, Adv ER+, Adv ER-
                                # Expected ER+/- is 30%/70%
                                # Rough from literature: early/adv is 15%/85%
                                #  if missings are missing at random (MAR)
                             )


############################################################
# Other-cause mortality
############################################################

ocd_HR = 1

############################################################
# Run model
############################################################

source(file.path(rootdir, '/screentreatGlobal/code/run_file.R'))

