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
    model_version <- 'usa_but_uga_stage-surv'
    base_path <- '~/Documents/jbirnbau/screentreatGlobal/examples'
}

############################################################
# Input data files
############################################################

treat_file = file.path(base_path, model_version, 'input', 'input.csv')
incidence_file = file.path(rootdir, '/screentreat/data/bc_1975-1979_incidence.csv')
library_file = file.path(rootdir, '/screentreat/code/screentreat_library.R')
#life_table_file = '~/cantrance/diagnostics/data/life_table_Reed2012.csv'

############################################################
# Simulation features
############################################################
nsim = 10
times = c(10,25)
pop_size = 100000
study_year = 2000

############################################################
# Population features
############################################################

pop_chars = 
    list(age=data.frame(age=c(50), prop=c(1)),
         male=data.frame(male=c(0), prop=c(1)))

# Is age in the data age at clinical incidence? 
# If not, provide incidence table
age_is_ageclin = FALSE
if (!age_is_ageclin) {
    inc_table = read.csv(incidence_file, header=TRUE, 
                         stringsAsFactors=FALSE)
}

############################################################
# Screening, treatment and cancer mortality
############################################################

# Stage shift
HR_advanced = 0.85

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
control_notreat = data.frame(stage=c(rep('Early',4),
                                     rep('Advanced',4)),
                             subgroup=rep(c('ER+HER2+',
                                            'ER+HER2-',
                                            'ER-HER2+',
                                            'ER-HER2-'),2),
                             mortrate=c(rep(.01992,4),rep(0.21, 4)),
                             prop=c(0.11*2*c(0.04, 0.38, 0.02, 0.06),
                                    0.89*2*c(0.06, 0.34, 0.03, 0.07))
                             )


############################################################
# Other-cause mortality
############################################################

ocd_HR = 1

############################################################
# Run model
############################################################

source(file.path(rootdir, '/screentreat/code/run_file.R'))

