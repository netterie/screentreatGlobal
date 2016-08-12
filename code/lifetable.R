
############################################################################
## File Name: lifetable.R

## File Purpose: Extract and format country-specific life tables
## Author: Jeanette Birnbaum
## Date: 07/07/16
## Edited on: 

## Additional Comments: 
#############################################################################

## EDIT BEFORE RUNNING
############################################################
# Country and life table selection
#ctry <- c('United States')
#source <- c('ihme')

# Files
ihme_file <- c('~/Dropbox/BCModel/data/lifetables/lifetable_IHME_GBD_2013_LIFE_TABLE_1990_2013_Y2014M12D17/IHME_GBD_2013_LIFE_TABLE_1990_2013_Y2014M12D17.csv')
############################################################

############################################################
# Setup
############################################################
setwd('~/Documents/jbirnbau/screentreatGlobal')
source('code/screentreat_library.R')

############################################################
# Format and save lifetables for each country
############################################################

usa <- format_lifetable(ihme_file, 'United States')
uga <- format_lifetable(ihme_file, 'Uganda')

write.csv(uga, 'data/uga_lifetable.csv')

tza <- format_lifetable(ihme_file, 'Tanzania')

write.csv(tza, 'data/tza_lifetable.csv')


