#############################################################################
## File Name: run_file.R

## File Purpose: Main code for CANTRANce area (c) (screening)
## Author: Jeanette Birnbaum
## Date: 10/14/2014

## Additional Comments: 
#############################################################################

############################################################
# Setup
############################################################

source(library_file)
library(cantrance)
set.seed(98103)

############################################################
# Prepare inputs
############################################################
cat('\nPreparing inputs...')

# Identify the trial arms 
trials <- gsub('prop_','',
               grep('prop', colnames(treat_chars), value=TRUE))
if (length(trials)==0) stop('Fix trial names in treatment input file')

# Identify subgroups
subgroups <- as.character(unique(control_notreat$subgroup))

# Prepare numeric and text group IDs:
# First for stage-subgroups (SS)
control_notreat <- transform(control_notreat,
                             SSno=1:nrow(control_notreat),
                             SSid=paste(stage, subgroup, 
                                           sep='.'))
treat_chars <- transform(treat_chars,
                         SSid=paste(stage, subgroup, sep='.'))
                         
treat_chars <- merge(treat_chars, 
                     subset(control_notreat, 
                            select=-c(stage, subgroup)),
                     by='SSid',
                     all.x=TRUE, sort=FALSE)

# Now for treatment-stage-subgroups
treat_chars <- transform(treat_chars,
                         txSSid=paste(SSid,tx,sep='.'),
                         txSSno=1:nrow(treat_chars))

############################################################
# Prepare population of clinically incident cases
############################################################
cat('\nSimulating clinical incidence...')

# First, simulate the indepenent characteristics given in
# pop_chars. Right now it assumes that each element in the list
# refers to a single variable rather than a joint distribution
# of variables. Very easy to generalize to return the row #
# of the original dataset rather than the value of a single
# variable. It would be easy to instead use Leslie's 
# create_pop_list() function, and that would work with
# a more complex age pattern, too

pop_chars_rows <- lapply(pop_chars, function(x, Npop, Nsim) {
                        sim_multinom(nsims=Npop, 
                                     nreps=Nsim,
                                     probs=x$prop,
                                     names=1:nrow(x))
                        }, pop_size, nsim)

# Ages at entry 
# Since we are not using cohort life tables, the birth_year
# is roughly specified for convenience, to use the cohort 
# lifetable code
pop_chars[['age']] <- transform(pop_chars[['age']],
                                    birth_year = study_year-50)

ageentry <- return_value_from_id(ids=pop_chars_rows[['age']],
                                 df=pop_chars[['age']],
                                 value='age')

# Ages at other-cause death (load alternative life table if desired)
if ('life_table_file'%in%ls()) {
    if (life_table_file!='') life_table <- read.csv(life_table_file)
    life_table$BirthCohort <- study_year-50
}

ageOC <- calc_ac_lifespan_pop(popdata=pop_chars,
                                bootrows=pop_chars_rows,
                                results_as_matrix=TRUE, 
                                survHR=ocd_HR) 

# Age at clinical incidence
if (!age_is_ageclin) {

    # Format the incidence data
    inc_data <- format_clinical_incidence(incidence_file,
                                          source=inc_source)

    # Simulate
    ageclin <- sim_clinical_incidence(popdata=pop_chars,
                               bootrows=pop_chars_rows,
                               incidence=inc_data,
                               results_as_matrix=TRUE)
} else {
    ageclin <- ageentry
}

############################################################
# Simulate subgroup and stage
############################################################
cat('\nSimulating subgroup, stage, and stage shift...')

# For control arms (for screening arms, stage will change)
# Numbers refer to rows of control_notreat
control_notreat_rows <- sim_multinom(nsims=pop_size, 
                              nreps=nsim, 
                              probs=control_notreat$prop, 
                              names=1:nrow(control_notreat))

############################################################
# Apply stage shift within subgroups
############################################################

# We can still refer to control_notreat, but now we will 
# have a distribution of rows shifted towards the early 
# stages

# First, identify Early-Advanced stage row pairs within
# subgroups
stage_pairs <- sapply(subgroups, function(x, df) which(df$subgroup==x), 
                      control_notreat)
dimnames(stage_pairs) <- list(control_notreat$stage[stage_pairs[,1]],
                              subgroups)

# Now, for each subgroup, use the HR_advanced stage-shift 
# parameter to determine whether a person gets shifted 
# (1=yes, 0=no). We will generate this for everyone, but we
# will only use it for the advanced stage cases.
shift <- matrix(rbinom(nsim*pop_size, 1, prob=1-HR_advanced),
                nrow=pop_size, ncol=nsim)

# Shift stage
screen_notreat_rows <- control_notreat_rows
for (s in subgroups) {
    adv_cases <- screen_notreat_rows==stage_pairs['Advanced',s]
    screen_notreat_rows[adv_cases & shift==1] <- 
        stage_pairs['Early',s]
    #Check result: table(screen_notreat_rows[adv_cases])
}

############################################################
# Simulate treatment received in control and screening
# arms
############################################################
cat('\nSimulating treatment received...')

# For control arms first
# For each trial, simulate treatment (referring to treat_chars)
# by stage-subgroups
control_treatments <- sapply_withnames(trials, funX=function(
                          x, 
                          control_notreat_rows,
                          treat_chars,
                          pop_size,
                          nsim) {
                            # Define which trial's proportions to use
                            thisprop <- paste('prop', x, sep='_')
                            treat_results <- 
                                sim_treatment_by_subgroup(treat_chars,
                                                          control_notreat_rows,
                                                          thisprop,
                                                          pop_size,
                                                          nsim)
                         }, 
                         control_notreat_rows, 
                         treat_chars,
                         pop_size, 
                         nsim)

# For screening arms, we need to change treatment only for those
# who got stage shifted, so shift==1 and stage==Advanced. 

# First create indicator of needing to change treatment
shift_treatment <- shift==1 & 
                   control_notreat_rows%in%stage_pairs['Advanced',]

# Now get new treatment assignments for all early-stage cases
screen_treatments <- sapply_withnames(trials, funX=function(
                          x, 
                          screen_notreat_rows,
                          treat_chars,
                          pop_size,
                          nsim) {
                            # Define which trial's proportions to use
                            thisprop <- paste('prop', x, sep='_')
                            treat_results <- 
                                sim_treatment_by_subgroup(treat_chars,
                                                          screen_notreat_rows,
                                                          thisprop,
                                                          pop_size,
                                                          nsim)
                         }, 
                         screen_notreat_rows, 
                         subset(treat_chars, stage=='Early'),
                         pop_size, 
                         nsim)

# Replace screen_treatments with control_treatments for non-shifted 
# early-stage cases
for (t in trials) {
    screen_treatments[[t]][!shift_treatment] <- 
        control_treatments[[t]][!shift_treatment]
}


############################################################
# Simulate time from ageclin to cancer death 
# according to stage-subgroup and treatment
############################################################
cat('\nSimulating mortality...')

if (surv_distr=='exponential') {
    # Baseline mortality rates by stage-subgroup
    control_baserate <- return_value_from_id(id=control_notreat_rows, 
                                         df=control_notreat,
                                         value='mortrate')
    screen_baserate <- return_value_from_id(id=screen_notreat_rows, 
                                         df=control_notreat,
                                         value='mortrate')
} else if (surv_distr=='weibull') {

    # Helper function: given paramaters for a weibull distribution and a hazard ratio, 
    # outputs new scale parameter incorporating the hazard ratio
    weibRRcalc <- function(shape, scale, RR){
      return(scale/(RR^(1/shape)))
    }

    #Baseline mortality weibull paramaters by stage-subgroup
    control_baseshape <- return_value_from_id(id=control_notreat_rows, 
                                              df=control_notreat,
                                              value='mortshape')
    screen_baseshape <- return_value_from_id(id=screen_notreat_rows, 
                                             df=control_notreat,
                                             value='mortshape')
    control_basescale <- return_value_from_id(id=control_notreat_rows, 
                                              df=control_notreat,
                                              value='mortscale')
    screen_basescale <- return_value_from_id(id=screen_notreat_rows, 
                                             df=control_notreat,
                                             value='mortscale')
}

# Treatment HRs
control_HRs <- sapply_withnames(control_treatments, 
                                funX=return_value_from_id, 
                                treat_chars, 
                                'HR')
screen_HRs <- sapply_withnames(screen_treatments, 
                               funX=return_value_from_id, 
                               treat_chars, 
                               'HR')

if (surv_distr=='exponential') {
    # Final mortality rate
    control_rate <- sapply_withnames(control_HRs, 
                                     funX=function(x, rate) { x*rate }, 
                                     control_baserate)
    screen_rate <- sapply_withnames(screen_HRs, 
                                    funX=function(x, rate) { x*rate }, 
                                    screen_baserate)
} else if (surv_distr=='weibull') {
    # Final mortality scale parameter
    control_scale <- sapply_withnames(control_HRs, 
                                      funX=function(x, shape, scale){weibRRcalc(shape, scale, x)}, 
                                      control_baseshape,
                                      control_basescale)
    screen_scale <- sapply_withnames(screen_HRs, 
                                     funX=function(x, shape, scale){weibRRcalc(shape, scale, x)}, 
                                     screen_baseshape,
                                     screen_basescale)
}

# Simulate time from ageclin to cancer death for the 1st control group
if (surv_distr=='exponential') {
    clin2cd <- matrix(rexp(n=rep(1,pop_size*nsim), rate=control_rate[[1]]),
               nrow=pop_size, ncol=nsim)
} else if (surv_distr=='weibull') {
    clin2cd <- matrix(rweibull(n=rep(1,pop_size*nsim), shape=control_baseshape[[1]], scale=control_scale[[1]]),
                      nrow=pop_size, ncol=nsim)
}

# Now for other groups
# For the first control group, we will just get back the same times,
# i.e. ttcd==control_ttcd[[1]]
if (surv_distr=='exponential') {
    control_clin2cd <- sapply_withnames(control_rate, 
                                     funX=sim_same_qexp, 
                                     oldtime=clin2cd, 
                                     oldrate=control_rate[[1]], 
                                     prefix='control')
    screen_clin2cd <- sapply_withnames(screen_rate, 
                                    funX=sim_same_qexp, 
                                    oldtime=clin2cd, 
                                    oldrate=control_rate[[1]], 
                                    prefix='screen')
} else if (surv_distr=='weibull') {                      
    sim_same_qweib <- function(oldtime, oldscale, oldshape, newscale, prefix){
      newtime = qweibull(p = pweibull(oldtime, shape = oldshape, scale=oldscale), shape = oldshape, scale = newscale)
      colnames(newtime) = paste0(prefix, 1:ncol(newtime))
      newtime[abs(oldtime - newtime) < 0.001] = oldtime[abs(oldtime - newtime) < 0.001]
      return(newtime)
    }

    control_clin2cd <- sapply_withnames(control_scale, 
                                     funX=sim_same_qweib, 
                                     oldtime=clin2cd, 
                                     oldscale=control_scale[[1]], 
                                     oldshape=control_baseshape,
                                     prefix='control')
    screen_clin2cd <- sapply_withnames(screen_scale, 
                                    funX=sim_same_qweib, 
                                    oldtime=clin2cd, 
                                    oldscale=control_scale[[1]], 
                                    oldshape=control_baseshape,
                                    prefix='screen')
}
# For shifted cases, add lead time?
if ('lead_time'%in%ls()) {
if (lead_time) {
    cat('\nAdding lead times...')
    lts <- matrix(rexp(n=rep(1,pop_size*nsim), rate=1/lt_mean),
               nrow=pop_size, ncol=nsim)
    # We can again use the shift_treatment indicator
    # to determine who needs a lead time
    for (i in 1:length(trials)) {
        screen_clin2cd[[i]][shift_treatment] <- 
            screen_clin2cd[[i]][shift_treatment] + 
            lts[shift_treatment]
    }
}
}


############################################################
# Tally age at cancer death, cause of death, and
# time from trial start to mortality
############################################################
cat('\nTabulating results...')

# Age at cancer death
control_ageCD <- sapply_withnames(control_clin2cd, 
                                  funX=function(x, ageclin) { x + ageclin }, 
                                  ageclin)
screen_ageCD <- sapply_withnames(screen_clin2cd, 
                                 funX=function(x, ageclin) {x + ageclin }, 
                                 ageclin)

# Time from study start to cancer death
control_ttcd <- sapply_withnames(control_ageCD, 
                                  funX=function(x, ageentry) { x - ageentry}, 
                                  ageentry)
screen_ttcd <- sapply_withnames(screen_ageCD, 
                                 funX=function(x, ageentry) {x - ageentry}, 
                                 ageentry)

# Cause of death
control_CoD <- sapply_withnames(control_ageCD,
                                funX=function(x, ageOC) {
                                    ifelse(x<ageOC,1,0)
                                },
                                ageOC)
screen_CoD <- sapply_withnames(screen_ageCD,
                                funX=function(x, ageOC) {
                                    ifelse(x<ageOC,1,0)
                                },
                                ageOC)

# Time from trial start to all-cause death
control_ttd <- sapply_withnames(control_ageCD,
                                funX=function(x, ageOC, ageentry) {
                                    ifelse(x<ageOC, x-ageentry, ageOC-ageentry)
                                },
                                ageOC, ageentry)
screen_ttd <- sapply_withnames(screen_ageCD,
                                funX=function(x, ageOC, ageentry) {
                                    ifelse(x<ageOC, x-ageentry, ageOC-ageentry)
                                },
                                ageOC, ageentry)

############################################################
# Summarize case structure
############################################################

ids <- sapply_withnames(control_treatments,
                        funX=return_value_from_id,
                        treat_chars, 'txSSid')
names(times) <- as.character(times)
cases <- sapply_withnames(times, funX=function(x) {
                    cases <- ageclin<=(ageentry+x)
                    df <- ldply(sapply_withnames(ids, 
                                       funX=function(y, c) {
                                           round(table(y[c])/ncol(c))
                                       }, cases), 
                            rbind)
                    df$times <- x
                    return(df)
                    #round(table(SSids[cases])/ncol(cases))
                                })
cases <- ldply(cases, rbind)
cases[is.na(cases)] <- 0

mcases <- melt(cases, id.vars=c('.id', 'times'))
colnames(mcases)[which(colnames(mcases)=='.id')] <- 'trial'
mcases <- data.frame(mcases, 
                     do.call(rbind, strsplit(as.character(mcases$variable), '\\.')))
colnames(mcases) <- c('Trial', 'Time', 'id', 'N', 'Stage', 'Subgroup', 'Treatment')

rcases <- cast(mcases, Time+Trial+Treatment~Stage+Subgroup, drop='id', value='N')
rcases$Note <- ''
rcases$Note[1]  <- paste('Numbers reflect population size of',
                     nrow(ageclin),
                     'women')

# Save
write.csv(rcases,
          file.path(base_path, model_version, 'output', 
                    'cases.csv'),
          row.names=FALSE)

############################################################
# Summarize mortality across arms and trials
############################################################
cat('\nConstructing results tables...')

# New results - 8/20/16

    if(!'denom'%in%ls()) denom=10000

    control_years <- tally_years_simple(times, 
                                      control_ttd,
                                      per=denom)
    screen_years <- tally_years_simple(times, 
                                      screen_ttd,
                                      per=denom)

# New results - 7/13/15

    if(!'denom'%in%ls()) denom=10000

    control_cuminc <- tally_cuminc_simple(times, 
                                      control_ttcd,
                                      control_CoD,
                                      per=denom)
    screen_cuminc <- tally_cuminc_simple(times, 
                                      screen_ttcd,
                                      screen_CoD,
                                      per=denom)

# New results - 10/3/16 - survival among clinically incident
# Use these to make the graph

    # Average cumulative inc across sims, per denom
    control_cumincTmp <- llply(control_cuminc, .fun=function(x){colSums(x)/nrow(x)})
    screen_cumincTmp <- llply(screen_cuminc, .fun=function(x){colSums(x)/nrow(x)})
    # Compile into a data frame
    survamongIncC <- transform(ldply(control_cumincTmp, rbind),Group='Control')
    survamongIncS <- transform(ldply(screen_cumincTmp, rbind),Group='Screen')
    survamongInc <- rbind(survamongIncC,survamongIncS)
    survamongInc <- cast(melt(survamongInc), .id+variable~Group)
    survamongInc <- rename(survamongInc, c('.id'='Time', 'variable'='Trial'))
    # Now add the average # of incident cases across sims, per denom
    incCases <- ddply(rcases,.(Time,Trial),function(x) {
                          thesecols <- 
                              !colnames(x)%in%c('Time', 'Trial', 'Treatment', 'Note')
                          sum(rowSums(x[,thesecols],na.rm=TRUE))
                                      })
    incCases <- rename(incCases, c('V1'='Incidence'))
    incCases <- transform(incCases, Incidence=Incidence*(denom/nrow(ageclin)))
    survInc <- merge(survamongInc, incCases, all=TRUE)
    survInc <- transform(survInc,
                         Control=round(100*Control/Incidence),
                         Screen=round(100*Screen/Incidence))
    survInc <- subset(melt(survInc), variable!='Incidence')
    survInc$Note <- ''
    survInc$Note[1] <- 'Out of 100 incident cases, number surviving'

# Save
write.csv(survInc,
          file.path(base_path, model_version, 'output', 
                    'simpleSurvival.csv'),
          row.names=FALSE)

# Construct rows of the table

    # Cumulative incidence
    r1 <- lapply(control_cuminc, 
                              summarize_over_sims, 
                              funX='mean',
                              onecell=TRUE,
                              numdec=0)
    r2 <- lapply(screen_cuminc, 
                            summarize_over_sims, 
                            funX='mean',
                            onecell=TRUE,
                            numdec=0)

    # Within-trial MRR
    wtmrr <- sapply_withnames(names(control_cuminc),
                              funX=function(x) {
                                screen_cuminc[[x]]/
                                  control_cuminc[[x]]
                              })
    r3 <- lapply(wtmrr,                             
                 summarize_over_sims, 
                 funX='mean',
                 onecell=TRUE,
                 numdec=2)

    # Across-trial MRR, 
    atmrr_noscreen <- sapply_withnames(names(control_cuminc),
                                       funX=function(x){
                                         control_cuminc[[x]]/
                                           control_cuminc[[x]][,1]
                                       })
    atmrr_screen <- sapply_withnames(names(control_cuminc),
                                       funX=function(x){
                                         screen_cuminc[[x]]/
                                           control_cuminc[[x]][,1]
                                       })
    r4 <- lapply(atmrr_noscreen,                             
                 summarize_over_sims, 
                 funX='mean',
                 onecell=TRUE,
                 numdec=2)
    r5 <- lapply(atmrr_screen,                             
                 summarize_over_sims, 
                 funX='mean',
                 onecell=TRUE,
                 numdec=2)

    # Within-trial ARR
    wtarr <- sapply_withnames(names(control_cuminc),
                              funX=function(x) {
                                control_cuminc[[x]]-
                                  screen_cuminc[[x]]
                              })
    r6 <- lapply(wtarr,                             
                 summarize_over_sims, 
                 funX='mean',
                 onecell=TRUE,
                 numdec=1)
    
    # Across-trial MRR, 
    atarr_noscreen <- sapply_withnames(names(control_cuminc),
                                       funX=function(x){
                                         replicate(length(trials),
                                                   control_cuminc[[x]][,1]) -
                                         control_cuminc[[x]]
                                       })
    atarr_screen <- sapply_withnames(names(control_cuminc),
                                     funX=function(x){
                                       replicate(length(trials),
                                                 control_cuminc[[x]][,1]) -
                                         screen_cuminc[[x]]
                                     })
    r7 <- lapply(atarr_noscreen,                             
                 summarize_over_sims, 
                 funX='mean',
                 onecell=TRUE,
                 numdec=1)
    r8 <- lapply(atarr_screen,                             
                 summarize_over_sims, 
                 funX='mean',
                 onecell=TRUE,
                 numdec=1)

    # Across-trial years of live saved, 
    atyears_noscreen <- sapply_withnames(names(control_years),
                                       funX=function(x){
                                         replicate(length(trials),
                                                   control_years[[x]][,1]) -
                                         control_years[[x]]
                                       })
    atyears_screen <- sapply_withnames(names(control_years),
                                     funX=function(x){
                                       replicate(length(trials),
                                                 control_years[[x]][,1]) -
                                         screen_years[[x]]
                                     })
    r9 <- lapply(atyears_noscreen,                             
                 summarize_over_sims, 
                 funX='mean',
                 onecell=TRUE,
                 numdec=1)
    r10 <- lapply(atyears_screen,                             
                 summarize_over_sims, 
                 funX='mean',
                 onecell=TRUE,
                 numdec=1)

# Compile across follow-up times

new_table <- lapply(names(control_cuminc), 
                function(x) {
                  empty_row = rep('', length(trials))
                  tab = data.frame(rbind(empty_row,
                        r1[[x]],r2[[x]],r3[[x]],
                        empty_row,
                        r4[[x]],r5[[x]],r6[[x]],
                        empty_row,
                        r7[[x]], r8[[x]],
                        empty_row,
                        r9[[x]], r10[[x]]
                        ),row.names=NULL)
                  tab = data.frame(`Follow-up`=c(as.character(x),
                                               rep('', nrow(tab)-1)),
                                   Measure=c('Cumulative breast cancer mortality',
                                             '', '', 
                                             'MRRs within trials',
                                             'MRRs across trials', '', '', 
                                             'ARRs within trials', 
                                             'ARRs across trials', '', '',
                                             'Years saved across trials', '', ''),
                                   SubMeasure=c('','No screening', 'Screening', 
                                                '', '', 'No screening', 'Screening', 
                                                '', '', 'No screening', 'Screening',
                                                '', 'No screening', 'Screening'),
                                   tab,
                                   check.names=FALSE)
                })

new_table_full <- do.call('rbind', new_table)

new_table_full$Note <- paste('Note: Results per', denom, 'women')

# Save
write.csv(new_table_full,
          file.path(base_path, model_version, 'output', 
                    'cuminc_mrr_newtable.csv'),
          row.names=FALSE)


# Old results - attempting to save cuminc plot

if (1==0) {
    ############################################################
    # Graph cumulative incidence
    ############################################################
    cat('\nConstructing results graphs...')
    
    cuminc_table <- subset(final_table, 
                           grepl('Cumulative Incidence', final_table$Measure) &
                           Effect=='Effect of screening, same treatment')
    cuminc_table <- transform(cuminc_table, check.names=FALSE,
                              Arm = ifelse(grepl('Group 0', Measure), 
                                           'Control', 
                                           'Screening'))
    
    # A hack
    if (model_version=='breast_ER-HER_5') {
        cuminc_table$Trial[cuminc_table$Trial=='Contemp1999'] <- '1999'
        cuminc_table <- transform(cuminc_table, check.names=FALSE,
                                  Trial=factor(Trial, 
                                               levels=c('Historical', 
                                                        '1999', 
                                                        'Perfect'),
                                               labels=c('Historical', 
                                                        '1999', 
                                                        'Perfect')))
    }
    
    cuminc_plot <- ggplot(cuminc_table,
                          aes(x=`Follow-Up Year`,
                              y=Estimate,
                              colour=Trial,
                              shape=Arm)) + 
                   geom_line() +
                   geom_point(size=3) + 
                   scale_y_continuous(name='Cumulative Incidence')
    
    ggsave(plot=cuminc_plot,
           filename=file.path(base_path, model_version, 'output',
                              'cuminc_plot.pdf'),
           width=6, height=5)
} # end commenting out old results


############################################################
# Validate model: Age-specific cancer mortality rates,
# incidence rate and prevalence percent
############################################################
cat('\nComputing statistics for validation...')

# Return age-specific mortality rates for each arm and trial
# We'll ignore confidence intervals

control_mortrates <- sapply_withnames(trials,
                                   funX=function(x,
                                                 ageentry,
                                                 control_ttd,
                                                 control_ageCD) {
                                       data.frame(
                                           Trial=x,
                                           Arm='Control',
                                           age_specific_eventrate(
                                                     ageentry, 
                                                     control_ttd[[x]],
                                                     control_ageCD[[x]]))
                                   }, 
                                   ageentry,
                                   control_ttd,
                                   control_ageCD)
screen_mortrates <- sapply_withnames(trials,
                                   funX=function(x,
                                                 ageentry,
                                                 screen_ttd,
                                                 screen_ageCD) {
                                       data.frame(
                                           Trial=x,
                                           Arm='Screen',
                                           age_specific_eventrate(
                                                     ageentry, 
                                                     screen_ttd[[x]],
                                                     screen_ageCD[[x]]))
                                   }, 
                                   ageentry,
                                   screen_ttd,
                                   screen_ageCD)

historic_incidence <- data.frame(Trial='Historical', 
                             Arm='Incidence', 
                             age_specific_eventrate(ageentry, 
                                                    control_ttd[[1]],
                                                    ageclin,
                                                    equals=FALSE))
historic_prevalence <- data.frame(Trial='Historical', 
                             Arm='Prevalence', 
                             age_specific_eventrate(ageentry, 
                                                    control_ttd[[1]],
                                                    ageclin,
                                                    equals=FALSE,
                                                    prevalence=TRUE))
mortrates <- do.call('rbind', c(control_mortrates, screen_mortrates))
mortrates <- rbind(mortrates, historic_incidence, historic_prevalence)
mortrates_wide <- cast(mortrates, Trial + Age.Group ~ Arm, value='Rate')

# Save
write.csv(mortrates_wide,
          file.path(base_path, model_version, 'output', 
                    'cancer_mort_per100K.csv'),
          row.names=FALSE)

############################################################
# Validate model: Mean survivals
############################################################

if (1==1) {
    # Baseline mean survivals
    if (surv_distr=='exponential') {
        base_msurvs <- subset(
                              transform(control_notreat, 
                                        MeanSurv=1/mortrate,
                                        Arm=SSid,
                                        Subgroup=c('Baseline survivals',
                                                   rep('',
                                                       nrow(control_notreat)-1))),
                              select=c(Subgroup, Arm, MeanSurv))
    } else if (surv_distr=='weibull') {
        stop('Not coded')
    }
    # Indicator of stage-subgroup
    control_stage <- return_value_from_id(control_notreat_rows,
                                          control_notreat,
                                          'SSid')
    screen_stage <- return_value_from_id(screen_notreat_rows,
                                          control_notreat,
                                          'SSid')

    # Mean net survivals - Control
    msurvs <- vector('list')
    for (s in control_notreat$SSid) {
        msurvs[[s]] <- sapply_withnames(control_clin2cd,
                                             mean_matrix_subgroup,
                                             control_stage==s)
    }
    # Mean net survivals - Screen
    msurvsScr <- vector('list')
    for (s in control_notreat$SSid) {
        msurvsScr[[s]] <- sapply_withnames(screen_clin2cd,
                                             mean_matrix_subgroup,
                                             screen_stage==s)
    }

    msurv_df <- data.frame(do.call('rbind', msurvs))
    msurvScr_df <- data.frame(do.call('rbind', msurvsScr))
    msurv_table <- data.frame(base_msurvs[,2:3],
                         msurv_df, msurvScr_df)
    colnames(msurv_table) <- c('Subgroup', 'Expected Mean',
                               paste0('Control_', trials),
                               paste0('Screen_', trials))
    for (i in 2:ncol(msurv_table)) {
        msurv_table[,i] <- unlist(msurv_table[,i])
    }
                         

    # Save
    write.csv(cbind(msurv_table),
              file.path(base_path, model_version, 'output', 
                        'mean_survivals.csv'),
              row.names=FALSE)

}



