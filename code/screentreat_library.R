 #############################################################################
## File Name: screentreat_library.R 

## File Purpose: Contains functions unique to the screening + treatment
## project
## Author: Jeanette Birnbaum
## Date: 5/2/2013
## Edited on: 9/9/2014

## Additional Comments: 
#############################################################################

############################################################
## sim_multinom
############################################################

sim_multinom <- function(
    nsims,
        ### Number of size=1 sims from the multinomial
    nreps,  
        ### Number of times to replicate the nsims
        ### Put the smaller number here, the bigger one as nsims
    probs,
        ### Vector of probabilities for each of the K categories
    names
        ### Vector of names associated with each of the categories
        ### Must be of same length as probs
        ### Can be character or numeric. Use numeric row IDs if 
        ### you're referring to a group defined by multiple
        ### variables whose groupings are defined in a separate
        ### data frame
) {

    # Helper function to do one rep
    one_rep <- function(id, nsims, probs, names) {

        # Do the draw from rmultinom
        draws <- rmultinom(nsims, size=1, probs)

        # Match to names
        indices <- apply(draws, 2, function(x) which(x==1))
        return(t(names[indices]))
    }

    # Now do the replicates
    all_reps <- sapply(1:nreps, FUN=one_rep, nsims, probs, names)

    return(all_reps)
}

############################################################
## format_age
############################################################
# The gsize variable indicates the age group intervals

format_age <- function(age_file, gsize=5, minAge=20, maxAge=85) {

    grouped <- read.csv(age_file, header=TRUE)

    if (!is.numeric(grouped$Population)) stop('Error in format_age: character pop sizes')

    # Split into single-year age groups
    split <- lapply(grouped$Age, function(a,data=grouped,g=gsize) {
                        NewPop <- data[which(data$Age==a), 'Population']/g
                        data.frame(age=a:(a+4), pop=rep(NewPop,g))
                         })
    split <- do.call('rbind', split)
    
    # Verify
    if (sum(split$pop)!=sum(grouped$Population)) stop('Division error in format_age')

    # Return proportions after limiting the minimum age
    split <- subset(split, age>=minAge & age<=maxAge)
    split <- transform(split, prop=pop/sum(pop))

    return(split)

}

############################################################
## format_lifetable
############################################################

format_lifetable <- function(
    lifefile,
    ctry,
        ### Match 'location_name' field for IHME
    source='ihme'
) {

    # Read in file
    lt <- read.csv(lifefile)
    
    # Conditional survival (SurvY)
    switch(source,
           ihme = {
                ltc <- subset(lt, location_name==ctry)
                ltc <- subset(ltc, year==max(ltc$year) & metric=='Deaths' & 
                                   age_group_name!='All Ages' & 
                                   sex_name!='Both sexes',
                              select=c(age_group_name, sex_id, year, mean))
                ltc <- within(ltc, {
                                  Age=gsub(" to [0-9]*", '', age_group_name)
                                  Age[Age=="<1 year"] <- '0'
                                  Age[Age=="80 plus"] <- '80'
                                  Age=as.numeric(Age)
                                  Male=sex_id
                                  Male[sex_id==2] <- 0
                                  TotQ=mean/100000
                                  SurvY=exp(-TotQ*5) #5-year intervals!
                                  Year=year
                              })
            })

    # Order
    ltc <- ltc[order(ltc$Male, ltc$Age),]

    # Single years
    age_max <- max(ltc$Age)
    mortM <- data.frame(Age=0:age_max, TotQ = with(subset(ltc,Male==1), 
                                                         approx(Age, TotQ, 
                                                                xout=0:age_max))$y)
    mortF <- data.frame(Age=0:age_max, TotQ = with(subset(ltc,Male==0), 
                                                         approx(Age, TotQ, 
                                                                xout=0:age_max))$y)
    # Cumulative survival (Survival)
    #ltc <- ltc[order(ltc$Male, ltc$Age),]
    ltl <- data.frame(Male=c(rep(0,nrow(mortF)),rep(1,nrow(mortM))),
                      rbind(mortF, mortM))
    ltl <- transform(ltl, 
                     SurvY=exp(-TotQ),
                     Year=max(ltc$year))

    seqProd <- function(x) sapply(1:length(x), FUN=function(y) prod(x[1:y]))
    ltl$Survival <- c(seqProd(subset(ltl, Male==0)$SurvY),
                      seqProd(subset(ltl, Male==1)$SurvY))

    # All die at age 100 - how to do? Handle it on the simulation side for now

    return(subset(ltl, select=c(Year, Male, Age, Survival)))
}

############################################################
## format_clinical_incidence
############################################################

# Assumes the cantrance package is loaded

format_clinical_incidence <- function(
    incfile,
        ### Path to SEER clinical incidence table
    source='seer'
        ### Source: seer, globocan
) {

    switch(source, 
           seer = {
                # Read in data
                incdat <- read.csv(incfile, stringsAsFactors=FALSE)
                incdat <- within(incdat, {
                               Age <- gsub('[0-9]+=', '', Ages)
                               Age <- as.character(gsub(' years', '', Age))
                               Age[Age=='85+'] <- 87
                               Age <- suppressWarnings(sapply(strsplit(Age, '-'),
                                             function(x) mean(as.numeric(x))))
                        })
                incdat <- incdat[!is.nan(incdat$Age),]
                incdat <- incdat[!is.na(incdat$Age),]

           },
           globocan = {
                # Read in data
                incdat <- read.csv(incfile, stringsAsFactors=FALSE, header=FALSE)
                colnames(incdat)  <- c('sex', 'site', 'ageGroup', 
                                       'cases', 'persYears', 'source')
                # Select female breast cancer, construct Crude.Rate and 
                # Age (transform 1-19 into age lower bound; ignore 19=unknown age)
                incdat <- subset(incdat, site==113 & sex==2 & ageGroup<19)
                incdat <- within(incdat, {
                                    Crude.Rate=100000*(cases/persYears)
                                    Age=(ageGroup-1)*5
                                    Age[Age==85] <- 87
                                })

                # Simply reconcile cases but no person-years in older ages:
                # hold the rate constant at the last observed rate
                zero <- which(incdat$persYears==0)
                incdat[zero,'Crude.Rate'] <- incdat[min(zero)-1,'Crude.Rate']
                warning('Graph age-specific incidence to check assumptions')

           }
        )

    age_max <- max(incdat$Age)
    incidence <- 
        data.frame(age = 0:age_max, 
                   incidence_rate = with(incdat, 
                                          approx(Age, 
                                                 Crude.Rate, 
                                                 xout=0:age_max))$y/100000)

    incidence$incidencefree_rate <- 1-incidence$incidence_rate
    incidence_free_survival <- c()
    for(i in 1:(age_max+1)){
      incidence_free_survival[i] <- 
          prod(incidence$incidencefree_rate[1:i])
    }
    incidence$incidencefree_survival <- incidence_free_survival
    rm(incidence_free_survival)

    return(incidence)
}

############################################################
## sim_clinical_incidence
############################################################

sim_clinical_incidence = function# Use sim_KM to efficiently simulate ages at clinical incidence for nsim populations

##description<< Given a population with a clinical incidence-free  
## curve and ages, returns random draws of their ages at incidence

(popdata, 
    ### Data frame where individuals are rows, with columns 
    ### birth_year, age, and male
bootrows,
    ### Matrix/df of row indicators that can be applied to
    ### the data to recover different bootstraps of the
    ### data. Each column is a different bootstrap of the
    ### data
incidence,
    ### Matrix/df with columns 'age' and 'incidencefree_survival'
results_as_matrix=FALSE
    ### Convert the results from a data frame to a matrix?
) {

    # Determine number of simulations
    nsim = ncol(bootrows)
    if (is.null(nsim)) nsim = ncol(bootrows[[1]])

    # Estimate for each simulation separately
    newpop = sapply(1:nsim,
        function(y) {
            # Construct the dataset
            if (!is.data.frame(popdata)) {
                thisboot = subset(build_simpop(datalist=popdata,
                                               rowlist=bootrows,
                                               sim=y),
                                  select=c(age))
            } else {
                thisboot = popdata[bootrows[, y], 
                                   c("age",
                                     "birth_year",
                                     "male")]
            }
            thisboot$tempid = 1:nrow(thisboot)
            
            # Simulate ages at incidence
            inc = ddply(thisboot,
                        .(age),
                        function(x) {
                            # Extract age
                            Age = min(x$age)

                            # Nsims
                            if (is.data.frame(x))
                                nrow=nrow(x) else
                                    nrow=1

                            # Max draw based on current age
                            maxu <- incidence[incidence$age==Age,
                                              'incidencefree_survival']
                            ageinc <- 
                                with(incidence, 
                                     sim_KM(incidencefree_survival,
                                            age,
                                            smalltimes=Age,
                                            bigtimes=121,
                                            nsims=nrow,
                                            maxdraw=maxu))
                            ageinc <- 
                                data.frame(x, 
                                           ageinc=matrix(ageinc,
                                                         nrow=nrow,
                                                         ncol=1))
                            return(ageinc)
                        })
            
            # Return to original sort order, and return 
            # ageOC
            inc = inc[order(inc$tempid), ]
            return(inc$ageinc)
        })
    
    # Return results as a matrix?
    if (results_as_matrix) {
        colnames(newpop) = 1:nsim
        newpop = as.matrix(newpop)
    }

    return(newpop)

### A data frame or matrix with nsim columns of randomly
### drawn ages at other-cause death for individuals (rows)
}

############################################################
## sim_treatment_by_subgroup
############################################################

sim_treatment_by_subgroup = function# Simulate treatment by subgroups

##description<< Given the proportions of treatments by subgroups, 
## simulate treatment

(treat_chars,
    ### Data frame indicating stage-subgroup #s in the SSno 
    ### variable and stage-subgroup-treatment #s in txSSno, and 
    ### prop_* columns where * is the trial name. 
stage_subgroup_rows,
    ### Matrix indicating, for each person-sim, which stage-subgroup
    ### they are in (using the SSid from treat_chars)
thisprop,
    ### Character indicating the prop_* column name from which
    ### to take treatment proportions
pop_size,
    ### Population size
nsim
    ### Number of sims
) {
    
    # Matrix for results
    results <- matrix(NA, ncol=nsim, nrow=pop_size)

    # Loop through stage-subgroups
    for (g in sort(unique(treat_chars$SSno))) { 
        # Subset to this stage-subgroup (SS)
        tempid <- subset(treat_chars, SSno==g)        
        # Simulate treatment for everyone
        # using this group's probs
        treat_choice <- 
            sim_multinom(nsims=pop_size, 
                         nreps=nsim,
                         probs=tempid[,thisprop],
                         names=tempid[,'txSSno'])
        # Assign values to those in group g
        results[stage_subgroup_rows==g] <-
            treat_choice[stage_subgroup_rows==g]
    }

    return(results)
}

############################################################
## return_value_from_id
############################################################

return_value_from_id = function# Use an ID no to look up a value

##description<< Using a matrix of IDs, look up a value in a dataframe

(ids,
    ### Matrix of ids
df,
    ### Data frame containing the values of interest, sorted
    ### by ID (i.e. row #s == ID #s)
value
    ### Column name of value to return from the df
) {

    return(apply(ids, 2, function(x) df[,value][x]))
    
}

############################################################
## sapply_withnames
############################################################

sapply_withnames = function# Use sapply a certain way

##description<< sapply with USE.NAMES=TRUE and simplify=FALSE

(list_object,
    ### List to sapply
funX,
    ### Function to pass to FUN=...
...
    ### Other arguments to pass to sapply
) {

    sapply(list_object,
           FUN=as.function(funX),
           ...,
           USE.NAMES=TRUE,
           simplify=FALSE)
    
}


############################################################
## summarize_over_sims
############################################################

summarize_over_sims = function# Return the mean and 95% interval over sims

##description<< Return the mean and 95% interval over sims. Sims are 
##represented as columns

(data,
    ### Matrix, with columns (default) or rows as different sims
 funX,
    ### Function to apply to columns, e.g. 'mean' or 'sum'
 ...,
    ### Optional arguments to mean and quantile, e.g. 'na.rm=TRUE'
 applyto=2,
    ### 2=columns, 1=rows
 onecell=FALSE,
    ### Old format - UI in separate column
 numdec=4
    ### Rounding
) {
    
    fixdec <- function(x, k) format(round(x, k), nsmall=k)
  
    estimates <- apply(data, applyto, funX)
  
    if (onecell)  {
      ci <- apply(data, applyto, quantile, probs=c(0.025, 0.975))
      returnvector <- paste0(fixdec(estimates, numdec), 
                             ' (',
                             paste(fixdec(ci[1,],numdec),
                                    fixdec(ci[2,], numdec),
                                   sep=','),
                             ')')
      names(returnvector) <- names(estimates)
      return(returnvector)
    } else {
      return(data.frame(Estimate=round(mean(estimates, ...),numdec),
                        UI=paste0('(',
                                  paste(round(quantile(estimates, 
                                                       probs=c(0.025, 0.975),
                                                       ...),numdec),
                                        collapse=','),
                                  ')')))
      
    }
}

############################################################
## tally_years_simple
############################################################

tally_years_simple = function# Return years of life lived
 
##description<< For specified follow-up times, return years lived
##patterned after tally_cuminc_simple
 
(followup,
    ### Vector of follow-up times
 etimes,
    ### List of matrices of times to event
 per=10000
    ### Denominator for cumulative incidence
) {
    years <- lapply(followup,
                     function(x) {
                       sapply(names(etimes), function(y) {
                         timesmat <- matrix(x,ncol=ncol(etimes[[y]]),
                                            nrow=nrow(etimes[[y]]))
                         cd_by_x <- pmin(timesmat, etimes[[y]])
                         (colSums(cd_by_x)/nrow(cd_by_x))*per
                       })
                     })  
    
    names(years) <- followup
    return(years)
}
 
############################################################
## tally_cuminc_simple
############################################################

tally_cuminc_simple = function# Return cumulative incidence
 
##description<< For specified follow-up times, return cumulative
##incidence of event across all trials and arms and sims
 
(followup,
    ### Vector of follow-up times
 etimes,
    ### List of matrices of times to event
 event,
    ### List of matrices of event indicators (1==yes, 0=no)
 per=10000
    ### Denominator for cumulative incidence
) {
    cuminc <- lapply(followup,
                     function(x) {
                       sapply(names(etimes), function(y) {
                         cd_by_x <- etimes[[y]]<x & event[[y]]==1
                         (colSums(cd_by_x)/nrow(cd_by_x))*per
                       })
                     })  
    
    names(cuminc) <- followup
    return(cuminc)
}
 
############################################################
## tally_cuminc
############################################################

tally_cuminc = function# Return cumulative incidence

##description<< For specified follow-up times, return cumulative
##incidence of an event

(followup,
    ### Vector of follow-up times
 etimes,
    ### Matrix of times to event
 event,
    ### Matrix of event indicator (1==yes, 0=no)
 etimes2=NULL,
    ### Optional matrix of times to event. Will be considered
    ### a screening arm compared to etimes as the control arm
 event2=NULL,
    ### Optional matrix of event indicator for 
    ### etimes2
 numdec=4 
    ### Decimals
) {

    cuminc <- lapply(followup,
               function(x) {
                   cd_by_x <- etimes<x & event==1
                   summarize_over_sims(cd_by_x, 'sum')
               })

    if (is.null(etimes2)) {
        return(data.frame(`Follow-Up Year`=followup,
                          do.call('rbind', cuminc)))

    } else {

        cuminc2 <- lapply(followup,
                   function(x) {
                       cd_by_x <- etimes<x & event==1
                       cd_by_x2 <- etimes2<x & event2==1

                       # Compute MRR, ARR and NNS for each sim
                       mrr <- matrix(colSums(cd_by_x2)/colSums(cd_by_x),
                                      nrow=1)
                       arr <- matrix(colMeans(cd_by_x)-colMeans(cd_by_x2),
                                      nrow=1)
                       nns <- 1/arr
                       if (sum(is.nan(mrr))>(ncol(etimes)/2))
                           warning('Warning: more than half of sims have no deaths in the control group')

                       # Now summarize over sims
                       toreturn <- data.frame(
                            `Follow-Up Year`=x,
                            Measure=c('Cumulative Incidence, Group 0',
                                      'Cumulative Incidence, Group 1',
                                      'MRR',
                                      'ARR',
                                      'NNS'),
                            rbind(
                               summarize_over_sims(cd_by_x, 'mean'),
                               summarize_over_sims(cd_by_x2, 'mean'),
                               summarize_over_sims(mrr, 'mean', na.rm=TRUE),
                               summarize_over_sims(arr, 'mean', na.rm=TRUE),
                               summarize_over_sims(nns[, is.finite(nns), 
                                                   drop=FALSE], 
                                                    'mean', na.rm=TRUE)
                            ))

                       # Replace summary means for MRR, ARR and NNS with
                       # means derived from the mean deaths in control
                       # and screening arms

                           # 1/23/15 edit: since we may not want to do this,
                           # I'm creating duplicate columns to save the original
                           # Estimates
                           toreturn = transform(toreturn, Estimate_Orig=Estimate)

                       screenrow <- which(toreturn$Measure==
                                          'Cumulative Incidence, Group 1')
                       controlrow <- which(toreturn$Measure==
                                           'Cumulative Incidence, Group 0')
                       toreturn[which(toreturn$Measure=='MRR'), 'Estimate'] <- 
                           round(
                               toreturn[screenrow, 'Estimate']/
                               toreturn[controlrow, 'Estimate'],
                           numdec)
                       toreturn[which(toreturn$Measure=='ARR'), 'Estimate'] <- 
                           round(
                               toreturn[controlrow, 'Estimate']-
                               toreturn[screenrow, 'Estimate'],
                           numdec)
                       toreturn[which(toreturn$Measure=='NNS'), 'Estimate'] <- 
                           round(
                                 1/toreturn[which(toreturn$Measure=='ARR'),
                                            'Estimate'],
                           numdec)
                    return(toreturn)
                   })
        return(do.call('rbind', cuminc2))
    }
    
}

############################################################
## age_specific_mortrate
############################################################

age_specific_eventrate = function# Return age-specific cancer mortality rates

##description<< Given age at entry and times to cancer and all-cause death,
##calculate age-specific cancer mortality rates

(ageentry,
    ### Matrix of ages at entry
 timealive,
    ### Matrix of time to all-death from age at entry
 ageEvent,
    ### Matrix of age at event
 equals=TRUE,
    ### Must ageEvent==ageentry+timealive to be 
    ### a qualifying event? This is TRUE when
    ### we are evaluating cancer-specific mortality
    ### and this properly identifies those who 
    ### died of cancer
 denominator=100000,
    ### Denominator for population rates
 prevalence=FALSE
    ### Return prevalence percent instead of rate
) {

    # Determine 5-year age-groups in the population
    agesatdeath <- ageentry + timealive
    age_groups <- sort(unique(floor(c(round(agesatdeath/5,0))*5)))
    age_groups <- c(age_groups, max(age_groups)+5)

    # For each age group, determine person-years lived in each sim
    rates <- sapply(as.character(age_groups),
                           function(a, agesatdeath, ageEvent,
                                    prev=prevalence, interval=5,
                                    n=nrow(ageEvent), p=ncol(ageEvent)) { 
                              a <- as.numeric(a)
                              upper <- a+interval

                              # Determine the subset of eligible individuals
                              # For cause-specific mort rates: 
                              # ageEvent==agesatdeath for cancer deaths
                              # For other rates:
                              # ageEvent must happen before death
                              if (equals) 
                                  subset <- (ageEvent==agesatdeath)
                              else 
                                  subset <- (ageEvent<=agesatdeath)
                              # Start out with inInterval=0
                              inInterval <- matrix(0, nrow=n, ncol=p)
                              # Evaluate interval events
                              happened <- ageEvent<upper 
                              inInterval[subset & ageEvent>=a & happened] <- 1
                              # check: 
                              # summary(c(agesatdeath[as.logical(inInterval)]))

                              # Now determine years alive
                              yrsalive <- agesatdeath-a
                              yrsalive[yrsalive<0] <- 0
                              if (prev) alive <- (yrsalive>0)
                              yrsalive[yrsalive>interval] <- interval

                              # Compute rates 
                              if (!prev) 
                                  rates <- colSums(inInterval)/colSums(yrsalive) 
                              else 
                                  rates <- colSums(subset & happened & alive)/colSums(alive)
                             
                              return(rates)
                                        
                           }, agesatdeath, ageEvent, prevalence,
                           USE.NAMES=TRUE)

    if (prevalence) denominator <- 100
    rates <- data.frame(`Age Group`=age_groups,
                        Rate=colMeans(rates)*denominator)

    return(rates)
}

############################################################
## mean_matrix_subgroup
############################################################

mean_matrix_subgroup = function# Return the mean of a numeric matrix

##description<< Optional subgroup indicator

(mmatrix,
    ### Matrix, with columns as different sims
 subset=NULL
    ### Optional logical matrix of same dimensions indicating
    ### a subset over which to compute the mean
) {

    if (is.null(subset)) subset <- matrix(TRUE, nrow=nrow(mmatrix),
                                          ncol=ncol(mmatrix))
    return(mean(mmatrix[subset]))
}

