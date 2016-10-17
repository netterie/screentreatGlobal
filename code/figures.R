
############################################################
# Compile results into one Adjuvant-O style graph
############################################################


##############################
# Setup
##############################

base.path <- '/Users/jeanette/Documents/jbirnbau/screentreatGlobal/examples'

# Models and their labels
ids <- data.frame(models=c('tza_3_chemo1b',
                           'tza_3_chemo1b',
                           'tza_3_chemo1c'
                           ),
                  subset=c('Control',
                           'Screen',
                           'Screen'),
                  labels=c('ER Screening Only',
                           'ER Screening plus downstaging\nto 60% advanced',
                           'ER Screening plus downstaging\nto 35% advanced'
                           ))



# Compile results
results <- list()
for (a in c('agesAll', 'ages30-49', 'ages50-69')) {

    if (a=='agesAll') pasteA = '' else pasteA <- a

    for (i in 1:nrow(ids)) {
        # Read in results
        results[[paste0(a,i)]] <- transform(subset(read.csv(file.path(base.path, 
                                           paste0(ids$models[i],pasteA),
                                           'output',
                                           'simpleSurvival.csv'),
                                           header=TRUE),
                               variable==ids$subset[i]),
                                  Scenario=ids$labels[i],
                                  Ages=gsub('ages','',a))
    }
}

resultsA <- do.call('rbind', results)

##############################
# Compute incremental results
##############################

# First, repeat the Standard of Care results for the downstaging scenarios
standard <- subset(resultsA, variable=='Control' & Trial=='withoutER')
standard1 <- standard
standard$Scenario='ER Screening plus downstaging\nto 60% advanced'
standard1$Scenario='ER Screening plus downstaging\nto 35% advanced' 

resultsA <- rbind(resultsA, standard, standard1)

# Compute the incremental results
resultsA <- within(resultsA, {
                       Incremental=''
                       Incremental[variable=='Control'&
                                 Trial=='withoutER'] <- 'Base'
                       Incremental[Trial=='withER'] <- 'PlusTam'
                       Incremental[Trial=='withERchemo'] <- 
                           'PlusChemo'
                                  })

## We're not interested in downstaging without targeted treatment
resultsA <- subset(resultsA, !is.na(Incremental))

resultsC <- cast(resultsA, Ages+Time+Scenario~Incremental)
resultsC <- transform(resultsC, 
                      IncrementalChemo=PlusChemo-PlusTam,
                      IncrementalTam=PlusTam-Base)

resultsM <- subset(melt(resultsC, id.vars=c('Ages', 'Time', 'Scenario')),
                   variable!='PlusChemo' & variable!='PlusTam')

# Edit labels
resultsM <- within(resultsM, {
                       order=NA
                       Treatment=''
                       Treatment[variable=='Base'] <- 'Standard of Care'
                       order[variable=='Base'] <- 1
                       Treatment[variable=='IncrementalTam'] <- 'Tamoxifen for ER+'
                       order[variable=='IncrementalTam'] <- 2
                       Treatment[variable=='IncrementalChemo'] <- 
                           'Tamoxifen for ER+\nChemo for ER-'
                       order[variable=='IncrementalChemo'] <- 3
                   })
resultsM <- within(resultsM, {
                       Treatment=NA
                       Treatment[variable=='Base'] <- 1
                       Treatment[variable=='IncrementalTam'] <- 2
                       Treatment[variable=='IncrementalChemo'] <- 3
                       Treatment <- factor(Treatment, levels=c(1,2,3),
                                           labels=c('Standard of Care',
                                                    'Tamoxifen for ER+', 
                                                    'Tamoxifen for ER+\nChemo for ER-'))
                                  })

resultsM <- resultsM[order(resultsM$Treatment),]

# Factor orders
scenario.labels=c('ER Screening Only',
           'ER Screening plus downstaging\nto 60% advanced',
           'ER Screening plus downstaging\nto 35% advanced')
resultsM <- transform(resultsM,
                      Scenario=factor(as.character(Scenario),
                                      levels=scenario.labels,
                                      labels=scenario.labels))


##############################
# Graphs
##############################

# Compute labels and label positions
resultsM = ddply(resultsM, .(Ages, Time,Scenario), transform, 
                 pos = (cumsum(value) - 0.5 * value))
#resultsM$label = paste0(sprintf("%.0f", resultsM$value), "%")
resultsM$label = as.character(resultsM$value)

# Plotting function
plotfun <- function(df) {
    g <- ggplot(df, aes(x = factor(Time), y = value, fill = Treatment)) +
       geom_bar(stat = "identity", width = .7) +
       geom_text(aes(y = pos, label = label), size = 2) +
       theme(text = element_text(size=10)) + 
       scale_x_discrete(name='Years after intervention') + 
       scale_y_continuous('Percent Surviving',limits=c(0,100)) + 
       theme_bw() + facet_grid(.~Scenario)
    g + theme(legend.position='bottom') + theme(legend.title=element_blank())
}

g <- ggplot(subset(resultsM,Ages=='All'), 
            aes(x = factor(Time), y = value, fill = Treatment)) +
   geom_bar(stat = "identity", width = .7) +
   geom_text(aes(y = pos, label = label), size = 2) +
   theme(text = element_text(size=10)) + 
   scale_x_discrete(name='Years after intervention') + 
   scale_y_continuous('Percent Surviving',limits=c(0,100)) + 
   theme_bw() + facet_grid(.~Scenario)
g + theme(legend.position='bottom') + theme(legend.title=element_blank()) +
    ggtitle('All Ages')

plotfun(subset(resultsM, Ages=='All')) + ggtitle('All Ages')
plotfun(subset(resultsM, Ages=='30-49')) + ggtitle('Ages 30-49')
plotfun(subset(resultsM, Ages=='50-69')) + ggtitle('Ages 50-69')



