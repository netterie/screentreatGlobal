# Setup: libraries and functions
# A note about deploying: 
# rsconnect::deployApp('/Users/jeanette/Documents/jbirnbau/screentreatGlobal/app/', account='screeningandtreatment', appName='main')
library(shiny)
library(ggplot2)
library(reshape2)
library(plyr)
source('code.R')

shinyServer(function(input, output, session) {
  
  # DEBUGGER
  output$debug <- renderPrint({
    as.character(print_vec(a0t.reactive()))
  })

  # TREATMENT-TUMOR SUBGROUP PROPORTIONS
  a0t.reactive <- reactive({
    treatvec <- treattumor_props(as.numeric(input$prop_ERpos),
                                             input$tam.elig.control,
                                             as.numeric(input$tam.prop.control),
                                             input$chemo.elig.control,
                                             as.numeric(input$chemo.prop.control))
    return(treatvec)
  })
  output$a0t <- renderUI({
    textInput('prop.a0.t', 'Advanced cases, control', 
                paste(as.character(a0t.reactive()),collapse=','))
  })
  e0t.reactive  <- reactive({
    treatvec <- treattumor_props(as.numeric(input$prop_ERpos),
                                 input$tam.elig.control,
                                 as.numeric(input$tam.prop.control),
                                 input$chemo.elig.control,
                                 as.numeric(input$chemo.prop.control))
    return(treatvec)
  })
  output$e0t <- renderUI({
    textInput('prop.e0.t', 'Early cases, control', 
              paste(as.character(e0t.reactive()),collapse=','))
  })
  a1t.reactive <- reactive({
    treatvec <- treattumor_props(as.numeric(input$prop_ERpos),
                                 input$tam.elig.interv,
                                 as.numeric(input$tam.prop.interv),
                                 input$chemo.elig.interv,
                                 as.numeric(input$chemo.prop.interv))
    return(treatvec)
  })
  output$a1t <- renderUI({
    textInput('prop.a1.t', 'Advanced cases, intervention', 
              paste(as.character(a1t.reactive()),collapse=','))
  })
  e1t.reactive <- reactive({
    treatvec <- treattumor_props(as.numeric(input$prop_ERpos),
                                 input$tam.elig.interv,
                                 as.numeric(input$tam.prop.interv),
                                 input$chemo.elig.interv,
                                 as.numeric(input$chemo.prop.interv))
    return(treatvec)
  })
  output$e1t <- renderUI({
    textInput('prop.e1.t', 'Early cases, intervention', 
              paste(as.character(e1t.reactive()),collapse=','))
    
  })
  prop_s <- reactive({ 1-(input$prop_a1/input$prop_a0) })
  
  # PARAMETER SUMMARY TABLES
  output$paramsum1 <- renderTable({
      data.frame(Parameter=c('Annual incidence per 100,000',
                             'Proportion ER+',
                             'Proportion surviving k years, advanced stage',
                             'Proportion surviving k years, early stage',
                             'Proportion presenting in advanced stage'),
                 Control=c(input$incidence,
                           input$prop_ERpos,
                           input$surv.adv,
                           input$surv.early,
                           input$prop_a0),
                 Intervention=c(NA,
                                NA,
                                NA,
                                NA,
                                input$prop_a1))

  }, NA.string='-')
  output$paramsum2 <- renderTable({
      data.frame(`ER Status`=c('ER+', '', '', '', '',
                               'ER-', '', '', '', ''),
                 Treatment=
                     c(rep(c('', 'None', 'Endocrine', 'Chemo', 'Endocrine+Chemo'),2))
                 ,
                 Control=c(NA,
                           a0t.reactive()[c('ERpos.None', 'ERpos.Tam', 'ERpos.Chemo', 
                                            'ERpos.TamChemo')],
                           NA,
                           a0t.reactive()[c('ERneg.None', 'ERneg.Tam', 'ERneg.Chemo')],
                           0),
                 Intervention=c(NA,
                           a1t.reactive()[c('ERpos.None', 'ERpos.Tam', 'ERpos.Chemo', 
                                            'ERpos.TamChemo')],
                           NA,
                           a1t.reactive()[c('ERneg.None', 'ERneg.Tam', 'ERneg.Chemo')],
                           0),
                 check.names=FALSE)

  }, NA.string='')
  output$paramsum3 <- renderTable({
      data.frame(`ER Status`=c('ER+', '', '', '', '',
                               'ER-', '', '', '', ''),
                 Treatment=
                     c(rep(c('', 'None', 'Endocrine', 'Chemo', 'Endocrine+Chemo'),2))
                 ,
                 Control=c(NA,
                           e0t.reactive()[c('ERpos.None', 'ERpos.Tam', 'ERpos.Chemo', 
                                            'ERpos.TamChemo')],
                           NA,
                           e0t.reactive()[c('ERneg.None', 'ERneg.Tam', 'ERneg.Chemo')],
                           0),
                 Intervention=c(NA,
                           e1t.reactive()[c('ERpos.None', 'ERpos.Tam', 'ERpos.Chemo', 
                                            'ERpos.TamChemo')],
                           NA,
                           e1t.reactive()[c('ERneg.None', 'ERneg.Tam', 'ERneg.Chemo')],
                           0),
                 check.names=FALSE)

  }, NA.string='')
  # Later, use this thread to improve formatting in the table
  # https://groups.google.com/forum/#!topic/shiny-discuss/2jlYOYFp2-A
  output$hazards <- renderTable({
      data.frame(`ER Status`=c('ER+', '', '', '',
                               'ER-', ''),
                 Treatment=
                     c('', 'Endocrine', 'Chemo', 'Endocrine+Chemo', '', 'Chemo')
                 ,
                 `Hazard Ratio`=
                     c(NA, 0.7, 0.775, 0.5425, NA, 0.775),
                 `Implied percent improvement in survival`=
                     c(NA, 30, 22.5, 45.75, NA, 22.5),
                 check.names=FALSE)

  }, NA.string='')
  # RESULTS TABLES
  output$resultsTable1 <- renderTable({
      subset(
      mrr.annualInc(N=100000, 
                    p.inc=as.numeric(input$incidence)/100000, 
                    p.a=as.numeric(input$prop_a0), 
                    p.s=as.numeric(prop_s()), 
                    m.a=exp.rate(as.numeric(input$surv.adv),
                                 year=as.numeric(input$year.surv)), 
                    m.e=exp.rate(as.numeric(input$surv.early),
                                 year=as.numeric(input$year.surv)),
                    k=5, 
                    h=as.numeric(unlist(strsplit(input$treat.hr,","))), 
                    p.a0.t=as.numeric(unlist(strsplit(input$prop.a0.t,","))), 
                    p.e0.t=as.numeric(unlist(strsplit(input$prop.e0.t,","))), 
                    p.a1.t=as.numeric(unlist(strsplit(input$prop.a1.t,","))), 
                    p.e1.t=as.numeric(unlist(strsplit(input$prop.e1.t,",")))),
             select=-Year)
  }, digits=0)
  output$resultsTable2 <- renderTable({
      subset(
      mrr.annualInc(N=100000, 
                    p.inc=as.numeric(input$incidence)/100000, 
                    p.a=as.numeric(input$prop_a0), 
                    p.s=as.numeric(prop_s()), 
                    m.a=exp.rate(as.numeric(input$surv.adv),
                                 year=as.numeric(input$year.surv)), 
                    m.e=exp.rate(as.numeric(input$surv.early),
                                 year=as.numeric(input$year.surv)),
                    k=10, 
                    h=as.numeric(unlist(strsplit(input$treat.hr,","))), 
                    p.a0.t=as.numeric(unlist(strsplit(input$prop.a0.t,","))), 
                    p.e0.t=as.numeric(unlist(strsplit(input$prop.e0.t,","))), 
                    p.a1.t=as.numeric(unlist(strsplit(input$prop.a1.t,","))), 
                    p.e1.t=as.numeric(unlist(strsplit(input$prop.e1.t,",")))),
             select=-Year)
  }, digits=0)
  
  # RESULTS GRAPHS
  output$resultsGraph <- renderPlot({
    s <- 
      rbind(
        mrr.annualInc(N=100000, 
                      p.inc=as.numeric(input$incidence)/100000, 
                      p.a=as.numeric(input$prop_a0), 
                      p.s=as.numeric(prop_s()), 
                      m.a=exp.rate(as.numeric(input$surv.adv),
                                   year=as.numeric(input$year.surv)), 
                      m.e=exp.rate(as.numeric(input$surv.early),
                                   year=as.numeric(input$year.surv)),
                      k=5, 
                      h=as.numeric(unlist(strsplit(input$treat.hr,","))), 
                      p.a0.t=as.numeric(unlist(strsplit(input$prop.a0.t,","))), 
                      p.e0.t=as.numeric(unlist(strsplit(input$prop.e0.t,","))), 
                      p.a1.t=as.numeric(unlist(strsplit(input$prop.a1.t,","))), 
                      p.e1.t=as.numeric(unlist(strsplit(input$prop.e1.t,",")))
        ),
        mrr.annualInc(N=100000, 
                      p.inc=as.numeric(input$incidence)/100000, 
                      p.a=as.numeric(input$prop_a0), 
                      p.s=as.numeric(prop_s()), 
                      m.a=exp.rate(as.numeric(input$surv.adv),
                                   year=as.numeric(input$year.surv)), 
                      m.e=exp.rate(as.numeric(input$surv.early),
                                   year=as.numeric(input$year.surv)),
                      k=10, 
                      h=as.numeric(unlist(strsplit(input$treat.hr,","))), 
                      p.a0.t=as.numeric(unlist(strsplit(input$prop.a0.t,","))), 
                      p.e0.t=as.numeric(unlist(strsplit(input$prop.e0.t,","))), 
                      p.a1.t=as.numeric(unlist(strsplit(input$prop.a1.t,","))), 
                      p.e1.t=as.numeric(unlist(strsplit(input$prop.e1.t,","))))
      )
    s <- transform(s, `Gained by Intervention`=Intervention-Control,
                   check.names=FALSE)
    sl <- subset(melt(s,id.vars=c('Year', 'Statistic')),
                 Statistic=='Percent surviving' & variable!='Intervention')
    sl <- transform(sl, 
                    Percent=round(value))
    
    sl = ddply(sl, .(Year), transform, pos = (cumsum(Percent) - 0.5 * Percent))
    #sl$label = paste0(sprintf("%.0f", sl$Percent), "%")
    sl$label = as.character(sl$Percent)
    
    
    g <- ggplot(sl, aes(x = factor(Year), y = Percent, fill = variable)) +
      geom_bar(stat = "identity", width = .7) +
      geom_text(aes(y = pos, label = label), size = 4) +
      theme(text = element_text(size=10)) + 
      scale_x_discrete(name='Years after intervention') + 
      scale_y_continuous('Percent Surviving',limits=c(0,100)) + 
      theme_bw()
    g + theme(legend.position='top') + theme(legend.title=element_blank())
  })
})

