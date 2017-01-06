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

  # CONVERSION FROM PERCENTS TO PROPORTIONS
  prop_ERpos <- reactive({ input$prop_ERpos/100 })
  surv.adv <- reactive({ input$surv.adv/100 })
  surv.early <- reactive({ input$surv.early/100 })
  prop_a0 <- reactive({ input$prop_a0/100 })
  prop_a1 <- reactive({ input$prop_a1/100 })
  tam.prop.control <- reactive({ input$tam.prop.control/100 })
  chemo.prop.control <- reactive({ input$chemo.prop.control/100 })
  tam.prop.interv <- reactive({ input$tam.prop.interv/100 })
  chemo.prop.interv <- reactive({ input$chemo.prop.interv/100 })

  # TREATMENT-TUMOR SUBGROUP PROPORTIONS
  a0t.reactive <- reactive({
    treatvec <- treattumor_props(as.numeric(prop_ERpos()),
                                             input$tam.elig.control,
                                             as.numeric(tam.prop.control()),
                                             input$chemo.elig.control,
                                             as.numeric(chemo.prop.control()))
    return(treatvec)
  })
  output$a0t <- renderUI({
    textInput('prop.a0.t', 'Advanced cases, control', 
                paste(as.character(a0t.reactive()),collapse=','))
  })
  e0t.reactive  <- reactive({
    treatvec <- treattumor_props(as.numeric(prop_ERpos()),
                                 input$tam.elig.control,
                                 as.numeric(tam.prop.control()),
                                 input$chemo.elig.control,
                                 as.numeric(chemo.prop.control()))
    return(treatvec)
  })
  output$e0t <- renderUI({
    textInput('prop.e0.t', 'Early cases, control', 
              paste(as.character(e0t.reactive()),collapse=','))
  })
  a1t.reactive <- reactive({
    treatvec <- treattumor_props(as.numeric(prop_ERpos()),
                                 input$tam.elig.interv,
                                 as.numeric(tam.prop.interv()),
                                 input$chemo.elig.interv,
                                 as.numeric(chemo.prop.interv()))
    return(treatvec)
  })
  output$a1t <- renderUI({
    textInput('prop.a1.t', 'Advanced cases, intervention', 
              paste(as.character(a1t.reactive()),collapse=','))
  })
  e1t.reactive <- reactive({
    treatvec <- treattumor_props(as.numeric(prop_ERpos()),
                                 input$tam.elig.interv,
                                 as.numeric(tam.prop.interv()),
                                 input$chemo.elig.interv,
                                 as.numeric(chemo.prop.interv()))
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
                             'Percent ER+',
                             'Percent surviving k years, advanced stage',
                             'Percent surviving k years, early stage',
                             'Percent presenting in advanced stage'),
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
                           100*a0t.reactive()[c('ERpos.None', 'ERpos.Tam', 'ERpos.Chemo', 
                                            'ERpos.TamChemo')],
                           NA,
                           100*a0t.reactive()[c('ERneg.None', 'ERneg.Tam', 'ERneg.Chemo')],
                           0),
                 Intervention=c(NA,
                           100*a1t.reactive()[c('ERpos.None', 'ERpos.Tam', 'ERpos.Chemo', 
                                            'ERpos.TamChemo')],
                           NA,
                           100*a1t.reactive()[c('ERneg.None', 'ERneg.Tam', 'ERneg.Chemo')],
                           0),
                 check.names=FALSE)

  }, NA.string='', digits=0)
  output$paramsum3 <- renderTable({
      data.frame(`ER Status`=c('ER+', '', '', '', '',
                               'ER-', '', '', '', ''),
                 Treatment=
                     c(rep(c('', 'None', 'Endocrine', 'Chemo', 'Endocrine+Chemo'),2))
                 ,
                 Control=c(NA,
                           100*e0t.reactive()[c('ERpos.None', 'ERpos.Tam', 'ERpos.Chemo', 
                                            'ERpos.TamChemo')],
                           NA,
                           100*e0t.reactive()[c('ERneg.None', 'ERneg.Tam', 'ERneg.Chemo')],
                           0),
                 Intervention=c(NA,
                           100*e1t.reactive()[c('ERpos.None', 'ERpos.Tam', 'ERpos.Chemo', 
                                            'ERpos.TamChemo')],
                           NA,
                           100*e1t.reactive()[c('ERneg.None', 'ERneg.Tam', 'ERneg.Chemo')],
                           0),
                 check.names=FALSE)

  }, NA.string='', digits=0)
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
                    p.a=as.numeric(prop_a0()), 
                    p.s=as.numeric(prop_s()), 
                    m.a=exp.rate(as.numeric(surv.adv()),
                                 year=as.numeric(input$year.surv)), 
                    m.e=exp.rate(as.numeric(surv.early()),
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
                    p.a=as.numeric(prop_a0()), 
                    p.s=as.numeric(prop_s()), 
                    m.a=exp.rate(as.numeric(surv.adv()),
                                 year=as.numeric(input$year.surv)), 
                    m.e=exp.rate(as.numeric(surv.early()),
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
                      p.a=as.numeric(prop_a0()), 
                      p.s=as.numeric(prop_s()), 
                      m.a=exp.rate(as.numeric(surv.adv()),
                                   year=as.numeric(input$year.surv)), 
                      m.e=exp.rate(as.numeric(surv.early()),
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
                      p.a=as.numeric(prop_a0()), 
                      p.s=as.numeric(prop_s()), 
                      m.a=exp.rate(as.numeric(surv.adv()),
                                   year=as.numeric(input$year.surv)), 
                      m.e=exp.rate(as.numeric(surv.early()),
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

