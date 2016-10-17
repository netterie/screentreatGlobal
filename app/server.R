# Setup: libraries and functions
# A note about deploying: 
# rsconnect::deployApp('/Users/jeanette/Documents/jbirnbau/screentreatGlobal/app/', account='screeningandtreatment', appName='main')
library(shiny)
source('code.R')

shinyServer(function(input, output, session) {
  
  #TO MAKE DYNAMIC
  output$debug <- renderPrint({
    as.character(treattumor_props2(input$prop_ERpos,
                                  input$prop_Nodepos,
                                  input$tam.elig.control,
                                  input$tam.prop.control,
                                  input$chemo.elig.control,
                                  input$chemo.prop.control))
  })
  output$a0t <- renderUI({
    treatvec <- treattumor_props(as.numeric(input$prop_ERpos),
                                             as.numeric(input$prop_Nodepos),
                                             input$tam.elig.control,
                                             as.numeric(input$tam.prop.control),
                                             input$chemo.elig.control,
                                             as.numeric(input$chemo.prop.control))
    #textInput('prop.a0.t', 'Advanced cases', "0.06,0.14,0.24,0.56")
      textInput('prop.a0.t', 'Advanced cases', 
                paste(as.character(treatvec),collapse=','))
  })
  output$e0t <- renderUI({
    treatvec <- treattumor_props(as.numeric(input$prop_ERpos),
                                 as.numeric(input$prop_Nodepos),
                                 input$tam.elig.control,
                                 as.numeric(input$tam.prop.control),
                                 input$chemo.elig.control,
                                 as.numeric(input$chemo.prop.control))
    #textInput('prop.e0.t', 'Early cases', "0.06,0.14,0.24,0.56")
    textInput('prop.e0.t', 'Advanced cases', 
              paste(as.character(treatvec),collapse=','))
  })
  output$a1t <- renderUI({
    treatvec <- treattumor_props(as.numeric(input$prop_ERpos),
                                 as.numeric(input$prop_Nodepos),
                                 input$tam.elig.interv,
                                 as.numeric(input$tam.prop.interv),
                                 input$chemo.elig.interv,
                                 as.numeric(input$chemo.prop.interv))
    #textInput('prop.a1.t', 'Advanced cases', "0.3,0,0,0.7")
    textInput('prop.a1.t', 'Advanced cases', 
              paste(as.character(treatvec),collapse=','))
  })
  output$e1t <- renderUI({
    treatvec <- treattumor_props(as.numeric(input$prop_ERpos),
                                 as.numeric(input$prop_Nodepos),
                                 input$tam.elig.interv,
                                 as.numeric(input$tam.prop.interv),
                                 input$chemo.elig.interv,
                                 as.numeric(input$chemo.prop.interv))
    #textInput('prop.e1.t', 'Early cases', "0.3,0,0,0.7")
    textInput('prop.e1.t', 'Advanced cases', 
              paste(as.character(treatvec),collapse=','))
    
  })
  
  output$resultsTable <- renderPrint({
    rbind(
      mrr.annualInc(N=100000, 
                    p.inc=as.numeric(input$incidence)/100000, 
                    p.a=as.numeric(input$prop_a), 
                    p.s=as.numeric(input$prop_s), 
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
                    p.a=as.numeric(input$prop_a), 
                    p.s=as.numeric(input$prop_s), 
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
  })
  
  output$resultsGraph <- renderPlot({
    s <- 
      rbind(
        mrr.annualInc(N=100000, 
                      p.inc=as.numeric(input$incidence)/100000, 
                      p.a=as.numeric(input$prop_a), 
                      p.s=as.numeric(input$prop_s), 
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
                      p.a=as.numeric(input$prop_a), 
                      p.s=as.numeric(input$prop_s), 
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
                 Statistic=='Surv' & variable!='Intervention')
    sl <- transform(sl, 
                    Percent=round(100*value))
    
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

