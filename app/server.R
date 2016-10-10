# Setup: libraries and functions
library(shiny)
source('code.R')

shinyServer(function(input, output, session) {
  
  #make dynamic slider
  output$slider <- renderUI({
    sliderInput("inSlider", "Slider", min=input$min_val, max=input$max_val, value=2000)
  })
  
  #TO MAKE DYNAMIC
  output$a0t <- renderUI({
      textInput('prop.a0.t', 'Advanced cases', "0.06,0.14,0.24,0.56")
  })
  output$e0t <- renderUI({
    textInput('prop.e0.t', 'Early cases', "0.06,0.14,0.24,0.56")
  })
  output$a1t <- renderUI({
    textInput('prop.a1.t', 'Advanced cases', "0.3,0,0,0.7")
  })
  output$e1t <- renderUI({
    textInput('prop.e1.t', 'Early cases', "0.3,0,0,0.7")
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
