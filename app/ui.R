library(shiny)

shinyUI(fluidPage(
  
  titlePanel("Welcome to the cancer early detection and treatment impact model"),
  
  navlistPanel(
    tabPanel("Introduction",
             h4('Overview'),
             p('This interface allows you to model the survival benefit of a 
                screening and/or treatment intervention in a population of your choosing. 
                The model posits a simple stage-shift mechanism of screening benefit. 
                Specify inputs using the navigation sidebar.'),
             h4('Default Model: Breast Cancer in Tanzania'),
             p('The input defaults reflect an breast cancer example in Tanzania:'),
             h5('Standard of Care: "Control" Scenario'),
             p('Because chemotherapy treatment is logistically and financially difficult
                in Tanzania, women typically either receive no adjuvant treatment (80%) or
                endocrine therapy (20%). However, the endocrine therapy is not targeted to the 50% of 
                women who are ER+. This means that some ER- receive endocrine therapy and but do
               not benefit from it.'),
             h5('Standard of Care: "Intervention" Scenario'),
             p('The intervention modeled is ER screening, which leads to endocrine therapy being
               administered only to the 50% of women who are ER+. All of these women
               receive a survival benefit from the endocrine therapy.')
    ),
    "Natural History",
    tabPanel("Cancer Incidence",
             h4('Enter the expected annual incidence per 100,000 women and the proportion 
                  of cancers that are ER positive.'),
             p('Annual incidence in Tanzania is patterned after Uganda, 
                where ther are about 60 cases per year among 100,000 women ages 30-49. 
                About 50% are ER positive.'),
             # Incidence
             numericInput('incidence', label='Annual incidence per 100,000', 
                          60, min = 0, max = 100000, step = NA, 
                          width = NULL),
             # Proportion ER+
             sliderInput("prop_ERpos", label = "Proportion ER positive",
                         min=0, max=1, step=0.01, value=0.50)
             ),
    tabPanel("Stage-Specific Survival",
              h4('Select a year, k, by which you will specify the proportion of cases 
                surviving at k years after diagnosis, or "baseline survival."'),
              p('Baseline survival should typically be survival in the absence of 
                systemic treatment, but it can be treated survival if the 
                intervention does not impact treatment. You must specify k-year 
                survival for advanced- and early-stage cases separately.'),
              p('Data from Uganda suggest that advanced-stage cases in Tanzania may 
                have a 5-year survival rate of 35%. Given the low incidence of 
                early-stage cancer in Tanzania and similar countries, survival data 
                are sparse for early stage caes. We approximate early-stage survival 
                using historical data from the US in the 1940s, which puts 5-year survival 
                of localized cases at 80%.'),
             selectInput('year.surv', label='Year of survival statistic, k', 
                         choices=c(5,10), selected=5),
             
             sliderInput('surv.adv', label='Advanced cases: baseline survival at k years', 
                         0.35, min = 0, max = 1, step = .01),
             
             sliderInput('surv.early', label='Early cases: baseline survival at k years', 
                         0.80, min = 0, max = 1, step = .01)
             ),
    "Impact of Early Detection",
    tabPanel("Stage Distributions",
              h4('Enter the proportion of advanced-stage cases in the control and intervention scenarios.'),
              p('In the control scenario, this is the typical proportion of cases who are 
                advanced-stage at the time of clinical diagnosis. In the 
                intervention scenario, the proportion presenting in advanced stage may
                decrease due to early detection efforts.'),
              p('The default intervention involves no formal early detection strategy,
                only changes to the treatments available (see next two panels).'),
             # Proportion advanced
              h5('CONTROL SCENARIO'),
             sliderInput("prop_a0", label = "Proportion advanced",
                         min=0, max=1, step=0.01, value=0.85),
             # Proportion advanced
              h5('INTERVENTION SCENARIO'),
             sliderInput("prop_a1", label = "Proportion advanced",
                         min=0, max=1, step=0.01, value=0.85)
             ),
    "Treatments Available",
    tabPanel("Control Scenario",
              h4('Specify who is eligible for each treatment, 
                and what proportion of eligible cases receive it.'),
              p('The current standard of care in Tanzania is that about 20% of all women 
                  receive endocrine therapy, even though only the cases who are 
                  estrogen-receptor positive (ER+) can actually benefit from it.'),
              br(),
              h5('ENDOCRINE THERAPY'),
              radioButtons("tam.elig.control", "Who is eligible for endocrine therapy?",
                           c("All" = 'All',
                             "ER+ only" = 'ERpos')),
              sliderInput('tam.prop.control', label='What proportion of eligible women receive endocrine therapy?', 
                          0.20, min = 0, max = 1, step = .01),
              br(),
              h5('CHEMOTHERAPY'),
              radioButtons("chemo.elig.control", "Who is eligible for chemotherapy?",
                          c("All" = 'All',
                            "ER- only" = 'ERneg',
                            "ER- and advanced-stage ER+" = 'ERnegERposAdv'
                            )),
              sliderInput('chemo.prop.control', label='What proportion of eligible women receive chemotherapy?', 
                          0.0, min = 0, max = 1, step = .01)
              ),
    tabPanel("Intervention Scenario",
             h4('Specify who is eligible for each treatment, 
                and what proportion of eligible cases receive it.'),
             p('In the default intervention, all ER+ women receive endocrine therapy.'),
             br(),
             h5('ENDOCRINE THERAPY'),
             radioButtons("tam.elig.interv", "Who is eligible for endocrine therapy?",
                          c("All" = 'All',
                            "ER+ only" = 'ERpos'),
                          selected='ERpos'),
             sliderInput('tam.prop.interv', label='What proportion of eligible women receive endocrine therapy?', 
                         1.0, min = 0, max = 1, step = .01),
             br(),
             h5('CHEMOTHERAPY'),
             radioButtons("chemo.elig.interv", "Who is eligible for chemotherapy?",
                          c("All" = 'All',
                            "ER- only" = 'ERneg',
                            "ER- and advanced-stage ER+" = 'ERnegERposAdv'
                          ),
                          selected='All'),
             sliderInput('chemo.prop.interv', label='What proportion of eligible women receive chemotherapy?', 
                         0.0, min = 0, max = 1, step = .01)
    ),
    "Parameter Summary",
    tabPanel("User-Defined Parameters",
              h4('Review the selected parameter values'),
             p('The parameters specified on the previous pages are summarized below.
               To make changes, revisit the previous pages.'),
             br(),
             h5('INCIDENCE AND TUMOR CHARACTERISTICS'),
             tableOutput('paramsum1'),
             br(),
             h5('ADVANCED-STAGE TREATMENTS'),
             em('Values represent percents falling into each group. Columns
                should sum to 100%'),
             tableOutput('paramsum2'),
             br(),
             h5('EARLY-STAGE TREATMENTS'),
             em('Values represent percents falling into each group. Columns
                should sum to 100%'),
             tableOutput('paramsum3'),

             br(), br(), br(), br(), br(), br(),
             h4('ADVANCED CONTROLS - DO NOT EDIT'),
             ## Can bring back the displays below to regain control beyond
             ## what is allowed in the user options on previous tabs
             br(),
#             h5('TUMOR-TREATMENT GROUPS'),
             textInput('treat.vec', 'Treatment names (separate by commas)', 
                        'ERneg.Tam, ERneg.Chemo, ERneg.None,
                        ERpos.TamChemo, ERpos.Tam, ERpos.Chemo, ERpos.None'
             ),
             br(),
#             h5('EFFICACIES BY TUMOR-TREATMENT GROUP'),
             textInput('treat.hr', 'Treatment hazard ratios (separated by commas)', 
                       '1, 0.775, 1, 0.5425, 0.7, 0.775, 1'
                       ),
             br(),
#             h5('CONTROL TREATMENT PROPORTIONS'),
             uiOutput('a0t'),
             uiOutput('e0t'),
             br(),
#             h5('INTERVENTION TREATMENT PROPORTIONS'),
             uiOutput('a1t'),
             uiOutput('e1t')
             ),
    tabPanel("Fixed Parameters",
             h4('The benefits of treatments are sourced from the literature'),
             p('Meta-analyses from the Early Breast Cancer Trialists Collaborative Group
               have summarized the benefits of endocrine therapy and chemotherapy across
               regimens. The benefits are quantified as "hazard ratios" on survival, which
               indicate the percent improvement in survival due to treatment'),
             br(),
             h5('SURVIVAL BENEFIT OF TREATMENT'),
             em('All treatments not specified in the table have a hazard ratio of 1,
                i.e. there is no impact on survival. So for example, endocrine therapy
                for ER- tumors has a hazard ratio of 1 and thus is not in the table.'),
             tableOutput('hazards'),
             p('References:'),
             em('1.
                Early Breast Cancer Trialists’ Collaborative Group (EBCTCG), Davies C, Godwin J, Gray R, Clarke M, Cutter D, et al. Relevance of breast cancer hormone receptors and other factors to the efficacy of adjuvant tamoxifen: patient-level meta-analysis of randomised trials. Lancet. 2011 Aug 27;378(9793):771–84. 
                '),
                br(),
             em('2.
                Early Breast Cancer Trialists’ Collaborative Group (EBCTCG), Peto R, Davies C, Godwin J, Gray R, Pan HC, et al. Comparisons between different polychemotherapy regimens for early breast cancer: meta-analyses of long-term outcome among 100,000 women in 123 randomised trials. Lancet. 2012 Feb 4;379(9814):432–44. 
                ')
             ),
    "Results",
    tabPanel("Results Tables",
#              verbatimTextOutput('debug'),
             h4('Results after 5 years'),
             tableOutput('resultsTable1'),
             br(),
             h4('Results after 10 years'),
             tableOutput('resultsTable2')
             ),
    tabPanel("Results Plots",
             h4('Percent Surviving'),
              plotOutput('resultsGraph')
             )
  )
))


