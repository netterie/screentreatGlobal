library(shiny)

shinyUI(fluidPage(
  
  titlePanel("Welcome to the cancer screening and treatment impact model"),
  
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
                Tamoxifen (20%). However, the Tamoxifen is not targeted to the 30% of 
                women who are ER+. This means that some ER- receive Tamoxifen and but do
               not benefit from it.'),
             h5('Standard of Care: "Intervention" Scenario'),
             p('The intervention modeled is ER screening, which leads to Tamoxifen being
               administered only to the 30% of women who are ER+. All of these women
               receive a survival benefit from the Tamoxifen.')
    ),
    "Natural History Parameters",
    tabPanel("Cancer Incidence",
             h4('Enter the expected annual incidence per 100,000 women and the proportion 
                  of cancers that are advanced-stage (stage III/IV), ER positive, and 
                  node positive at the time of clinical diagnosis.'),
             p('Annual incidence in Tanzania is approximated as the same as Uganda, 
                which is 40 cases per year in an age-standardized population of 100,000. 
                About 85% of cases in Tanzania are advanced stage at the time of clinical 
                diagnosis, i.e. diagnosis in the absence of screening. About 30% are ER
                positive.'),
             # Incidence
             numericInput('incidence', label='Annual incidence per 100,000', 
                          40, min = 0, max = 100000, step = NA, 
                          width = NULL),
             # Proportion advanced
             sliderInput("prop_a", label = "Proportion advanced",
                         min=0, max=1, step=0.01, value=0.85),
             # Proportion ER+
             sliderInput("prop_a", label = "Proportion ER positive",
                         min=0, max=1, step=0.01, value=0.30),
             # Proportion node+
             sliderInput("prop_a", label = "Proportion node positive",
                         min=0, max=1, step=0.01, value=0),
             #Numeric Inputs
             numericInput("min_val", "Enter Minimum Value", 1993),
             numericInput("max_val", "Enter Maximum Value", 2013)
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
    "Treatment Parameters",
    tabPanel("Control Scenario",
              h4('Specify who is eligible for each treatment, 
                and what proportion of eligible cases receive it.'),
              p('The current standard of care in Tanzania is that about 20% of all women 
                  receive Tamoxifen, even though only the cases who are 
                  estrogen-receptor positive (ER+) can actually benefit from it.'),
              br(),
              h5('TAMOXIFEN'),
              radioButtons("tam.elig.control", "Who is eligible for Tamoxifen?",
                           c("All" = 'All',
                             "ER+ only" = 'ERpos')),
              sliderInput('tam.prop.control', label='What proportion of eligible women receive Tamoxifen?', 
                          0.20, min = 0, max = 1, step = .01),
              br(),
              h5('CHEMOTHERAPY'),
              radioButtons("chemo.elig.control", "Who is eligible for chemotherapy?",
                          c("All" = 'All',
                            "ER- only" = 'ERneg',
                            "ER- and ER+Node+" = 'ERnegERposNodepos'
                            )),
              sliderInput('chemo.prop.control', label='What proportion of eligible women receive chemotherapy?', 
                          0.0, min = 0, max = 1, step = .01)
              ),
    tabPanel("Intervention Scenario",
             h4('Specify who is eligible for each treatment, 
                and what proportion of eligible cases receive it.'),
             p('In the default intervention, all ER+ women receive Tamoxifen.'),
             br(),
             h5('TAMOXIFEN'),
             radioButtons("tam.elig.interv", "Who is eligible for Tamoxifen?",
                          c("All" = 'All',
                            "ER+ only" = 'ERpos'),
                          selected='ERpos'),
             sliderInput('tam.prop.interv', label='What proportion of eligible women receive Tamoxifen?', 
                         1.0, min = 0, max = 1, step = .01),
             br(),
             h5('CHEMOTHERAPY'),
             radioButtons("chemo.elig.interv", "Who is eligible for chemotherapy?",
                          c("All" = 'All',
                            "ER- only" = 'ERneg',
                            "ER- and ER+Node+" = 'ERnegERposNodepos'
                          ),
                          selected='All'),
             sliderInput('chemo.prop.interv', label='What proportion of eligible women receive chemotherapy?', 
                         0.0, min = 0, max = 1, step = .01)
    ),
    "Screening Parameters",
    tabPanel("Intervention Stage Shift",
              h4('Enter the proportion of advanced-stage cases that are stage-shifted 
                  by the intervention. If the intervention only changes treatment, 
                  not stage at diagnosis, you may enter a stage shift of 0.'),
              p('There is no stage-shift associated with the default intervention.'),
              sliderInput("prop_s", label = "Proportion stage-shifted",
                         min=0, max=1, step=0.05, value=0.0)
             ),
    "Optional Parameters",
    tabPanel("Tumor Type and Treatment Details",
              h4('OPTIONAL: Define the relevant treatments and tumor subgroups, 
                and the hazard ratios of treatment for each group.'),
             p('The inputs on this page reflect the tumor characteristics and
                the treatment proportions specified on the previous pages. However,
                the resulting treatment-tumor subgroups and proportions may be edited
                here directly if finer control of the model is desired.'),
             
             #display dynamic UI
             uiOutput("slider"),
             br(),
             h5('TUMOR-TREATMENT GROUPS'),
             textInput('treat.vec', 'Treatment names (separate by commas)', 
                       "Tam.ERp,Tam.ERn,None.ERp,None.ERn"),
             br(),
             h5('EFFICACIES BY TUMOR-TREATMENT GROUP'),
             p('Hazard ratios (HRs) indicate the survival benefit of treatment 
                in that group, e.g., HR=1 means no benefit, whereas HR=0.7 indicates 
                a 30% survival benefit.'),
             textInput('treat.hr', 'Treatment hazard ratios (separated by commas)', 
                       "0.7,1,1,1"),
             br(),
             h5('CONTROL TREATMENT PROPORTIONS'),
             uiOutput('a0t'),
             uiOutput('e0t'),
             br(),
             h5('INTERVENTION TREATMENT PROPORTIONS'),
             uiOutput('a1t'),
             uiOutput('e1t')
             ),
    "Results",
    tabPanel("Summary Table",
              verbatimTextOutput('resultsTable')
             ),
    tabPanel("Survival Graph",
              plotOutput('resultsGraph')
             )
  )
))


