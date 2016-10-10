
# Convert survivals to rates
exp.rate = function(cum.surv, year=10, haz=1) {log(cum.surv)/(-year*haz)}
# Incidence is 39.8 per 100,000
# 85% advanced stage
# Two treatments: none, or tamoxifen at HR=0.7, but it only 
#    applies to ER+ who comprise 30%. Under control scenario,
#    only 30%*0.2= get HR=0.7. But 
#    under intervention, 30% get HR=0.7 due to perfect targeting
#    of the ER+

# General mortality function
mort.k = function(inc, p.a, p.s, m, k, p.t, h,
                  advanced=TRUE, shift=FALSE) {
  # Number of treatments
  n = length(p.t)
  # Proportion dying in treatment group i, year k
  prop.mort.k = 1-exp(rep(-k*m,n)*h)
  # Weight by treatment group proportions
  mort.k = sum(p.t*prop.mort.k)
  # Scale by incidence, stage proportions and stage shift
  if (advanced) {
    s = p.a 
    c = -1
  } else {
    s = 1-p.a
    c = 1
  }
  if (shift) s = s + c*p.a*p.s
  return(inc*s*mort.k)
}

# Instead of cumulative mortality, the # dying IN year k
mort.pdf = function(k, inc, p.a, p.s, m, p.t, h,
                    advanced=TRUE, shift=FALSE) {
  # Number of treatments
  n = length(p.t)
  # Proportion dying in treatment group i, year k
  cum.mort.k = 1-exp(rep(-k*m,n)*h)
  # Proportion who died in k-1 or before
  cum.mort.kminus1 = 1-exp(rep(-(k-1)*m,n)*h)
  # Proportion IN year k
  prop.mort.k = cum.mort.k-cum.mort.kminus1
  # Weight by treatment group proportions
  mort.k = sum(p.t*prop.mort.k)
  # Scale by incidence, stage proportions and stage shift
  if (advanced) {
    s = p.a 
    c = -1
  } else {
    s = 1-p.a
    c = 1
  }
  if (shift) s = s + c*p.a*p.s
  return(inc*s*mort.k)
}

# Years of life lived, based on annual inc and k-vector
yll.annualInc = function(inc,mort.vec) {
  k=length(mort.vec)
  years.tab <- data.frame(cohort=k:1,
                          mort.vec,
                          full.years=1:k,
                          inc=inc,
                          not.k.years=sapply(1:k,
                                             FUN=function(x) {sum(mort.vec[1:x])}))
  years.tab <- transform(years.tab,
                         yll.full=full.years*(inc-not.k.years),
                         yll.notfull=cohort*mort.vec)
  return(sum(years.tab$yll.full+years.tab$yll.notfull))
}


# Trying a function that uses annualized incidence
mrr.annualInc = function(N, p.inc, p.a, p.s, m.a, m.e, k, h, 
                         p.a0.t, p.e0.t, p.a1.t, p.e1.t) {
  # Annual incidence
  inc = N*p.inc
  # Control
  # Multiply by k:1 because of how many cohorts have the opportunity
  # to die at 1,2,...k years. If k=10, only 1 cohort h
  mort.a0.k = sapply(1:k, mort.pdf, inc, p.a, p.s, m.a, p.a0.t, h, 
                     advanced=TRUE, shift=FALSE)
  mort.e0.k = sapply(1:k, mort.pdf, inc, p.a, p.s, m.e, p.e0.t, h, 
                     advanced=FALSE, shift=FALSE)
  # Intervention
  mort.a1.k = sapply(1:k, mort.pdf, inc, p.a, p.s, m.a, p.a1.t, h, 
                     advanced=TRUE, shift=TRUE)
  mort.e1.k = sapply(1:k, mort.pdf, inc, p.a, p.s, m.e, p.e1.t, h, 
                     advanced=FALSE, shift=TRUE)
  # Denominator and numerator of MRR 
  denominator = sum((k:1)*mort.a0.k) + sum((k:1)*mort.e0.k)
  numerator = sum((k:1)*mort.a1.k) + sum((k:1)*mort.e1.k)
  # Years of life
  control.yll <- yll.annualInc(inc*p.a,mort.a0.k) +
    yll.annualInc(inc*(1-p.a),mort.e0.k)
  intervention.yll <- yll.annualInc(inc*p.a,mort.a1.k) +
    yll.annualInc(inc*(1-p.a),mort.e1.k)
  
  
  # Results
  dfres <- data.frame(Year=k,
                      Statistic=c('Mort',
                                  'Women Alive',
                                  'Surv',
                                  'MRR',
                                  'ARR',
                                  'Years of Life Saved'),
                      Control=c(denominator,
                                inc*k-denominator,
                                (inc*k-denominator)/(inc*k),
                                1,
                                0,
                                0),
                      Intervention=c(numerator,
                                     inc*k-numerator,
                                     (inc*k-numerator)/(inc*k),
                                     round(numerator/denominator,3),
                                     round(denominator-numerator,3),
                                     round(intervention.yll-control.yll,3))
  )
  return(dfres)
}

# OLD
# I think this is a less precise version
# Input needs to be sort of half cumulative incidence
mrr = function(N, p.inc, p.a, p.s, m.a, m.e, k, h, 
               p.a0.t, p.e0.t, p.a1.t, p.e1.t) {
  # Incidence
  inc = N*p.inc
  # Control
  mort.a0.k = mort.k(inc, p.a, p.s, m.a, k, p.a0.t, h,
                     advanced=TRUE, shift=FALSE)
  mort.e0.k = mort.k(inc, p.a, p.s, m.e, k, p.e0.t, h,
                     advanced=FALSE, shift=FALSE)
  # Intervention
  mort.a1.k = mort.k(inc, p.a, p.s, m.a, k, p.a1.t, h,
                     advanced=TRUE, shift=TRUE)
  mort.e1.k = mort.k(inc, p.a, p.s, m.e, k, p.e1.t, h,
                     advanced=FALSE, shift=TRUE)
  # Denominator and numerator of MRR 
  denominator = mort.a0.k + mort.e0.k
  numerator = mort.a1.k + mort.e1.k
  # Results
  dfres <- data.frame(Year=k,
                      Statistic=c('Mort',
                                  'Women Alive',
                                  'Surv',
                                  'MRR',
                                  'ARR',
                                  'Years of Life Saved'),
                      Control=c(denominator,
                                inc-mort.a0.k-mort.e0.k,
                                (inc-mort.a0.k-mort.e0.k)/inc,
                                1,
                                0,
                                0),
                      Intervention=c(numerator,
                                     inc-mort.a1.k-mort.e1.k,
                                     (inc-mort.a1.k-mort.e1.k)/inc,
                                     round(numerator/denominator,3),
                                     round(denominator-numerator,3),
                                     round((k/2)*(denominator-numerator),3))
  )
  return(dfres)
}
