rm(list=ls())
library(EMC2)

# This script introduces the 4 basic evidence-accumulation models (EAMs)
# supported by EMC2. 
#
# # Most of the sampling here is relatively quick as we fit only a small data set
# # from a single subject, although the DDM is somewhat slower. To work through
# # without having to run the sampling for all models load pre-computed samples
#  
# # print(load("BasicEAMs/FlankerSamples.RData"))

# First, as in 1-BasicEAMs.R, the data needs to be re-loaded.
print(load("BasicEAMs/FlankerData.RData"))
# Format the data to be suitable for analysis by the EMC2 package.
dat <- data_ajh
names(dat)[1] <- "subjects"
dat <- dat[,c("subjects","S","CI","E","R","rt")]

# There were some fast responses due to double button presses removed by the
# following line (these are all ~50ms which is the measurement resolution of
# the toy experiment written in RStudio). 
dat <- dat[dat$rt>.2,]

#### The "full" DDM ----

# First let's fit the DDM by also estimating the three between-trial variability
# parameters (they were previously set to constants).

# As before, we will assume a conventional model in which E selectively affects 
# thresholds and CI selectively affects rates.

designDDM <- design(model=DDM,data=dat,
  formula=list(a~E,v~0+S/CI,Z~1,t0~1,st0~1,sv~1,SZ~1)
)

# The priors for the extra parameters are again based on Matzke & Wagenmakers.
pmean = c(a=log(0.17),a_Espeed=0,
          v_Sleft=-2.25,v_Sright=2.25,
         'v_Sleft:CIincongruent'=0,'v_Sright:CIincongruent'=0,
         t0=log(.45),Z=0,st0=qnorm(.03),sv=log(1.2),SZ=qnorm(0.38))
psd = c(.7,.5,2.5,2.5,1,1,.4,.4,0.9,0.7,.75)
priorDDM <- prior(designDDM,pmean=pmean,psd=psd,type="single")
plot_prior(priorDDM,designDDM,layout=c(2,6))

# Fitting is much slower (~20x) than for the WDM because numerical integration
# over the st0 and SZ parameters is required to obtain the likelihood (although 
# EMC2 has its own C++ code for this likelihood see ?ddiffusion from the rtdists 
# package, which as the basis of EMC2's version).

# Hence, dropping st0 and/or SZ speeds fitting greatly. The former is often not 
# required (and we will drop it in later lessons) and the latter only when 
# errors are faster relative to correct responses (as occurs when fast 
# responding is emphasized, as we look at this manipulation we continue to 
# estimate this parameter throughout). 
sDDM <-  make_emc(dat,designDDM,type="single",rt_resolution=.05,prior=priorDDM)
sDDM <- fit(sDDM)

# Although slower, sampling works quite well thanks to EMC2's robust and 
# efficient sampler.
check(sDDM)

# Looking at chain plots we see t0 shows some evidence of bi-modality, with
# sampling jumping between a lower and higher mode.
plot(sDDM)

# The bi-modality is even more evident in a pairs plot. This makes 
# interpretation of t0 estimates difficult. 
pairs_posterior(sDDM)

# It is also clear that the data do not strongly constraint on the v, sv and SZ 
# parameters.
plot_pars(sDDM,layout=c(2,6))
plot_pars(sDDM,layout=c(2,6),use_prior_lim=FALSE,map=TRUE)

# We will look at goodness of fit after introducing the other EAMS. As the 
# basics of posterior parameter testing were covered in the last lesson, we will
# largely leave posterior testing of the parameters of the models fit in this
# lesson as an exercise, but make a few comments on parameter values to 
# provide some orientation.

# Briefly looking at DDM parameter estimates relative to the WDM we see rates are 
# larger and and thresholds greater. sv and SZ being greater than 
# zero causes extra errors, so to maintain the same accuracy faster rates and 
# higher thresholds are required.
posterior_summary(sDDM)

#### Race models ----

# EMC2 provides three "race" models, the Racing Diffusion Model (RDM), the 
# Linear Ballistic Accumulator (LBA) and the Lognormal Race (LNR).

# For race models the Rlevels argument automatically makes available a "latent
# response" (lR) factor that indexes each accumulator by its corresponding
# response. 

# We discuss below how lR can be used to allow for response bias as race models  
# do not have a separate response-bias parameter as does the DDM (i.e., Z). 

# Note that unlike the DDM, race models can accommodate more than two responses 
# by having more than two Rlevels, although here we focus on only the binary 
# response case.

# A new argument to design used with race model, matchfun, is supplied with 
# a function that indicates which accumulator (indexed by the lR factor) 
# corresponds to the a correct response for different stimuli. 

# In simple designs where the levels of a "stimulus" factor (S here) and the
# levels of lR (i.e., the levels of the R factor) are the same the matchfun
# has the following form.
matchfun=function(d)d$S==d$lR

# This automatically makes available a second way of addressing the 
# accumulators, a "latent match" factor (lM) with levels that can be used to 
# set different rates for matching accumulators (i.e., corresponding to a 
# correct response in a given condition) and mismatching accumulators
# (corresponding to incorrect responses).

# Race models have one accumulator for each response, each with its own 
# parameters. For rates it is useful to recode these parameters in terms of
# the average rate across accumulators (intercept), and d = difference 
# between accumulators. 

# This re-coding is achieved using a contrast matrix that is supplied to the
# contrasts argument of design, so we create it here for use in the race
# models examined below. The values of 1/2 and -1/2 are chosen so the "d" 
# parameter magnitude corresponds to the magnitude of the difference between
# accumulators (AD = average/difference.
ADmat <- matrix(c(-1/2,1/2),ncol=1,dimnames=list(NULL,"d")); ADmat

# When applied to rates this contrast is typically associated with the lM factor.
# In this case, d is positive when responding is above chance and equals 
# the match - mismatch rate. The latter is analogous to the WDM/DDM rate. The 
# intercept is an index of decision urgency and stimulus magnitude effects. 

# We will call the d parameter rate "quality" as it indicates how well a
# participant can discriminate between correct and incorrect responses. We will
# call the intercept either "urgency" or more generically rate "quantity" as
# it quantifies the overall level of evidence driving a response

# When applied to thresholds (B), the ADmat contrast is typically linked to
# the lR factor. d corresponds to response bias (i.e., the difference between 
# thresholds for each accumulator. Here that would be right - left thresholds
# (as the R factor has levels left and right), and so d < 0 corresponds to a 
# bias to respond left and d > 0 a bias to respond right, with d = 0 unbiased.

# Although we generally recommend this parameterization of thresholds here we
# use the default contrast for B in order to illustrate differences between
# different types of contrast coding.

# NB: Beyond being able to address choices among more than two options, an 
#     advantage of race models is that they can separate urgency/quantity and 
#     quality effects. However, this can make them harder to estimate and lead
#     to large posterior correlations as we illustrate below.


### RDM ----

# First we fit the RDM. All RDM parameters are positive and so estimated on a
# log scale. Negative rates, although conceptually possible, mean that sometimes
# and accumulator may never reach its threshold, and it also invalidates the
# analytic solution for the likelihood, requiring a much slower calculation.

# As in the DDM we set moment-to-moment variability to one and we also use the 
# simple form with no start-point variability (A) by fixing it to zero (although
# there is a fast analytic solution when A>0 so fitting is still fast). 

# As before, we will assume a conventional model in which E selectively affects 
# thresholds and CI selectively affects rates.

# Matching the WDM/DDM models we assume:
# a) Thresholds (B) are affected by E and allow for response bias that is the 
#    same for speed and accuracy by an additive combination with lR. 
designRDM <- design(data=dat,model=RDM,
                  formula=list(B~E+lR),
                  matchfun=function(d)d$S==d$lR,
                  contrasts=list(lM=ADmat))

# Here the intercept (B) estimates the threshold in accuracy-left, which acts as
# a baseline condition.

# B_Espeed is the change from accuracy to speed (equally for both left and
# right responses due to the additive linear model). Hence, it quantifies the
# effect of speed emphasis, with negative values indicating a lower threshold
# in speed and positive values a higher threshold for speed. The conventional
# expectation is the former.

# B_lRright is the change from left to right (equally for speed and
# accuracy conditions due to additivity). Hence, it quantifies the bias effect
# (negative values indicating a bias to left and positive a bias to right).

# Given accuracy-left is the baseline note how both B_Espeed and B_lRright are 
# added to determine speed-right.

# b) We assume rates vary with CI, both in terms of the difference between 
#    match and mismatch (analogous to the DDM) and in terms of urgency (which
#    does not have an analogue in the standard DDM). Here we use the "*" 
#    operator, which is shorthand for CI + lM + CI:lM.
designRDM <- design(data=dat,model=RDM,
                  formula=list(v~lM*CI),
                  matchfun=function(d)d$S==d$lR,
                  contrasts=list(lM=ADmat))

# In this case the v parameter is the average rate (over match and mismatching
# accumulators, i.e., urgency or quantity) for the congruent condition.

# v_CIincongruent is the incongruent condition effect on urgency/quantity (i.e.,
# the difference from the congruent condition).

# v_lMd is the average quality for the congruent condition, here we would expect
# a positive value because congruent accuracy is well above chance. 

# v_lMd:CIincongruent is the incongruent condition effect on quality 
# (i.e., the difference from the congruent condition). Here we would expect a
# negative value because accuracy is lower in the incongruent condition.

# You can see the last contrast is the product of the preceding two, which is
# typical for interaction (":") parameters.

# c) We assume the same non-decision time (t0) for all conditions, and so the full
# design is:

designRDM <- design(data=dat,model=RDM,
                  formula=list(B~E+lR,v~lM*CI,t0~1),
                  matchfun=function(d)d$S==d$lR,
                  contrasts=list(lM=ADmat))

# NB: By assigning ADmat to lM as is done in the contrasts argument above 
#     means that lM will use ADmat for every parameter. In the hierarchical
#     lesson, we will show how to associate different design matrices to the same 
#     factor for different model parameters.


## More on priors ----

# Priors are based on our informal experience with this model to be 
# be reasonable but fairly vague.
pmean <- c(B=log(2),B_Espeed=log(1),B_lRright=log(1),v=log(2),
  v_lMd=log(2),v_CIincongruent=log(1),'v_lMd:CIincongruent'=log(1),t0=log(.2))

# With the aid of plot_prior we choose standard deviations so that the prior 
# distributions cover a reasonable range.
psd <-  c(1,.5,.5,1,.5,.5,.5,.5)
priorRDM <- prior(designRDM,pmean=pmean,psd=psd,type="single")
plot_prior(priorRDM,designRDM)

# Note that the priors have long tails, but these are truncated here 
# automatically by plot_prior to make distributions easy to visualize.

# Particularly when parameters are not on the natural scale it can be hard to
# anticipate how your choices of prior means map to the design. The following 
# function helps you do that.
mapped_par(pmean,designRDM)

# For the rate parameters, for example in the match (TRUE) case we sum the 
# intercept (v) and add half of the quality (v_lMd) then exponentiate to get
# back to the natural scale.
exp(log(2) + log(2)*0.5)

# For the mismatch (FALSE) case we subtract half of v_lMd
exp(log(2) + log(2)* -0.5)

# Suppose you wished the prior to reflect findings that thresholds are lower
# in speed than accuracy conditions by setting B_Espeed=log(.75), i.e., 25%
# less than B in the accuracy condition.
pmean['B_Espeed'] <- log(.75)
mapped_par(pmean,designRDM)

# For the accuracy conditions we see as before that the threshold is given by 
# exponentiating the intercept (B)
exp(log(2))

# But now for the speed condition adding B_Espeed produces a 25% mlower value
exp(log(2) + log(.75))

# Although it is useful to understand these basic computations, when 
# parameterizations become more complicated it is easier and less error prone
# to simply use mapped_par to guide your choices.

## Fitting ----

# Returning to fitting we see that it is fairly fast.
sRDM <-  make_emc(dat,designRDM,type="single",rt_resolution=.05,prior=priorRDM)
sRDM <- fit(sRDM)

# The results look good.
check(sRDM)
plot(sRDM)

# Priors are dominated and estimates are fairly well localized, although less
# so for the "FALSE" rates due to the low error rate. As is typical, t0 
# estimates are less than the DDM.
plot_pars(sRDM,layout=c(2,4))
plot_pars(sRDM,layout=c(2,5),use_prior_lim=FALSE,map=TRUE)

# There are quite strong correlations between the v parameters. As for the DDM
# there is a negative correlation between the B and t0 parameters, but no 
# bi-modality in t0. 
pairs_posterior(sRDM)

### LBA ----

# The LBA is very similar to the RDM, except that
# 1) A is always fit
# 2) between-trial rate variability (sv) replaced with trial rate variability (s)
# 3) Mean rates (v) are unbounded (i.e., can be negative) so are estimated on 
#    the natural scale. All other parameters are on a log scale.

# NB: Although the mean of v can be negative, the default setting for the LBA
#     uses a normal rate distribution truncated below at zero. As a consequence
#     rates on any given trial are always positive even when v is negative. 
#     This assumption ensures that the accumulators must eventually reach their
#     threshold so a response is guaranteed. Otherwise all rates on a trial 
#     might be negative and so no response would occur.

# By convention sv varies with lM, with the intercept (v) fixed for 
# identifiability. This typically results in sv for match being less than 
# mismatch.
designLBA <- design(data=dat,model=LBA,
                    matchfun=function(d)d$S==d$lR,
  formula=list(B~E+lR,v~lM*CI,A~1,sv~lM,t0~1),
  constants=c(sv=0),contrasts=list(lM=ADmat))

# Again, priors are based on our informal experience with this model to be 
# be reasonable but fairly vague.
pmean <- c(B=log(2),B_Espeed=log(1),B_lRright=log(1),v=3,
  v_lMd=2,v_CIincongruent=2,'v_lMd:CIincongruent'=0,A=log(.5),sv_lMd=log(1),t0=log(.2))
psd <- c(1,.25,.25,3,3,3,3,1,1,.5)

priorLBA <- prior(designLBA,pmean=pmean,psd=psd,type="single")
plot_prior(priorLBA,designLBA,layout=c(3,4))
           
# Note that when A is estimated the threshold is b = B+A (this is also true for
# the RDM). This is done so that b > A (i.e., you cant start above the 
# threshold) is easy to enforce by requiring B and A to be positive. 

# The addition is calculated after the mapping, so to see the implied prior we
# must look on the natural scale and ask plot_prior to also do the calculation.
# Both are achieved with the add_recalculated argument. We see the b prior
# is a little broader than the B prior due to the addition of A.
plot_prior(priorLBA,designLBA,layout=c(3,6),add_recalculated = TRUE)

# Fitting is quick and works well
sLBA <-  make_emc(dat,designLBA,type="single",rt_resolution=.05,prior=priorLBA)
sLBA <- fit(sLBA)
check(sLBA)
plot(sLBA)


#  We see that contraction is poor for rate parameters, particularly quality.
plot_pars(sLBA,layout=c(2,5))

# Looking at the mapped results on the natural scale we see that the mismatch
# rate is very uncertain, reflecting the lack of error responses that do most
# to constrain it. 
plot_pars(sLBA,layout=c(2,6),use_prior_lim=FALSE,map=TRUE)


# As commonly occurs when error rates are low, there are strong correlations 
# between the different v parameters. As with the LBA B and t0 are negatively 
# correlated. As is common with the LBA t0 is estimated as quite small.
pairs_posterior(sLBA)


### LNR ----

# Finally we fit the LNR model, a very simple race model that does not separately
# identify rates and thresholds. As in the LBA, between-trial variability is 
# conventionally allowed to vary with lM but nothing needs to be fixed for 
# identifiability. 

# Here we fit a model where the Lognormal mean (m) is an additive combination
# of lM, CI and E main effects and interactions between lM and CI and lM and E 
# (analogous to the d effects in the RDM and LBA). 

# Note that larger values of m result in slower performance, so they can be
# thought of as being directly proportional to thresholds and inversely
# proportional to rates.

# This means that if we interpret "d" parameters as analogous to quality they 
# will be negative when accuracy is above chance (as then the value for match
# should be LESS than for mismatch). To make interpretation easier we flip the 
# sign of the ADmat (i.e., we use -ADmat) so that 
designLNR <- design(data=dat,model=LNR,
  matchfun=function(d)d$S==d$lR,
  contrasts=list(lM=-ADmat),
  formula=list(m~lM+CI+E+lM:CI+lM:E,s~lM,t0~1))

# The m parameter is the mean over accumulators in the baseline congruent-accuracy 
# condition.

# m_CIincongruent is the difference between congruent and incongruent conditions.
# As responding is slower in the incongruent condition, this is expected to be
# negative.

# m_Espeed the difference between speed and accuracy conditions. As responding is 
# slower in the accuracy condition, this is expected to be negative.

# m_lMd is analogous to rate quality, where bigger implies better discrimination.
# It estimates quality in the baseline congruent-accuracy condition.

# m_lMd:CIincongruent is the difference in quality between congruent
# and incongruent conditions. As accuracy is less in the incongruent condition
# this parameter is expected to be negative. 

# m_lMd:Espeed is the difference in quality between accuracy and speed 
# conditions. As accuracy is less in the speed condition this parameter is 
# expected to be negative. 

# Priors are based on our informal experience with this model to be 
# be reasonable but fairly vague. 
pmean <- c(m=-.5,m_lMd=1,m_CIincongruent=0,m_Espeed=0,
  'm_lMd:CIincongruent'=0,'m_lMd:Espeed'=0,s=log(1),s_lMd=log(1),t0=log(.2))
psd <- c(2,2,2,2,2,2,1,1,1)
priorLNR <- prior(designLNR,pmean=pmean,psd=psd,type="single") 
plot_prior(priorLNR,designLNR,layout=c(3,4))

# Fitting is very quick as the LNR likelihood is very fast to compute, 
# and convergence is excellent.
sLNR <-  make_emc(dat,designLNR,type="single",rt_resolution=.05,prior=priorLNR)
sLNR <- fit(sLNR)
check(sLNR)
plot(sLNR)

# Contraction is quite good except for the quality (d) parameters, reflecting a
# lack of constraint due to low error rates. t0 is estimated at a similar level
# to the RDM.
plot_pars(sLNR,layout=c(2,5))
plot_pars(sLNR,layout=c(2,6),use_prior_lim=FALSE,map=TRUE)

# Consistent with our observations about v parameters in the other race models,
# there are high correlation among m parameters. 
pairs_posterior(sLNR)

### More on posterior testing ----

# One aspect of posterior testing we had not covered previously is directly
# comparing estimated parameters. For example, suppose we wanted to test 
# whether the conflict effect for quantity (m_CIincongruent) is bigger than the 
# speed effect (m_Espeed).

# We could do this with the "fun" method described previously, but it is
# easier to just supply parameter names
credible(sLNR,c("m_Espeed","m_CIincongruent"))

# In this case the difference is not credible. Note that it is not possible to 
# make the same test with the hypothesis function, so we cannot test equality
# this way, but the fun method still works. This test provides positive evidence
# for no difference.

# As another example, for quality the conflict effect is credibly larger in 
# magnitude, with strong evidence for the difference.
credible(sLNR,c("m_lMd:Espeed","m_lMd:CIincongruent"))
hypothesis(sLNR,fun=\(x) diff(x[c("m_lMd:Espeed","m_lMd:CIincongruent")]),selection="alpha", map = FALSE)

#### Comparing models and model fit ----

# EMC2 provides three methods of comparing models, all of which take account
# of both goodness-of-fit and simplicity. 

# Two, DIC and BPIC, are relatively insensitive to the prior, with BPIC having
# twice the complexity penalty of DIC and so preferring simpler models.

# The third, marginal deviance (MD) is more dependent on the prior. Bayes 
# Factors are MD ratios. 

# All model selection methods favor the LBA, with the LNR second. 
compare(list(WDM=sWDM,DDM=sDDM,RDM=sRDM,LNR=sLNR,LBA=sLBA))

# Note that this tables also provides measures of misfit (without any complexity
# penalty), minD and meanD. Use the smallest of the two (they are often very 
# similar as here), and smaller indicates a better fit. In this case the LNR
# fits best. Note, however, that these are rather approximate measures so are
# to be used as only a rough guide.

# The table also provides "EffectiveN", which counts parameters while taking 
# account of correlations among them to come up with something analogous to ESS
# for samples. It identifies the DDM as being most complex, but again is rather
# approximate and can give negative values which are hard to interpret, as it 
# does for the LBA here.

# To examine model fit we calculate post-predictives.
ppWDM <- predict(sWDM,n_cores=4)
ppDDM <- predict(sDDM,n_cores=4)
ppRDM <- predict(sRDM,n_cores=4)
ppLNR <- predict(sLNR,n_cores=4)
ppLBA <- predict(sLBA,n_cores=4)

# Although the improvement in fit across models is clear even the best have
# some misfit evident.
plot_fit(dat,ppWDM,layout=c(2,4),lpos="bottomright",xlim=c(.4,1.4))
plot_fit(dat,ppDDM,layout=c(2,4),lpos="bottomright",xlim=c(.4,1.4))
plot_fit(dat,ppRDM,layout=c(2,4),lpos="bottomright",xlim=c(.4,1.4))
plot_fit(dat,ppLNR,layout=c(2,4),lpos="bottomright",xlim=c(.4,1.4))
plot_fit(dat,ppLBA,layout=c(2,4),lpos="bottomright",xlim=c(.4,1.4))

# It is often hard to judge the fine-grained details of misfit from CDF plots.
# To address this we can also focus on how well the model accounts for effects 
# on particular statistics like the mean and SD of RT and error rate.

# This can be done by supply plot_fit with a function to calculate each 
# statistic (the same function you would use on the data) through the stat
# argument.
mrt <- function(d)mean(d$rt)
sdrt <- function(d)sd(d$rt)
err <- function(d)mean(d$R!=d$S)

# The resulting plots specify the data value as a vertical line superimposed
# on the posterior predictive distribution of the statistics.

# In contrast to CDF fit plots, the default is to treat all of the data as 
# one. Here we use this to compare all models. The race models all look similar,
# whereas the WDM is clearly inadequate and the DDM in between
par(mfcol=c(3,5))
plot_fit(dat,ppWDM,stat=mrt,main="WDM",xlim=c(.64,.72)); plot_fit(dat,ppWDM,stat=sdrt)
plot_fit(dat,ppWDM,stat=err)
plot_fit(dat,ppDDM,stat=mrt,main="DDM",xlim=c(.64,.72)); plot_fit(dat,ppDDM,stat=sdrt,xlim=c(.08,.15))
plot_fit(dat,ppDDM,stat=err,xlim=c(0,.09))
plot_fit(dat,ppRDM,stat=mrt,main="RDM",xlim=c(.64,.72)); plot_fit(dat,ppRDM,stat=sdrt,xlim=c(.08,.15))
plot_fit(dat,ppRDM,stat=err,xlim=c(0,.09))
plot_fit(dat,ppLNR,stat=mrt,main="LNR",xlim=c(.64,.72)); plot_fit(dat,ppLNR,stat=sdrt,xlim=c(.08,.15))
plot_fit(dat,ppLNR,stat=err,xlim=c(0,.09))
plot_fit(dat,ppLBA,stat=mrt,main="LBA",xlim=c(.64,.72)); plot_fit(dat,ppLBA,stat=sdrt,xlim=c(.08,.15))
plot_fit(dat,ppLBA,stat=err,xlim=c(0,.09))


# To understand these results in more details we will focus on the DDM and the
# the LBA (i.e., the best models of each overall type) and break the fit results
# down by the two key manipulations, E and CI. 


# The LBA provides a fairly good fit to all conditions. 
par(mfrow=c(3,4))
plot_fit(dat,ppLBA,stat=mrt,factors=c("CI","E"),xlim=c(.57,.8),main="Mean RT")
plot_fit(dat,ppLBA,stat=sdrt,factors=c("CI","E"),xlim=c(.04,.2),main="SD RT")
plot_fit(dat,ppLBA,stat=err,factors=c("CI","E"),xlim=c(0,.2),main="Error Rate")

# The DDM struggles to match variability in the incongruent accuracy condition,
# but otherwise does a good job.
par(mfrow=c(3,4))
plot_fit(dat,ppDDM,stat=mrt,factors=c("CI","E"),xlim=c(.57,.8),main="Mean RT")
plot_fit(dat,ppDDM,stat=sdrt,factors=c("CI","E"),xlim=c(.04,.2),main="SD RT")
plot_fit(dat,ppDDM,stat=err,factors=c("CI","E"),xlim=c(0,.2),main="Error Rate")

# #  If you save your results choose a different file name.
# save(dat,sWDM,sDDM,sLNR,sRDM,sLBA,
#      ppWDM,ppDDM,ppRDM,ppLNR,ppLBA,file="BasicEAMs/FlankerSamples.RData")

#### Exercises ----

# 1) Try dropping the st0 parameter from fitting. 
# a) By what factor does it increase the speed of fitting?
# b) Use model selection to test the support for including st0.
# 
# 2) The RDM did less well than the LBA, but had two disadvantages, 
# i) it did not allow s to vary with lM (analogous to sv~lM in the LBA) 
# ii) it did not estimate A
# 
# a) Fit a model addressing (i), how long does it take relative to the original?
# b) Fit a model addressing (i) and (ii), how long does that take?
# c) Compare the three models and the LBA, which wins?
# d) Which is more beneficial appears more beneficial, adding s~lM or adding A~1
# as well?

# 3) Both RDM and LBA suffer from problems with identifying mismatch rates due
# to low error rates. Try to address this by setting the mismatch rate to a
# constant value of 1 and freely estimate the match rates as before. For the RDM
# use the sv~lM. Note that this makes the race models more like the DDM, as an
# increase in match v increases quantity and quality together.
# 
# a) Which model wins between the original and new RDM? 
# b) How is the fit of the new RDM model affected?
# c) How are posterior correlations and updating in the new RDM model compared 
#    to the old model?
# d) Which model wins between the original and new LBA?
# e) How is the fit of the new LBA model affected?
# f) How are posterior correlations and updating in the new RDM model compared 
#    to the old model?
#
# 4) Compare the mapped parameters for the new models that are analogous.
# a) Non-decision time
# b) Rate variability
# c) Thresholds
# d) Mean rates
# e) How does the picture of performance provided by each model differ?


