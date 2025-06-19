rm(list=ls()) # clear R environment

# This script introduces methods to work with evidence-accumulation 
# models (EAM) using as an example a simple EAM for binary choice, the Wiener 
# Diffusion Model (WDM)
#
# # Most of the sampling here is relatively quick as we fit only a small data set
# # from a single subject, although the DDM is somewhat slower. To work through
# # without having to run the sampling for all models load pre-computed samples
#  
# # print(load("BasicEAMs/FlankerSamples.RData"))

#### Getting started ----

# # If not already done, install the EMC2 package
# remotes::install_github("ampl-psych/EMC2",ref="dev")

library(EMC2)
# This script analyzes the data collected by one of the authors (ajh) in a 
# conflict task (the Flanker task) with an added speed vs. accuracy emphasis
# manipulation. The task is described in more detail, and instructions provided
# on how to collect your own data in this task, in the BasicEAMs/FlankerTask.pdf 
# document.

### Load the data ---- 

print( load("BasicEAMs/FlankerData.RData"))

# Format the data to be suitable for analysis by the EMC2 package.
dat <- data_ajh
names(dat)[1] <- "subjects"
dat <- dat[,c("subjects","S","CI","E","R","rt")]

# There were some fast responses due to double button presses removed by the
# following line (these are all ~50ms which is the measurement resolution of
# the toy experiment written in RStudio). 
dat <- dat[dat$rt>.2,]

# NB: Your data may also contain very slow outlying responses, check for and
#     remove these if you are analyzing your own data.

# We will now fit some standard evidence accumulation models to this data, the
# Diffusion Decision Model (DDM, and its simplified version the Wiener Diffusion
# model, WDM), and three race models, the LBA (Linear Ballistic Accumulator), 
# RDM (Racing Diffusion Model) and LNR (log-normal race). 

# The expectation is that these models, although appropriate for standard 
# decision tasks with speed vs. accuracy manipulations, will have shortcomings 
# with conflict data, enabling an illustration of how to critically evaluate
# a cognitive model. 


### Specifying a design ----

# A design and associated model are specified in EMC2 using the "design" 
# function. We will now work through how to specify 5 core arguments (model,
# factors, Rlevels, formula, and constants)
# see:
?design()

# 1) model: The DDM model is specified with this argument:
model=DDM

# DDM is a function exported by the EMC2 package. The function is for
# internal use so you don't need to understand its structure now (although
# it will be addressed later in the course), but if you are curious you can look 
# at it by printing it.
DDM
# The help file of the function provides an overview of the parameters
?DDM()

# 2) factors: Two further arguments in design() provide information about the data to be fit. The 
# first, factors (factors that can be used in the formula for the model specification) 
# specifies a list of the levels of the factors spanning all of the cells in the data frame to be 
# analyzed (i.e., emphasis E, congruency CI and stimulus S)  which can have any name 
# except for the following names: 
#
# "R", "rt" and "subjects"
#
# These are reserved for columns in the data frame giving the response factor (R), 
# a numeric response  time column (rt), and a factor column specifying the subjects, 
# with levels that must also be specified in the factors argument. 

# You should also avoid factors called "trials" as that name is used internally by EMC2.

# In this case, given there is only one subject we have
factors=list(subjects=levels(dat$subjects),
              S=levels(dat$S),CI=levels(dat$CI),E=levels(dat$E))

# Note that instead of the levels() call, you can also supply a character vector.

# 3) Rlevels: we also separately specify the levels of the response factor in the Rlevels
# (response levels) argument. In more advanced applications, this gives EMC2 the
# ability to move beyond binary choice. 
Rlevels=levels(dat$R)

# 4) formula: design() specifies the mapping between model parameters and the different 
# types of data defined by factors ("cells") using R's linear modeling language 
# through a list of formulas, one for each model-parameter type, given as a list 
# to the formula argument. We will use a conventional characterization
# where:
# a) Thresholds (a) are affected by E (i.e., emphasis) with parameters a (the intercept),
#    corresponding to the threshold in accuracy conditions and a_Espeed, the
#    "speed minus accuracy" difference, expected to be negative as 
#    speed emphasis reduces thresholds.
formula=list(a~E)

# b) Rates (v), are also affected by S (stimulus) and CI (congruency) in a fully crossed manner. 
#    To do this in an easily interpreted way we: 
#
#    1) remove the intercept ("0+"), which allows us to have separate terms for the left and 
#    right stimulus rates (v_Sleft and v_Sright). The left rate is expected to 
#    be negative as left responses are mapped to the lower boundary and right 
#    positive as right responses are mapped to the upper boundary, and 
#
#    2) we use "nesting", S/CI ("CI nested within S"), which allows us to estimate the 
#    incongruent-congruent rate difference separately for the left and right stimuli
#    (v_Sleft:CIincongruent and v_Sright:CIincongruent). In 
#    both cases we might expect the incongruent rate to be less than the 
#    congruent rate, and hence these terms to be positive and negative 
#    respectively.
formula=list(a~E, v~0+S/CI)

# c) Non-decision time (t0) and response bias (Z, as a proportion of a, so
#    a value of 0.5 is unbiased) that is the same for all conditions. 
formula=list(a~E, v~0+S/CI, t0~1, Z~1)

# NB: We do not have to provide formulas for all 9 of the DDM's parameter types. 
#     If any model parameters are not named they are set at default values.
#     If parameters are given in the formula they can also be set to a fixed
#     value using the constants argument. For example, if we include the
#     the remaining 5 DDM parameters,
formula=list(a~E, v~0+S/CI, t0~1, Z~1, st0~1, sv~1, SZ~1, s~1, DP~1)

#     we can then set them to their defaults.
constants=c(s=log(1),st0=log(0),DP=qnorm(0.5),SZ=qnorm(0),sv=log(0))

# NB: Parameter types in EMC2 have a default transformation used during sampling,
#     which ensures they are unbounded, a log transformation for positive parameters 
#     e.g., variability parameters like a, t0, s, st0 and sv, and a probit 
#     transformation for doubly bounded parameters, like Z, SZ, and DP. Only
#     the mean rate (v) parameter is estimated on the natural scale.

# Both the default values and scales of model parameters are given in the help 
# file of the corresponding model ?DDM() in this case. 


######## Model Fitting 

# As a homework exercise after the class, apply the following analyses to your own
# collected Flanker data. Alternatively, you could do that during class while you follow along.

#### The WDM ----

# The WDM is a special case of the more generally used Diffusion Decision Model
# (DDM) without the latter's between-trial variability parameters (ie., they are
# set to zero: st0=log(0),sv=log(0),and SZ=qnorm(0)).

# As is conventional with the WDM and DDM, we also fix the standard deviation of 
# moment-to-moment variability (s) to 1 (as it is estimated on a log scale 
# s = log(1)).

# NB: We also assume response-production time for left and right responses is
#     the same by fixing DP (the difference is estimated on a probit scale, 
#     with equality at qnorm(0.5)). This parameter has only been estimated in
#     a few cases, and is likely only necessary when the force required to make
#     the two responses (and hence response production time) differs markedly.

### Specifying a design ----

# Rather than using the constants argument, we simply omit the fixed parameters 
# from the formula and so they take on the default values corresponding to the
# WDM. 

# Below we specify a conventional model in which E selectively affects 
# thresholds and CI selectively affects rates.

designWDM <- design(model=DDM,
                    Rlevels=levels(dat$R),
                    factors=list(subjects=levels(dat$subjects),
                      S=levels(dat$S),CI=levels(dat$CI),E=levels(dat$E)),
                    formula=list(a~E, v~0+S/CI, t0~1, Z~1)
                  )

# Note that the specification of Rlevels and factors can also be automatically
# taken from the data (i.e., dat). To do so, provide a data frame as the data argument of
# design. 

# To use this method, the only factors apart from R in the data frame must be the 
# design factors. Typically these design factors will not be redundant and
# together will define every design cell (8 in this case). The only
# non-factor column must be "rt" (unless there are covariates, see later lesson).

# We will usually use this compact design specification. 
designWDM <- design(model=DDM,data=dat,formula=list(a~E,v~0+S/CI,Z~1,t0~1))

# Design prints out the parameters that will be sampled during model fitting
# along with the design or "contrast" matrices used to map the parameters to the 
# data cells. Here we used R's default treatment contrast (see ?contr.treatment).

# Below we will discuss how you can pass different contrasts to EMC2. For now we
# explain the effect of the default choices. For the two parameter types that
# have more than one level (and so require a contrast).

## Thresholds ----
designWDMa <- design(model=DDM,data=dat,formula=list(a~E))

# The "a" parameter is the threshold for the accuracy condition. This is because
# only the accuracy row of the desgin matrix has an entry for a.

# The "a_Espeed" is the amount added to the accuracy condition threshold to 
# obtain the speed condition threshold. Such differences are usually called 
# "effects" (here the effect of speed instructions relative to accuracy 
# instructions).

# NB1: Because a is estimated on a log scale adding an effect has a 
#      multiplicative effect on the natural scale. For example if a = log(2)
#      and a_Espeed = log(1) = 0, then the speed condition threshold is 
#      exp( log(2) + log(1) ) = 2 x 1 = 2 (i.e, the same). Similarly if 
#      a_Espeed = log(.5) = -.69 then it is 2 x 1/2 = 1 (i.e., less). If
#      a_Espeed = log(2) = .69 then it is 2 x 2 = 4 (i.e., more).

# NB2: Generally EMC2 makes up parameter names by pasting the parameter type 
#      and an underscore (e.g., "a_"), with factor names + factor level (e.g., 
#      Espeed). Hence, to make names compact and readable we will typically use 
#      factors specified by one or two upper case letters and have longer 
#      but still compact descriptive names for their levels.

## Rates ----
# The design matrix for v is more complicated, as it is a function of two 
# factors. 
designWDMv <- design(model=DDM,data=dat,formula=list(v~0+S/CI))

# The v_Sleft and v_Sright parameters specify the congruent condition rates for
# left and right stimuli respectively (again you can see that only these 
# parameters have entries in those rows)

# v_Sleft:CIincongruent and v_Sright:CIincongruent specify the effect of incongruent
# relative to congruent, again for left and right stimuli respectively. The ":"
# here is read as "interacts with" or simply "and" in the R linear model 
# language.  

#### Specifying a prior ----

# In order to fit the model, we need to make a prior distribution for each parameter. In
# EMC2, priors are always normal (which is one of the reasons why parameters are 
# transformed to be unbounded).

# To see the parameter names for which priors need to be specified you can 
# extract the parameters that will be estimated from the design object.
names(sampled_p_vector(designWDM))

# The following prior for intercept parameters is based on the survey of DDM 
# fits by Matzke, D., & Wagenmakers, E. J. (2009). Psychonomic Bulletin & Review, 
# 16(5), 798–817.
#
# For effects (i.e., differences between conditions) we assume a mean of zero
# (i.e., no direction assumed) and standard deviations with a vague
# commitment to smaller values. Note that values are specified on the same 
# scale as sampling (i.e., unbounded). 

# A prior for a single-subject fit consists of two elements, prior means 
# (theta_mu_mean) and a variance-covariance matrix (theta_mu_var). Here we 
# assume that the covariances are zero.

# First specify the means
pmean = c(a=log(.17),a_Espeed=0,v_Sleft=-2.25,v_Sright=2.25,
  'v_Sleft:CIincongruent'=0,'v_Sright:CIincongruent'=0,t0=log(.45),Z=qnorm(.5))

# To avoid typing in names you could instead use
pmean <- sampled_p_vector(designWDM,doMap=FALSE)
pmean[1:length(pmean)] <- c(log(.17),0,-2.25,2.25,0,0,qnorm(.5),log(.45))

# For the standard deviations there is no need to provide names unless you want
# to do so, but make sure the order follows that of pmean.
psd <- c(.7,.5,2.5,2.5,1,1,.4,.4)

# We can use the "prior" function to combine this information, indicting that
# the prior is at the single subject level.
priorWDM <- prior(designWDM,pmean=pmean,psd=psd,type="single")

# We can use the "plot_prior" function to sample from and plot the prior.
plot_prior(priorWDM,designWDM,map=FALSE,layout=c(2,4))

# You can also save the prior samples instead of (just) plotting them. 
# The prior samples are returned as a coda mcmc.list object
priorSamples <- plot_prior(priorWDM,designWDM,map=FALSE,layout=c(2,4), do_plot = FALSE)
head(priorSamples[[1]])

# NB: In EMC2, alpha just means that the plotted the parameter is a subject-level parameter 
#     (and not e.g. a population-level parameter).

# It is more informative to plot the prior on the natural scale and mapped to
# the design cells over which they vary, which is the default (map=TRUE). 
plot_prior(priorWDM,designWDM,layout=c(2,4))

# NB: You will sometimes have to tinker with psd to make sure that the main part
#     of the prior encompasses a reasonable range on the natural scale, using 
#     the last command to guide your judgement graphically.

# We end by noting a somewhat un-intuitive property of probit transformed 
# parameters: choosing a large standard deviation can cause the prior to be 
# more informative rather than less informative. To illustrate, increasing
# the sd for Z to 4 makes extreme values (i.e., strong left or right 
# bias) more likely!
psd1 <- psd
psd1[7] <- 4 # changing the sd of Z
priorWDM1 <- prior(designWDM,pmean=pmean,psd=psd1,type="single")
plot_prior(priorWDM1,designWDM,use_par="Z")

# Instead, the maximally uninformative prior, when the mean is zero on the 
# probit scale (i.e, Z = 0.5 on the natural scale) is given by sd = 1, as that
# produces a uniform distribution on the natural scale.
psd1[7] <- 1
priorWDM1 <- prior(designWDM,pmean=pmean,psd=psd1,type="single")
plot_prior(priorWDM1,designWDM,use_par="Z")


#### Fitting ----

# Using make_emc(), we make an object to be used to perform fitting which combines the data, 
# design and prior. The sampling type is "single" because we are fitting a 
# non-hierarchical model with just one participant. 

# We set the rt_resolution argument (the accuracy with which rt is measured) to match the
# coarse grain of measurement (by default rt_resolution=.02 to 
# match a typical experiment using a 50Hz monitor with no synchronization with 
# the screen refresh) in our toy experiment. 

# This speeds fitting because the (often expensive) likelihood calculation is 
# done only once for each unique combination of  parameters and rt discretised 
# at the given resolution, with the speedup factor (relative to calculating 
# likelihoods for each trial) reported along with the number of unique trials 
# left after this compression.
#
# By default, three independent chains are sampled, which is useful in order to
# check whether sampling has succeeded, see below). 
sWDM <-  make_emc(dat,designWDM,type="single",rt_resolution=.05,prior=priorWDM)

# We can now run the fit. This goes through several stages and if all goes well
# gets a default 1000 converged iterations (3000 samples) in the final stage. 
#
# 1) "preburn", default 150 iterations with settings that try to quickly search
# for a parameter region with good fit.

# 2) "burn", the sampler runs repeatedly, adding 100 iterations then checking
# the rhat (gelman.diag) for each parameter until the mean is < 1.1 by default.
# Older samples are thrown away if that makes the rhat better. 
#
# 3) "adapt", in this stage the sampler estimates an approximate model of the
# posterior in order to make sampling more efficient (less autocorrelated) 
# in the next stage, which is evident by the end of this stage, with the chains
# looking more like "fat flat hairy caterpillars".
#
# 4) Finally in the sample stage the efficient sampling model is used to get
# the samples we will use for further analysis, with rhat < 1.1.

# By default the fitting function will report progress in each stage and the total
# run time. The default is also to use parallel processing, allocating one 
# core for each chain, so with the default number of chains 3 cores will be 
# used. As the chains are independent this (almost) triples the speed of fitting.
sWDM <- fit(sWDM)

# NB1: You can adjust the number of cores used with the cores_for_chains
#     argument, but will get no speed-up using more cores than there are chains.
#     Each extra chain inherits some inefficiency as everything before the 
#     sample stage is not used in further analyses.

# NB2: Once a fit is done you can add more samples by simply running it again.
#      It will automatically start from where it left off (including if it 
#      stopped before finishing the sample stage). The number of extra 
#      samples added in the sample stage is determined by the iter argument 
#      (default 1000).

#### Inspecting sampling performance ----

# The following function checks whether the final sampling stage has converged,
# using "Rhat" values from coda (see ?gelman.diag). In its basic form Rhat
# test whether the chains are "mixed". EMC2 also splits chains into the first
# and second half before calculating Rhat (so by default it is based on 6 
# chains), which makes it sensitive to "stationarity", that is, the distribution
# of a chain being unchanging (flat) over the course of sampling (e.g., a 
# constant mean, variance etc.).

# Check() also reports the effective sample size (ESS) for each parameter 
# (see ?effectiveSize from coda), the number of independent samples obtained 
# after taking account of autocorrelation (typically less than the nominal 
# number of samples, ere default 3 x 1000). 


# By default check provides a chain plot of the parameter with the worst Rhat.
# Each chain (set of mcmc samples) is shown in a different colour. 
check(sWDM)

# With practice you will be able to recognize when chains are converged, with
# each chain looking like a "fat flat hairy caterpillar" (fat because most 
# samples are in the middle, hairy because the distributions have tails). 

# To look at chain plots for all parameters we can simply plot the samples.
# This provides a good example of fat,flat, hair catapillars and good chain
# mixing.
plot(sWDM, layout = c(2,4))

# We can see why EAMs are "sloppy" models (i.e., their parameters are often
# highly correlated) with a "pairs" plot. Some correlations are due to the
# way the model is parameterised (e.g., v_Sright and v_Sright:CIincongruent),
# others are more structural, e.g., higher thresholds (a) are associated with 
# lower non-decision time (t0)). 
pairs_posterior(sWDM)

# NB: The other models examined here often have even higher correlations, making 
#     them difficult for standard samplers to work with, although sloppiness is 
#     not itself an intrinsically bad thing.

# We can also plot the posterior (black) and prior (red) densities together to
# see how much we have learned from the data (i.e., how much the information 
# provided by the data has "updated" the prior). Visually updating is indicated
# by the posterior "dominating" (being much higher) and "contracting" (being
# less variable) than the prior.
plot_pars(sWDM,layout=c(2,4))

# A contraction statistic is also provided, with values near 1 indicating strong 
# updating of the prior. In this case we see that contraction was least for the
# rate parameters, although it is still quite high. Also, we see our prior 
# greatly underestimated the absolute magnitude of the posterior (i.e., rates
# for congruent were much higher, and the effect of conflict on reducing rates
# much greater) than assumed by the prior.

# By default this plot is on the scale of the prior to check whether the priors 
# encompass the posterior and hence are unlikely to have been very influential.
# Here we see that this was not the case for v_Sleft and v_Sright, but that the
# data largely "overwhelmed" the prior, so that this is likely to have lead to
# only a small bias.

# We can also use map = TRUE to map the posterior samples back to the cells of
# the experimental design and to the natural scale (the two must occur together),
# which makes interpritation easier.
plot_pars(sWDM,layout=c(2,4), map = TRUE)

# In this case (with scales that make it hard to see in detail) it can 
# also be useful to look at the plot on the scale of the posterior. Setting
# use_prior_lim = FALSE, we see that all estimates are fairly well localized.
plot_pars(sWDM,layout=c(2,4),use_prior_lim=FALSE,map=TRUE)

# Credible intervals (by default 95%) can also be obtained to see how well the 
# estimates are localized in tabular form.
posterior_summary(sWDM,map=TRUE)

# The generic summary function provides the same information along with 
# convergence and efficiency information.
summary(sWDM)

#### Posterior testing ----

# The credible function acts much like a Bayesian version of the standard t.test 
# function. For example, we might ask whether responding is unbiased, i.e., 
# whether Z = 0.5 (although this example is somewhat redundant with the output
# of summary). 

# As expected, we see bias is not credible as the credible interval encompasses 
# 0.5. 
credible(sWDM,"Z",map=TRUE,mu=.5)

# The same test can be made on the probit scale, using the default null value of
# mu = 0. As well as providing a credible interval, the output reports a 
# probability for a directional test, by default the probability of samples 
# less than the null value. Suppose we had reason to believe bias would be 
# towards the right response (i.e., positive) then we might be interested in the
# probability of samples greater than zero. These results indicate that
# credibility of that directional hypothesis is low.
credible(sWDM,"Z",alternative="greater")

# Credible intervals cannot prove the null, but that can be done using a
# Savage-Dickey test which approximates a Bayes Factor by the ratio of the 
# heights of prior to posterior densities at the null value.

# A large ratio (prior >  posterior) rejects the NULL (i.e., after being informed
# by the data the NULL value is less likely than under the prior). 

# A small ratio (prior <  posterior) favors the NULL (i.e., after being informed
# by the data the NULL value is more likely than under the prior). 

# Here the latter holds, as is shown graphically by the "hypothesis" function,
# with the Bayes Factor (i.e., the ratio) being returned. Note that there is some
# variability because every time you run hypothesis(), new samples are generated
# from the prior.
layout(1)
hypothesis(sWDM,"Z",selection="alpha", map = FALSE)

# It is easier to interpret the ratio when greater than one, which can be
# obtained by inverting the output. One would conclude that unbiased 
# responding is more than 11 times more likely than biased responding. 
1/hypothesis(sWDM,"Z",selection="alpha", map = FALSE)

# Conventionally Bayes Factors of greater than 10 (or less than 0.1) are 
# considered "strong" evidence, values between 1/3 and 3 are considered
# "equivocal" and in between as providing "positive" evidence.

# We can also test hypotheses not implicit in the way the model is 
# parameterized. For example, suppose we wished to test if the (absolute)
# conflict effects on rates differ for left vs. right arrows.

# To do this we need to provide a function that calculates the necessary
# difference on a parameter vector. This function is then applied to the 
# posterior samples and the credible interval of the difference returned.
fun <- \(x) diff(abs(x[c("v_Sleft:CIincongruent","v_Sright:CIincongruent")]))

# Here the results are equivocal, i.e., we do not have clear evidence either 
# for no difference or for the trend for a greater effect for right.
credible(sWDM,x_fun = fun)
1/hypothesis(sWDM,fun = fun,selection="alpha", map =FALSE)

#### Goodness of fit ----

# In order to test whether the model fits the data we need to obtain "posterior
# predictive" samples. These are obtained by randomly selecting a posterior 
# sample and using it to simulate data in the same design as the empirical data.

# By default 100 data sets are obtained. More samples might sometimes be 
# needed to get clear indications about the uncertainty of the model's 
# predictions about the data.
ppWDM <- predict(sWDM)

# Global fit is shown by "defective" cumulative distribution function (CDF) 
# plots. A CDF shows the cumulative probability RT. The CDF is not defective
# (i.e., it asymptotes on the right at a probability of 1) when only one 
# response is made, as is the case here (i.e., perfect accuracy) for all but the 
# incongruent speed condition. 

# Observed data are shown in black and posterior predictives in grey and in
# both cases 5 points are shown corresponding to the 10th, 30th, 50th (median),
# 70th and 90th percentiles. For the posterior predictives the large grey dot
# is the average of the 100 samples. 

# The small grey dots show each replication as a means of indicating uncertainty.
# Note that in contrast to frequentist approaches, where uncertainty is in 
# terms of hypothetical future replications, in the Bayesian approach "the data
# are the data" (so shown with no uncertainty) and it is the model that is 
# uncertain (i.e., its predictions have a distribution).
plot_fit(dat,ppWDM,layout=c(2,4))

# For the incongruent condition there are two CDFs, one for correct and one for
# error responses, and so data CDFs asymptote at values below 1 (i.e., the 
# probability of correct and error responses).

# Generally we see that the WDM provides a poor model of the data. Data 
# quantiles often barely overlap with predicted quantiles and the asymptotic
# levels are very different in the incongruent speed condition.

# In the next lesson we will look at the four main EAMs provided by EMC2 to 
# see if they do better.

