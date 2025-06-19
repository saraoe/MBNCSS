rm(list=ls())
library(EMC2)
# This lesson describes how to go beyond the standard models offered by EMC2, 
# It describes how to implement and test new models, model parameteriations and 
# ways of bounding parameters, the latter being sometimes necessary when there 
# are convergence problems. It also gives insight into the inner working of EMC.

# We will use the Basic lesson data to illustrate
load("BasicEAMs/FlankerData.RData")
dat <- data_ajh
names(dat)[1] <- "subjects"
dat <- dat[,c("subjects","S","CI","E","R","rt")]
dat <- dat[dat$rt>.2,-7]

# We will use this contrast later
ADmat <- matrix(c(-1/2,1/2),ncol=1,dimnames=list(NULL,"d"))

# # Load pre-computed sample 
# load("NewModels/RDMt0nL.RData")


#### Making a new model ----

# Suppose that you have reason to believe that non-decision time (t0) is slower
# for right than left responses (e.g., due to the right response button being
# stiffer than the left response button) and that non-decision time is faster
# under speed than accuracy emphasis, and crucially for this illustration, that
# these two effects are additive on the natural (i.e., seconds) scale. For 
# example, suppose that the two effects have equal and opposite effects of 0.1s. 
# Additively implies the two effects exactly cancel for right (0.1s slower) 
# speed emphasis (0.1s faster) responses. 

# Unfortunately there is no easy way to make the effects additive using the 
# standard RDM model as it estimates t0 on a  log scale, and an additive model 
# (t0 ~ S + E) on the log scale is not additive on the natural scale.
#
# To illustrate consider the standard model with an additive t0 effect.

designRDMt0 <- design(data=dat,model=RDM,
  matchfun=function(d)d$S==d$lR,contrasts=list(lM=ADmat),
  formula=list(B~E+lR,v~lM*CI,A~1,t0~S+E,s~1),constants=c(s=log(1),A=log(0)))

# Make a parameter vector where the baseline t0 = 0.2, right is 0.1s slower
# (so on a log scale we add log(1.5) which is the same as multiply by 1.5 on 
# the natural scale, an increase of 0.1 (.2*1.5=0.3) and speed is 0.1s faster, 
# (so on a log scale we add log(0.5) which is the same as multiply by 0.5 on 
# the natural scale, a decrease of 0.1 (.2*0.5=0.1).

p_vectort0 <- c(B=log(2),B_Espeed=log(1),B_lRright=log(1),v=log(2),
                v_lMd=log(2),v_CIincongruent=log(1),'v_lMd:CIincongruent'=log(1),
                t0=log(.2),t0_Sright=log(1.5),t0_Espeed=log(0.5))

# Examining t0 for right speed we that rather than being 0.2 (i.e., the two
# effects cancel) it is 0.15 ( exp(log(.2) + log(1.5) + log(.5)) = .2*1.5*.5)
mapped_par(p_vectort0,designRDMt0)

# We will now make a new model to solve this problem, which we will call
# RDMt0natural. Before doing that we will describe a method to check on 
# whether a model has been implemented correctly and some of the inner workings
# of EMC2.

### Profile plots and augmented data  ----

# When implementing a new model it is easy to make coding mistakes. A useful 
# check is to see if random and likelihood functions (both of which use the 
# transforms) are still consistent by simulating a large amount of data data and 
# making likelihood "profile" plots. These plots hold all parameters at the data 
# generating values except one, which is varied to check if the likelihood peaks 
# above the data generating value.

# To illustrate we use the standard RDM to simulate large data set so that 
# sampling noise is at a minimum. The result is 40,000 trials, 5000 in each of 
# the 8 design cells.
sdat <- make_data(p_vectort0,designRDMt0,n_trials=5000)

# We see the profiles peak above the data-generating values. By default the 
# calculates the likelihood at 100 points and prints out a table containing the
# true values, the value for which the likelihood was maximum (a crude estimate
# of the best fitting parameter) and the difference between the two. We see
# this "miss" is very small in all cases, consistent with the model being  
# correctly specified. 
profile_plot(sdat, designRDMt0, p_vector = p_vectort0, layout = c(2,5),n_cores=12)

# Profile plot is initially quite slow because it must first create a "dadm"
# object using the EMC2's internal "design_model" function, which requires 
# several steps
# 1) data augmentation
# 2) data compression
# 3) binding the augmented and compressed data with a design and model

# For a race model data augmentation consists of stacking as many copies of the
# data together as there are accumulators in the model. For the binary RDM
# example here this doubles the size. To illustrate we first make a 
# dadm without compression and with only minimal rounding of rt's: 
dadm1 <- EMC2:::design_model(sdat,designRDMt0,compress=FALSE, rt_resolution=.001)
nrow(sdat)
nrow(dadm1)

# We see that each row of sdat is now repeated twice and three columns added:
# a) lR factor = latent response (with the same levels as R)
# b) lM logical = latent match, indicating which accumulator corresponds to the
#       the correct response given the stimulus (created by matchfun) and 
# c) winner logical = which accumulator corresponds to the response
head(dadm1)


# Next the data are compressed as determined by the rt_resolution parameter. 
# The default is 0.02 seconds, but here we use 0.05, so that only 1219 trials
# are kept that have the potential for having different likelihoods given any
# possible parameter vector. Calculating the compression can be time 
# consuming.
dadm2 <- EMC2:::design_model(sdat,designRDMt0,rt_resolution=.02)

# New now see that RTs are rounded to the nearest 0.05 and the dadm is much
# smaller (2 x 1219 rows)
nrow(dadm2)
head(dadm2)

# We see that if we now pass the dadm to profile plot it runs much more quickly
# as it does not have to do the time consuming compression operation. Although
# you mist still pass sdat this is no longer used.
profile_plot(sdat, designRDMt0, p_vector = p_vectort0, layout = c(2,5),
             dadm=dadm2,n_cores=12)

# Passing a dadm is the only way to examine the effect of compression on the
# profile. Generally we have found that compression of .05 or less has little
# impact on parameter recovery for race models.

### EMC2 models ----

# Now lets return to making a new model. Models in EMC2 are functions returning 
# a list of values and functions. Lets look at the RDM model and go through its 
# components that we have to change to make RDMt0natural. 
RDM

# Standard models like the RDM are also implemented in C++ through Rcpp (see the
# model_RDM.h in the package src directory), but this is only used when 
# sampling, where speed is critical. Otherwise (when simulating data etc.). 

# The new model will make will not have a C++ version. The presence of a C++ 
# version is signaled by the presence of the following list entry, so we will
# simply omit this entry from our new model
RDM()$c_name

# Note that even without this EMC2 still uses the C++ versions of functions
# to generate random samples from (rfun in the model list), and the density 
# (dfun) and cumulative density (pfun) of a single RDM accumulator defined at 
# the bottom of the R model function, it is only the likelihood function call 
# that is not in C.

# Relevant to the aim here, the Ntransform function, which transforms parameters
# back to the natural scale after they have been mapped. In this case this 
# involves exponentiation for all RMD parameters as they must all be positive,
# and so are all estimated on the log scale
Ntransform=function(x) {
  exp(x)
}

# In RDMt0natural we change this to:
Ntransform=function(x) {
  x[,dimnames(x)[[2]]  != "t0"] <- exp(x[,dimnames(x)[[2]]  != "t0"])
  x
}

# Note that x is a matrix of parameters with column names given by the 
# names of the p_types vector of the model function, which we leave unchanged.

# Being on the natural scale causes a potential problem that the log transform 
# was used to address, keeping all parameters positive after mapping. That is 
# no longer guaranteed, but is taken care of by a further model element that is 
# the same in both cases:
Ttransform = function(pars,dadm) {
  attr(pars,"ok") <- (pars[,"t0"] > .05) & ((pars[,"A"] > 1e-6) | pars[,"A"] == 0) &
                      ((pars[,"v"] > 1e-3) | pars[,"v"] == 0)
  pars
}
# Here "pars" is a matrix of parameters, one for each row of the dadm. For any 
# cases where ok=FALSE the likelihood return a value very close to to zero
# (by default 1e-10) which tends to make the sampler reject such parameters. 

# Here we see that the standard model already did not allow implausibly small 
# values of t0, and also very small values of start-point noise (which can 
# sometimes cause numerical problems in likelihood calculation), although A=0 
# is handled as a special case in the likelihood code so is allowed. 

# NB1:  When the data are not very constraining highly variable sampling can
#       occur near transformation bounds, making it hard to achieve convergence. 
#       In such cases it can be useful to tweak the lower bound (e.g., you might 
#       use pars[,"A"] > 1e-3) to archive convergence.

# NB2: The Ttransform function is also used to combine sampled parameters once 
#      they are  mapped back to the natural scale in order to create parameters 
#      required by the likelihood. For example, in the LBA it has the line
#              pars <- cbind(pars,b=pars[,"B"] + pars[,"A"]),
#      which makes the b, the response threshold, by adding A (the range of the 
#      start-point distribution, to B, the gap from the top of the start-point 
#      distribution to the threshold (in the RDM this is done inside rfun,
#      dfun and pfun).

# NB3: Although we will not change the other elements here are some brief notes
#      1) type: RACE for race models, DDM for the DDM. 
#      2) p_types: vector of default values with the names of the model's 
#         parameter types.
#      3) transform: applied to the parameter vector before mapping
#      4) log_likelihood: the likelihood function to use with the model (for
#         RACE type models usually EMC2:::log_likelihood_race, for DDM type
#         models usually EMC2:::log_likelihood_ddm)
#      5) The dfun and pfun elements just directly call the underlying RDM 
#         density and cumulative density functions, but rfun also calculates the
#         same "ok" test as Ttransform. This is used with data simulation (where
#         inadmissible values can also cause problems). lR is the latent 
#         response factor, which is used to pass the names of responses to the
#         data simulation functions. 

# Finally, as the model is external to EMC2 we must add  "EMC2:::" to functions 
# in the list that are not exported. In this case the functions are rRDM, dRDM, 
# pRDM, and log_likelihood_race.

### The new model ----

# Here is the new model.
RDMt0natural <- function(){
  list(
    type="RACE",
    p_types=c("v" = log(1),"B" = log(1),"A" = log(0),"t0" = 0,"s" = log(1)),
    # Transform to natural scale
    Ntransform=function(x) {
      x[,dimnames(x)[[2]]  != "t0"] <- exp(x[,dimnames(x)[[2]]  != "t0"])
      x
    },
    # p_vector transform
    transform = function(x) x,
    # Trial dependent parameter transform
    Ttransform = function(pars,dadm) {
      attr(pars,"ok") <- (pars[,"t0"] > .05) & ((pars[,"A"] > 1e-6) | pars[,"A"] == 0) &
        ((pars[,"v"] > 1e-3) | pars[,"v"] == 0)
      pars
    },
    # Random function for racing accumulators
    rfun=function(lR=NULL,pars) {
      ok <- (pars[,"t0"] > .05) & ((pars[,"A"] > 1e-6) | pars[,"A"] == 0) &
        ((pars[,"v"] > 1e-3) | pars[,"v"] == 0)
      if (is.null(lR)) ok else EMC2:::rRDM(lR,pars,ok=ok)
    },
    # Density function (PDF) for single accumulator
    dfun=function(rt,pars) EMC2:::dRDM(rt,pars),
    # Probability function (CDF) for single accumulator
    pfun=function(rt,pars) EMC2:::pRDM(rt,pars),
    # Race likelihood combining pfun and dfun
    log_likelihood=function(p_vector,dadm,min_ll=log(1e-10))
      EMC2:::log_likelihood_race(p_vector=p_vector, dadm = dadm, min_ll = min_ll)
  )
}

# Now lets make a design 
designRDMt0natural <- design(data=dat,model=RDMt0natural,
  matchfun=function(d)d$S==d$lR,contrasts=list(lM=ADmat),
  formula=list(B~E+lR,v~lM*CI,A~1,t0~S+E,s~1),constants=c(s=log(1),A=log(0)))

# In this model t0 is no longer on a log scale, so we specify a parameter vector
# correspondingly
p_vectort0natural <- c(B=log(2),B_Espeed=log(1),B_lRright=log(1),v=log(2),
  v_lMd=log(2),v_CIincongruent=log(1),'v_lMd:CIincongruent'=log(1),
  t0=.2,t0_Sright=0.1,t0_Espeed=-0.1)

# We see that additivity on the natural scale now holds, i.e., 
#   left accuracy = t0 = .2
#   right accuracy = t0 + t0_Sright = 0.2 + 0.1 = 0.3
#   left speed = t0 + t0_Espeed = 0.2 - 0.1 = 0.1
#   right speed = t0 + t0_Sright + t0_Espeed = 0.2 + 0.1 - 0.1 = 0.2
mapped_par(p_vectort0natural,designRDMt0natural)

# Lets now check the profiles
sdat <- make_data(p_vectort0natural,designRDMt0natural,n_trials=5000)
profile_plot(sdat, designRDMt0natural, p_vector = p_vectort0natural, 
             layout = c(2,5),n_cores=12)

# We see that profiles and table are consistent with out implementation being
# correct. The step functions for the t0 parameters occur because smaller
# values trigger the "ok" mechanism causing some of the design cells to have
# mapped t0 values that are too small.

# Finally, note that we can reduce the range around the true value that is 
# plotted, giving us a finer grained set of profiles. At this finer resolution
# we more clearly see the small misses caused by sampling noise.
profile_plot(sdat, designRDMt0natural, p_vector = p_vectort0natural, 
             layout = c(2,5), range=.05,n_cores=12)

# Of course with smaller sample sizes such misses are more pronounced
sdat1 <- make_data(p_vectort0natural,designRDMt0natural,n_trials=100)
profile_plot(sdat1, designRDMt0natural, p_vector = p_vectort0natural, 
             layout = c(2,5),range=.05,n_cores=12)

#### Asymptotic parameter recovery study ----

# When you make a new model it is a good idea to do parameter recovery studies 
# to make sure it is identifiable in your design of interest. First we make a
# prior.

p_vectort0natural <- c(B=log(2),B_Espeed=log(1),B_lRright=log(1),v=log(2),
  v_lMd=log(2),v_CIincongruent=log(1),'v_lMd:CIincongruent'=log(1),
  t0=.4,t0_Sright=0.1,t0_Espeed=-0.1)
priorRDMt0 <- prior(designRDMt0natural,type="single",pmean=p_vectort0natural,
                    psd=c(1,.5,.5,1,.5,.5,.5,.1,.05,.05))
plot_prior(priorRDMt0,designRDMt0natural,layout=c(2,6),map=TRUE)

# Here we will just do a large sample single-subject recovery as a check further
# check that the implementation is correct. If you were going to apply a new
# model to your data you may also wish to run the other types of recovery
# studies described in the last lesson.

sRDMt0nL <-  make_emc(sdat,designRDMt0natural,type="single",rt_resolution=.001,
                      prior=priorRDMt0)

# Slow (see .Rout file) so run in batch mode!
# save(sRDMt0nL,file="NewModels/RDMt0nL.RData")

# Check recovery
plot_pars(sRDMt0nL,layout=c(2,5),true_pars=p_vectort0natural,
          use_prior_lim = FALSE)
recovery(sRDMt0nL,true_pars=matrix(p_vectort0natural,nrow=1))


#### Changing Ntransform bounds ----

# Implementing bounds with the "ok" mechanism (e.g., the 0.05 lower bound on t0)
# can be viewed as problematic because they violate the assumption that sampling 
# occurs in an unbound space. Although, as the previous example illustrates, 
# such violations cannot always be avoided for some model parameterizations, in 
# other cases there are more flexible ways of bounding values on the natural 
# scale while maintaining keeping sampling in an unbounded space.

# For example, in the following model puts a lower bound of 0.1s on t0 in the
# Ntranform function. That is, when -Inf is sampled on the natural scale that
# corresponds to 0.1 (i.e., 0.1 + exp(-Inf)). 


RDMt0lb <- function(){
  list(
    type="RACE",
    p_types=c("v" = log(1),"B" = log(1),"A" = log(0),"t0" = log(0),"s" = log(1)),
    # Transform to natural scale
    Ntransform=function(x) {
      x[,dimnames(x)[[2]]  != "t0"] <- exp(x[,dimnames(x)[[2]]  != "t0"])
      x[,dimnames(x)[[2]]  == "t0"] <- 0.1 + exp(x[,dimnames(x)[[2]]  == "t0"])
      x
    },
    # p_vector transform
    transform = function(x) x,
    # Trial dependent parameter transform
    Ttransform = function(pars,dadm) {
      attr(pars,"ok") <- (pars[,"t0"] > .05) & ((pars[,"A"] > 1e-6) | pars[,"A"] == 0)
      pars
    },
    # Random function for racing accumulators
    rfun=function(lR=NULL,pars) {
      ok <- (pars[,"t0"] > .05) & ((pars[,"A"] > 1e-6) | pars[,"A"] == 0)
      if (is.null(lR)) ok else EMC2:::rRDM(lR,pars,ok=ok)
    },
    # Density function (PDF) for single accumulator
    dfun=function(rt,pars) EMC2:::dRDM(rt,pars),
    # Probability function (CDF) for single accumulator
    pfun=function(rt,pars) EMC2:::pRDM(rt,pars),
    # Race likelihood combining pfun and dfun
    log_likelihood=function(p_vector,dadm,min_ll=log(1e-10))
      EMC2:::log_likelihood_race(p_vector=p_vector, dadm = dadm, min_ll = min_ll)
  )
}

# Lets make a model, with a simple t0~1 mapping and check the profile.
designRDMt0lb <- design(data=dat,matchfun=function(d)d$S==d$lR,
  formula=list(B~E+lR,v~lM*CI,A~1,t0~1,s~1),constants=c(s=log(1),A=log(0)),
  contrasts=list(lM=ADmat),model=RDMt0lb)

# Note that when we specify a value for t0 we have to take account of the 
# shfit by 0.1, here so that t0=0.3s
p_vectort0lb <- c(B=log(2),B_Espeed=log(1),B_lRright=log(1),v=log(2),
                  v_lMd=log(2),v_CIincongruent=log(1),'v_lMd:CIincongruent'=log(1),t0=log(.3 - .1))
mapped_par(p_vectort0lb,designRDMt0lb)

# The profiles show the model appears to work. 
sdat <- make_data(p_vectort0lb,designRDMt0lb,n_trials=5000)
profile_plot(sdat, designRDMt0lb, p_vector = p_vectort0lb,layout = c(2,4),n_cores=12)

#### The DDM's complicated scales and bounds ----

# Finally lets take a quick look at the standard DDM model, which has the most
# complicated set of transformation in EMC2.
DDM

# The transformations are of three types 
# 1) natural scale
#   v = (positive favors upper)
# 2) log scale (lower bound only)
#   t0 > 0: lower bound of non-decision time
#   st0 > 0: width of non-decision time distribution
#   a > 0: upper threshold, a
#   sv > 0: standard deviation v
#   s > 0: moment-to-moment standard deviation, s
# probit scale (upper and lower bounds)
#   0 < Z < 1: relative start point, absolute start point, z = Z*a
#   0 < SZ < 1: relative start-point variability, absolute start-point 
#               variability sz = 2*SZ*min(c(a*Z,a*(1-Z))
#   0 < DP < 1: d = t0(upper)-t0(lower) = (2*DP-1)*t0

# The complicated form for SZ ensures that the upper and lower bounds of the 
# start-point range can not fall outside of 0-a. When z approaches 0 or a this
# means that sz shrinks, inducing a dependence so other approaches may 
# sometimes be desirable (e.g, using ok below or an entirely different 
# parameterization). Similarly the transform for DP, which controls the 
# difference in t0 between responses at the upper and lower boundaries (due to 
# different response production times) ensures that t0 for both responses
# remains positive. 

# The DDM likelihood is obtained by numerical integration, which can sometimes
# fail with unusual values. The bounds here were set based on our experience 

ok <- !(abs(pars[, "v"]) > 20 | pars[, "a"] > 10 | pars[, "sv"] > 10 | 
        pars[, "SZ"] > 0.999 | pars[, "st0"] >  0.2 | pars[, "t0"] < 0.05)
if (pars[1, "sv"] != 0) attr(pars, "ok") <- attr(pars,"ok") & pars[, "sv"] > 0.001
if (pars[1, "SZ"] != 0) attr(pars, "ok") <- attr(pars,"ok") & pars[, "SZ"] > 0.001

# These may not always work or be appropriate so may need to be altered (e.g., if
# non-decision time was unusually variable the upper bound on st0 of 0.2 may not
# be appropriate). 
