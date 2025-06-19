rm(list=ls())
library(EMC2)

# This lesson illustrates methods to check whether parameter estimation and
# model selection are valid in a given design and region of parameter space.

# By validity we mean that results can be trusted in a given design conditional 
# on the assumption that the the model that is correct in the sense that it is
# an accurate representation of the data generation process. 

# Given this assumption we can use the model to generative data for the design
# of interest and then fit that simulated data and check that we recover 
# accurate estimates of the parameters used to generate the data or to identify
# the model used to generate the data from models that did not generate the 
# data (either models of the same type with different parameterizations or 
# models of a different type, e.g., a DDM vs. LBA). 

# We will first focus on parameter recovery, where "accuracy" means that 
#  i) The estimated parameters are close to the data generating parameters (i.e., 
#     that any BIAS, that is a systematic difference from the data generating 
#     values is relatively small.  
# ii) Importantly in a Bayesian context, estimates of parameter  uncertainty are 
#     well calibrated, that is, that credible intervals contain true values at a 
#     nominal level (i.e., 95% of the time for a 95% CI), a property that is 
#     called COVERAGE).

# Recovery performance (also called IDENTIFIABILITY) is not a property of a 
# model alone but rather also depends on the design to which it is fit, 
# particularly number of trials and the set of parameters that generate data 
# for that design, but also crucially for EAMs the difficulty of choices and 
# hence the error rate. 

# Hence, recoveries are often done by first fitting some data to determine the
# region of parameter space to focus on (i.e., what data generating parameters
# to use).

# Ideally, recoveries should be based on pilot data, where you first specify a
# potential design, including the number of trials per cell. After fitting that
# data you look at recovery performance and use that as a basis to refine the
# design (e.g., changing the number of trials per cell or perhaps the factorial
# structure of the design).

# Once you are happy that you have a design with sufficiently good recovery 
# properties you can then run the main study.

# NB1: Bayesian parameter estimation is intrinsically biased, in that it is
#      influenced by the prior, although that bias tends to diminish in designs 
#      with more observations.
#
# NB2: We will focus on coverage for 95% credible intervals, as they are the
#      the most commonly used measure of uncertainty. There is a method 
#      available that allow coverage to be addressed for all possible credible
#      intervals (SBC: simulation-based calibration) but that is not yet 
#      implemented in EMC2.
#

# Recovery studies are computationally intensive as they involve simulating 
# and fitting many data sets. To reduce the computational burden we will use
# the LNR model for our illustrations, as it is quick to simulate and calculate
# a likelihood for.

#### Load pre-computed samples, simulated data and results ----

print(load("Recovery/sLNR1.RData"))
print(load("Recovery/sdatsIndividual.RData"))
print(load("Recovery/sRDM.RData"))
load("Recovery/HierarchcialRecoveryData.RData")
load("Recovery/sfDDM.RData")
load("Recovery/sfLBA.RData")
load("Recovery/sfDDMLBA.RData")
load("Recovery/sfLBADDM.RData")

#### Simulating individual data ----

# First load the fits from the Basic lesson.
print(load("BasicEAMs/FlankerSamples.RData"))

# Set up the LNR model
ADmat <- matrix(c(-1/2,1/2),ncol=1,dimnames=list(NULL,"d")); ADmat
designLNR <- design(data=dat,model=LNR,
  matchfun=function(d)d$S==d$lR,
  contrasts=list(lM=-ADmat),
  formula=list(m~lM+CI+E+lM:CI+lM:E,s~lM,t0~1))

pmean <- c(m=-.5,m_lMd=1,m_CIincongruent=0,m_Espeed=0,
  'm_lMd:CIincongruent'=0,'m_lMd:Espeed'=0,s=log(1),s_lMd=log(1),t0=log(.2))
psd <- c(2,2,2,2,2,2,1,1,1)
priorLNR <- prior(designLNR,pmean=pmean,psd=psd,type="single") 

# We can use the "make_data" function to simulate data from a single subject in 
# several different ways (we will illustrate four). 

# 1) The simplest is to supply an emc object, which results in exactly the same 
#    number of trials in each cell as in the data fit with this object and
#    using posterior median parameters.
sdat <- make_data(sLNR)

# We can also explicitly supply a parameter vector, in which case we must also 
# supply the design. In the next examples this is the posterior mean.
pmean <- apply(parameters(sLNR,selection="alpha")[,-1],2,mean)

# 2) We can then simulate data from a balanced design with n_trials trials per 
#    design cell, e.g., 
sdat1 <- make_data(pmean,designLNR,n_trials=100)

# 3) Or by supplying data, in which case the same number of trials per cell as
#    in the data are simulated (this enable simulation of unbalanced designs)
sdat2 <- make_data(pmean,designLNR,data=dat)

# 4) In any of the three previous cases we can expand the trial numbers by a
#    multiple (useful when investigating the number of trials required for 
#    adequate recovery in a given unbalanced design), e.g., 
sdat3 <- make_data(sLNR,expand=2)

#### Parameter recovery for individual data ----

# Now make a samplers object (here we use low resolution as that is another
# feature fo this design).
simLNR <-  make_emc(sdat,designLNR,type="single",rt_resolution=.05,
                      prior=priorLNR)
# Fit
simLNR <- fit(simLNR)

# plot_pars allows allows us to superimpose true values as vertical lines.
# First lets extract the posterior medians
pmedian <- posterior_summary(sLNR)[[1]][,"50%"]

# We see that parameter recovery looks excellent.
plot_pars(simLNR,layout=c(2,5),true_pars=pmedian,
          use_prior_lim=FALSE,plot_prior=FALSE)

# In the case where we generated data from an emc object we can also plot 
# the full distribution of the data generating samples (in green) 
# superimposed on the recovered samples (in black) to provide a fuller picture.
plot_pars(simLNR,layout=c(2,5),true_pars=sLNR,use_prior_lim=FALSE,
          plot_prior=FALSE,true_plot_args = list(cols=c("green")))

# Another way to look at the results it to use the "recovery" function, which 
# for an individual fit plots the data generating and recovered parameters in 
# the same plot with 95% credible intervals (this only works well if all 
# parameters are on roughly the same scale).

# The recovery function also invisibly returns a corresponding table.
print(recovery(simLNR,true_pars=pmedian))

# The stats part here provides pearson and Spearman (ordinal) correlations, both
# of which confirm good recovery and two statistics, the root mean squared error
# (rmse, a measure of absolute disagreement between true and recovered 
# parameters) and the coverage. The latter here is 100%, but given only 9 values
# being tested this is not very meaningful. Note that which of these values are
# provided on the plot can be controlled by the correlation and stat arguments 
# of recovery.

# Similarly to plot_pars we can also provide an emc object to instead plot the
# data-generating samples (with credible intervals) on the x-axis.
print(recovery(simLNR,true_pars=sLNR))

#### Multiple individual parameter recoveries  ----

# A more informative method is to simulate and fit multiple data sets (here
# we will use 100). We begin by creating a matrix of parameters with one row 
# for each data set, where parameters from each data set are randomly drawn from 
# the posterior of the sLNR fit. 

# posterior without subject column sample rows with replacement
p1 <- parameters(sLNR,selection="alpha",resample=TRUE,N=100)[,-1]  

# Next we make a design with 100 subjects. We could write out the whole design
# again, but knowing the internal structure we can also copy and tweak the 
# original design.
designLNR100 <- designLNR
designLNR100$Ffactors$subjects <- 1:100

# The parameter matrices can then be directly passed to make_data, here 
# making a balanced design similar to the original (to make an unbalanced 
# design we would need to make a version of dat with 100 subjects)
sdats1 <- make_data(p1,designLNR100,n_trials=100)

# We will save off these example data sets for later reference
# save(sdats1,file="Recovery/sdatsIndividual.RData")

# Note that when data are simulated a "p_vector" attribute containing the data
# generating parameters is added to the data frame, e.g., 
head(round(attributes(sdats1)$p_vector,2))

# Make and save emc objects
sLNR1 <-  make_emc(sdats1,designLNR100,type="single",rt_resolution=.05,prior=priorLNR)
# save(sLNR1,file="Recovery/sLNR1.RData")

# Fit with runsLNR1.R. Note that in these scripts we use the "subset" command 
# to only keep the adapt and sample stages in order to minimize the
# size of the emc object. Note that although we only use the sample stage
# for analysis it is wise to keep adapt in case you want to add more samples
# (usually the adapt stage is not required but sometimes it can be).

# By default check will only provide a printout for the first subject, but we
# can save a list with results for all subjects
check1 <- check(sLNR1,plot_worst=FALSE)

# These can then be made into table and the worst cases picked out, all showing
# good convergence
gd_summary(sLNR1, selection = "alpha")

# # NB: Make sure the right p1 is being used to check recovery by getting it 
# # from the saved sdats1
# load("Recovery/sdatsIndividual.RData")
# p1 <- attributes(sdats1)$p_vector


# We first plot the true parameters and posterior median estimates along with
# correlations between them (the default) and "coverage", the percentage of 
# times that the true value falls within the 95% CI of the estimated value.
recovery(sLNR1,true_pars=t(p1),do_CI=FALSE,stat="coverage",selection="alpha")

# In some cases recovery is relatively good (e.g., s_lMd) but in others strong
# "shrinkage" can be seen (i.e., the variation of the estimates is much less 
# that the true variation, so the scatter plot is flattened and the estimates
# "pulled" (biased) towards the prior mode (e.g., m_CIincongruent). 

# The strong shrinkage is likely due to the low error rate in the data from 
# which these parameters were derived, and salutary of the need to design
# experiments in a way that produces higher error rates (ideally 10% or more)
# when fitting EAMs.

# However, it is also reassuring that the coverage is close to nominal in most
# cases, meaning that inferences about parameter values are appropriately 
# calibrated (although it is important to acknowledge that with only 100 
# replications the estimates of coverage are fairly crude). To see coverage
# graphically we can also plot the 95% CIs. 
recov <- recovery(sLNR1,true_pars=p1,correlation="spearman",selection="alpha")

# Here we also displayed two other summary statistics, the spearman (ordinal)
# correlation and root-mean squared error and save results. 

# The results are a list named for each parameter of summary statistics
head(recov$t0$stats)

# and quantiles
head(recov$t0$quantiles)

# The latter can be useful for calculating other statistics, such as the average
# bias, here demonstrating the under-estimation of t0 evident in the plot
mean(recov$t0$quantiles[,"miss"])

#### Individual model recovery ----

# In a model recovery study we fit data generated by one model with a different
# model. We then look at the ability of a model selection method to choose the
# correct data generating model.

# We first fit the RDM model from the Basic lesson to the simulated LNR data

designRDM <- design(data=dat,model=RDM,
                  formula=list(B~E+lR,v~lM*CI,t0~1),
                  matchfun=function(d)d$S==d$lR,
                  contrasts=list(lM=ADmat))

pmean <- c(B=log(2),B_Espeed=log(1),B_lRright=log(1),v=log(2),
  v_lMd=log(2),v_CIincongruent=log(1),'v_lMd:CIincongruent'=log(1),t0=log(.2))
psd <-  c(1,.5,.5,1,.5,.5,.5,.5)
priorRDM <- prior(designRDM,pmean=pmean,psd=psd,type="single")

sRDM <- make_emc(sdats1,designRDM,type="single",rt_resolution=.05,
                      prior=priorRDM)
# save(sRDM,file="Recovery/sRDM.RData")

# We then apply the "compare_subject" function, which calculates model selection 
# weights based on DIC and BPIC. 

LNR1RDM <- compare_subject(list(Correct=sLNR1,Wrong=sRDM))

# While this is running we note that as well as printing a table of model
# weights on completion (which shows that the correct model as selected in all
# cases by both DIC and BPIC) we can also save the output, which contains 
# detailed results.

# For example, for the first subject we see that the RDM is simpler (smaller
# EffectiveN) but fits much better (bigger Dmean/minD) so overall it is not
# selected. Thes same is true for the other simulated subjects.
LNR1RDM[[1]]

#### Hierarchical parameter recovery  ----

# We will check how well the parameters of the selected LBA model from the 
# 4-Hierarchical can be recovered in the forstmann design by simulating data based on 
# posterior median estimates and fitting the simulated data.

# Get data, design, prior, and samples for the original model fit.
print(load("Hierarchical/fits_DDM/DDM.RData"))
print(load("Hierarchical/fits_LBA/B.RData"))

Emat <- matrix(c(0,-1,0,-1,0,0),nrow=3)
Vmat <- matrix(c(-1,1),ncol=1,dimnames=list(NULL,"r"))
design_DDM <- design(data=forstmann,model=DDM,
  contrasts=list(S=Vmat,E=Emat),
  formula=list(v~S,a~E, Z~1, t0~1,sv~1,SZ~1)
)

pmean <- c(v=0,v_Sr=2.25,a=0.2,'a_Ea-n'=0,'a_Ea-s'=0,t0=-0.8,Z=0,sv=0.2,SZ=-0.3)
psd <- c(2,2.5,0.75,.25,.25,0.4,0.4,.5,0.75)
priorDDM <- prior(design_DDM,mu_mean=pmean,mu_sd=psd)

design_LBA <- design(matchfun=function(d)d$S==d$lR,
                   data=forstmann,model=LBA,
                   contrasts=list(lM=ADmat,E=Emat),
                   formula=list(v~lM,B~E,A~1,t0~1))

pmean <-  c(v=3,v_lMd=2,B=log(2),'B_Ea-n'=0,'B_Ea-s'=0,A=log(1),t0=log(.2))
psd <- c(2,2,1,.25,.25,0.5,.5)
priorLBA <- prior(design_LBA,mu_mean = pmean, mu_sd = psd)


# Simulate data from the same design as fit by these samples objects
fDDM <- make_data(sDDM)
fLBA <- make_data(sB)

# save(fDDM,fLBA,file="Recovery/HierarchcialRecoveryData.RData")

# When make_data is supplied with hierarchical samples it simulates based on
# posterior median alphas for each participant.
attr(fDDM,"p_vector")[1,]
posterior_summary(sDDM,selection="alpha",subject=1)[[1]][,"50%"]


sfDDM <- make_emc(fDDM,design_DDM,prior=priorDDM)
# save(sfDDM,file="Recovery/sfDDM.RData")

sfLBA <- make_emc(fLBA,design_LBA,prior=priorLBA)
# save(sfLBA,file="Recovery/sfLBA.RData")


pDDM <- attr(fDDM,"p_vector")
pLBA <- attr(fLBA,"p_vector")

# DDM is good for all but between-trial variability parameters
# which show strong shrinkage to the prior. The number of subjects (19) is 
# low for an accurate estimate of coverage, but seems fairly good.
recovery(sfDDM,true_pars=pDDM,selection="alpha",stat="coverage")

# LBA is good in most cases, with some shrinkage of high values of B_E2 evident
# and generally OK coverage.
recovery(sfLBA,true_pars=pLBA,selection="alpha",stat="coverage")

# Note that we can also plot the posterior uncertainty for both the data
# generating and recovered samples simultaneously. Note that we would expect
# would shrinkage in the recovered. 
recovery(sfDDM,true_pars=sDDM,do_CI=T,stat="coverage",selection="mu")
recovery(sfDDM,true_pars=sDDM,do_CI=T,stat="coverage",selection="sigma2")

# We can also plot the densities on top of each other
plot_pars(sfDDM, true_pars = sDDM, use_prior_lim = F)

#### Hierarchical Model Selection ----

# Now lets "cross-fit" the LBA to the DDM data and see how well model selection
# can pick the data generating model

sfDDMLBA <- make_emc(fDDM,design_LBA,prior=priorLBA)
# save(sfDDMLBA,file="Recovery/sfDDMLBA.RData")

# and the DDM to the LBA data.

sfLBADDM <- make_emc(fLBA,design_DDM,prior=priorDDM)
# save(sfLBADDM,file="Recovery/sfLBADDM.RData")

# All metrics strongly prefer the correct DDM model, largely because the correct
# model fits a lot better.
compare(list(Correct=sfDDM,Wrong=sfDDMLBA),cores_per_prop = 3)
#            MD wMD    DIC wDIC   BPIC wBPIC EffectiveN  meanD  Dmean   minD
# Correct -8835   1 -10591    1 -10373     1        218 -10809 -10936 -11027
# Wrong   -7114   0  -8136    0  -7948     0        188  -8324  -8437  -8513

# And all metrics strongly prefer the correct LBA model, which again fits 
# better.
compare(list(Correct=sfLBA,Wrong=sfLBADDM),cores_per_prop = 3)
#             MD wMD    DIC wDIC   BPIC wBPIC EffectiveN  meanD  Dmean   minD
# Correct -13273   1 -14377    1 -14202     1        175 -14552 -14679 -14727
# Wrong    -8178   0 -10027    0  -9790     0        237 -10264 -10376 -10500
