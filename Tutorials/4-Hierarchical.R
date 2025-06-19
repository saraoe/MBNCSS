rm(list=ls())
library(EMC2)

# In this lesson you will learn how to build "hierarchical" models, which
# simultaneously estimate parameters of an EAMs for each member of a group of 
# individuals and also a estimates a model of the population distribution of
# EAM parameters from which those individuals are drawn.

# You will likely mostly use hierarchical modeling in most real applications,
# and so it basic form, which assumes a multivariate normal population 
# distribution, is EMC2's "standard" type when making emc objects. 

# We will illustrate hierarchical modeling using the forstman data set that 
# comes with EMC2 (see ?forstmann and the next section). We will fit and 
# compare hierarchical DDM and LBA models (extension to the other basic EMC2
# EAMs, the RDM and LNR, are left as an exercise).

# Hierarchical sampling can be slow, so to work through this lesson without 
# having to run the sampling for all models load pre-computed samples and
# posterior predictives.
#
print(load("Hierarchical/fits_DDM/DDM.RData"))
print(load("Hierarchical/fits_DDM/DDMt0.RData"))
print(load("Hierarchical/fits_DDM/DDMv.RData"))
print(load("Hierarchical/fits_DDM/DDMvt0.RData"))
print(load("Hierarchical/fits_LBA/B.RData"))
print(load("Hierarchical/fits_LBA/BBa.RData"))
print(load("Hierarchical/fits_LBA/BBasv.RData"))
print(load("Hierarchical/fits_LBA/Eat0.RData"))
print(load("Hierarchical/fits_LBA/Eav.RData"))
print(load("Hierarchical/fits_LBA/Eavt0.RData"))
print(load("Hierarchical/fits_LBA/Eavt02.RData"))
print(load("Hierarchical/fits_DDM/ppDDM.RData"))
print(load("Hierarchical/fits_LBA/ppLBA.RData"))
print(load("Hierarchical/extras.RData"))


#### The "forstmann" data  ----

# 2 x 3 design

# The 2 level S factor is for left vs. right movement in a 
# random-dot display. 

# The three level factor is for instructions emphasizing speed, equally 
# emphasizing speed and accuracy (neutral) or emphasizing accuracy.
lapply(forstmann,levels)

# There are 19 participants, ~800 trials each (some less)
table(forstmann$subjects)

# Take a look at the first subject (here increasing the default smoothing on the
# density plot with adjust)
plot_defective_density(forstmann,layout=c(2,3),subject=1,adjust=1.5)

# The default is to plot results aggregated over subjects
plot_defective_density(forstmann,layout=c(2,3))

# Not much difference between accuracy and neutral in error rate, neutral a
# a little slower, errors the same or a little slower than correct. Speed less
# accurate and substantially faster, with errors faster than corrects, perhaps
# some bias to respond right (more accurate, slower).

#### DDM ----

# We will first analyze the forstmann data with the full DDM, except st0 is 
# fixed at zero (which greatly speeds sampling, see Basic Exercises)

### Contrasts  ----

# Set up contrast matrices to compare caution in neutral and speed to accuracy.

# Accuracy - Neutral, Accuracy - Speed (note the -1, if we used 1
# it would be e.g., Neutral - Speed, but we expect accuracy to be
# bigger so to keep the contrast positive used -1)
Emat <- matrix(c(0,-1,0,-1,0,0),nrow=3)
dimnames(Emat) <- list(NULL,c("a-n","a-s"))
Emat


# For rates, we use the intercept and rate to upper boundary contrast matrix
Vmat <- matrix(c(-1,1),ncol=1,dimnames=list(NULL,"r"))
Vmat

### A conventional DDM model ----

# We make the typical assumption that instructions related to speed and 
# accuracy (i.e., the E factor) affect only thresholds. We also allow rates 
# (both rate bias and sensitivity) to differ for left and right stimuli. 
design_DDM <- design(data=forstmann,model=DDM,
  contrasts=list(S=Vmat,E=Emat),
  formula=list(v~S,a~E, Z~1, t0~1,sv~1,SZ~1)
)

# This results in a 9 parameter model
sampled_p_vector(design_DDM,doMap=FALSE)

## Prior ----

# We use a (somewhat) informed prior for the population mean and standard 
# deviation. This is just a slightly rounded version of the prior used in the 
# Basic lesson for the DDM.
pmean <- c(v=0,v_Sr=2.25,a=0.2,'a_Ea-n'=0,'a_Ea-s'=0,t0=-0.8,Z=0,sv=0.2,SZ=-0.3)
psd <- c(2,2.5,0.75,.25,.25,0.4,0.4,.5,0.75)

# When making the prior we do not have to specify type as the default is 
# "standard"
priorDDM <- prior(design_DDM,mu_mean=pmean,mu_sd=psd)

# My default we get a plot of the prior for the population means.
plot_prior(priorDDM,design=design_DDM)

# The standard model also requires a prior for the population variance (sigma2).
# As we did not specify it in prior we are using the default form, a half-t 
# distribution (df) with 2 degrees of freedom.
plot_prior(priorDDM,design=design_DDM,selection="sigma2")

# The default has an uninformative (uniform between -1 and 1) distribution on
# correlations. 
plot_prior(priorDDM,design=design_DDM,selection="correlation",
           flatten=TRUE,layout=c(4,7))

# Which corresponds to a distribution on the covariance as a gamma mixture of 
# inverse Wishart distributions.
plot_prior(priorDDM,design=design_DDM,selection="covariance",
           flatten=TRUE,layout=c(4,7))

# The A (scale, default = .3) and v (degrees of freedom, default = 2) parameters 
# control the shape of the population variance-covariance prior. For example,
# increasing the scale produces a wider distribution.
priorDDM1 <- prior(design_DDM,mu_mean=pmean,mu_sd=psd,A=1)
plot_prior(priorDDM1,design=design_DDM,selection="sigma2",xlim=c(0,5))

## Sampling ----

# First make an emc object in the same way as before, binding together the data,
# design and prior, but as type="standard" is the default we dont need to 
# specify that. 
samplers <- make_emc(forstmann,design_DDM,prior=priorDDM)

# Because sampling takes a while to run we typically do not fit a hierarchical 
# model by running a command in the console. Instead we fit in "batch" mode.

# The first step is to save the emc object to disk. It is often useful to 
# do this (do this in a sub-directory, here we use "Hierarchical/fits_DDM" to
# contain all DDM model fits). 

# # NB: Don't run this line as it will overwrite the pre-computed samples!)
# save(samplers,file="Hierarchical/fits_DDM/DDM.RData")

# We then make a text file file (usually with a .R extension, here it was called
# runDDM.R) in the same directory as the file containing the emc object. This
# file which contains the commands required to run fitting, as follows"

# library(EMC2)
# load("DDM.RData")
# sDDM <- fit(samplers,verbose=TRUE,cores_per_chain = 4,fileName = "tmpDDM.RData")
# save(sDDM,file="DDM.RData")

# The "fileName" argument causes intermediate results to be saved to a file 
# after every fitting cycle, so that the results of a long run are not lost by
# some mishap.

# The "cores_per_chain" argument, which only works on Linux and Mac (it will be
# ignored on Windows), speeds fitting by using multiple cores for each chain.
# With the default 3 chains that means the setting above 12 cores will be used
# in total. 

# NB: Note that increasing cores_for_chains beyond 1 is not as efficient as 
#     allocating a core to each chain, at least for models with likelihoods
#     that are easy to compute.

# The batch file can be run by opening a console (Tools >> Terminal >> New
# Terminal in RStudio), navigating to the directory containing the files (the
# terminal opens in the project directory so type "cd Hierarchical/fits_DDM/")
# and typing the following (you can also use an Rscript version)
#      nohup R CMD BATCH runDDM.R &

# "nohup" mean "no hangups", it means fitting continues if the terminal closes
# It is not necessary if you will keep the terminal open. If you do use it a
# "nohup.out" file will be created, but its contents are not of interest.

# The "&" at the end of the command line causes the fit to run in background,
# so that after you hit enter (sometimes you will need to do it twice) the 
# command prompt reappears and you can enter further commands.

# When running in batch mode output that would normally be printed in the 
# console is saved in a ".Rout" file, in this case "runDDM.Rout". You can open 
# this file while fitting is ongoing to monitor progress. 

### Alternative DDM models ----

# In this lesson we will fit four DDM models, the standard one just specified 
# and three others that add to the conventional assumption that the E
# manipulation selectively affects thresholds additional flexibility, the first
# two by allowing E to affect t0 or v, and then both in the third.

# Letting E affect t0 produces an 11 parameter model. Note that t0 will use
# Emat just like a.
design_DDMt0 <- design(data=forstmann, model=DDM,
  contrasts=list(S=Vmat,E=Emat),
  formula=list(v~S,a~E, t0~E, Z~1, sv~1, SZ~1)
)

# Add in our usual vague effect prior and zero effect means.
pmean <- c(v=0,v_Sr=2.25,a=0.2,'a_Ea-n'=0,'a_Ea-s'=0,
           t0=-0.8,'t0_Ea-n'=0,'t0_Ea-s'=0,Z=0,sv=0.2,SZ=-0.3)
psd <- c(2,2.5,0.75,.25,.25,0.4,1,1,0.4,.5,0.75)
priorDDMt0 <- prior(design_DDMt0,mu_mean = pmean, mu_sd = psd)
plot_prior(priorDDMt0,design_DDMt0,layout=c(3,4))

samplers <- make_emc(forstmann,design_DDMt0,prior=priorDDMt0)
# save(samplers,file="Hierarchical/fits_DDM/DDMt0.RData")

# In the next model we illustrate how to apply a different contrast for the same
# factor when it affects different parameter types, keeping Emat for a  and 
# E2mat (defined below) for v. To archive this the assignments in contrasts list 
# are done in lists named for the parameters they apply to. 

# The new contrast simply reverses the signs in Emat. E2mat still uses the 
# accuracy condition as the baseline (so the v parameter corresponds to rate 
# quantity for the accuracy condition) and sets up "neutral" and "speed" 
# contrasts, which are the change from the accuracy baseline to these conditions.
E2mat <- cbind(neutral=c(0,1,0),speed=c(1,0,0)); E2mat

# When when v is affected by both S and E we have a 13 parameter model (the
# extra two parameters relative to the last model occur because we allow for
# both an rate bias (main effect) of E (v_Eneutral, v_Espeed) and different
# rate quality (Sr) effects for the between emphasis conditions (v_Eneutral:Sr,
# v_Eaccuracy:Sr and v_Espeed:Sr).
design_DDMv <- design(data=forstmann,model=DDM,
  contrasts=list(v=list(S=Vmat,E=E2mat),a=list(E=Emat)),
  formula=list(v~E/S,a~E, t0~1, Z~1, sv~1, SZ~1)
)

# For our prior we set effects to zero except for r (rate), which we expect to
# be positive, reflecting fairly accurate performance (using the same value as
# our previous prior). 
pmean <- c(v=0,v_Eneutral=0,v_Espeed=0,'v_Eaccuracy:Sr'=2.25,
           'v_Eneutral:Sr'=2.25,'v_Espeed:Sr'=2.25,a=0.2,'a_Ea-n'=0,'a_Ea-s'=0,
           t0=-0.8,Z=0,sv=0.2,SZ=-0.3)
psd <- c(2,2,2,2.5,2.5,2.5,0.75,.25,.25,0.4,0.4,.5,0.75)
priorDDMv <- prior(design_DDMv,mu_mean = pmean, mu_sd = psd)
plot_prior(priorDDMv,design_DDMv,layout=c(3,5))

samplers <- make_emc(forstmann,design_DDMv,prior=priorDDMv)
# save(samplers,file="Hierarchical/fits_DDM/DDMv.RData")


# Finally, we specify a 15 parameter model that adds E effects to both v and t0.
design_DDMvt0 <- design(data=forstmann,model=DDM,
  contrasts=list(v=list(S=Vmat,E=E2mat),a=list(E=Emat),t0=list(E=Emat)),
  formula=list(v~E/S,a~E, t0~E, Z~1, sv~1, SZ~1)
)

pmean <- c(v=0,v_Eneutral=0,v_Espeed=0,
           'v_Eaccuracy:Sr'=2.25,'v_Eneutral:Sr'=2.25,'v_Espeed:Sr'=2.25,
           a=0.2,'a_Ea-n'=0,'a_Ea-s'=0,
           t0=-0.8,'t0_Eneutral'=0,'t0_Espeed'=0,Z=0,sv=0.2,SZ=-0.3)
psd <- c(2,2,2,2.5,2.5,2.5,0.75,.25,.25,0.4,1,1,0.4,.5,0.75)
priorDDMvt0 <- prior(design_DDMvt0, mu_mean = pmean, mu_sd = psd)
plot_prior(priorDDMvt0,design_DDMvt0,layout=c(3,5))

samplers <- make_emc(forstmann,design_DDMvt0,prior=priorDDMvt0)
# save(samplers,file="Hierarchical/fits_DDM/DDMvt0.RData")

### Checking fits ----

print(load("Hierarchical/fits_DDM/DDM.RData"))
# We now have new checks for the population mean (mu), variance (sigma2) and 
# correlations as well as the individual estimates (alpha). By default the check
# function does not do correlations 

check(sDDM)
# The first model is well converged, SZ and sv inefficient, which is usually the 
# case, with other parameters between 1/3 and 1/2 efficiency (recall nominally 
# there are 3000 samples).

# We can also look at correlations, which are in line with the previous results
check(sDDM,selection = "correlation")

# Using the selection argument we can look at any subset, shown using the 
# remaining two possible values for selection (redundantly as Sigma contains the 
# covariance).
check(sDDM,selection = c("covariance","Sigma"))

# Lets now look at convergence of the population mean and variance, and 
# individual parameters for the other models (we leave looking at the 
# correlation convergence as an exercise), which are similarly well converged.
print(load("Hierarchical/fits_DDM/DDMt0.RData"))
check(sDDMt0)
print(load("Hierarchical/fits_DDM/DDMv.RData"))
check(sDDMv)
print(load("Hierarchical/fits_DDM/DDMvt0.RData"))
check(sDDMvt0)

### Model comparison ----

# A comparison among the four models tests the conventional assumption that
# speed emphasis affects only thresholds.

# Computing model comparisons for hierarchical models can be slow so we use 
# more cores (remembering that the default here is 4 chains, so the following 
# uses 12 cores in total). 

# IC supports the most complex model, but MD supports the DDMt0.
compare(list(DDM=sDDM,DDMt0=sDDMt0,DDMv=sDDMv,DDMvt0=sDDMvt0),cores_per_prop = 3)
#            MD wMD    DIC wDIC   BPIC wBPIC EffectiveN  meanD  Dmean   minD
# DDM    -11439   0 -13024    0 -12802     0        222 -13247 -13404 -13469
# DDMt0  -12191   1 -14434    0 -14148     1        286 -14721 -14926 -15007
# DDMv   -11106   0 -13152    0 -12787     0        364 -13516 -13768 -13880
# DDMvt0 -11650   0 -14482    1 -14053     0        428 -14910 -15223 -15338

# The same comparisons can be made on a per-subject basis with a summary
# table of winning models. This function does not yet support MD or multicore 
# so can be a bit slow. 

# As at the group level, the most complex model wins in all cases.  
compare_subject(list(DDM=sDDM,DDMt0=sDDMt0,DDMv=sDDMv,DDMvt0=sDDMvt0))
# ...
# Winners
#      DDMvt0
# DIC      19
# BPIC     19

### Priors and posteriors updating ----

# Lets look at the most complex model's parameters.

# First we look at population means (mu, default selection for the following 
# functions): priors strongly dominated, even sv and SZ
plot_pars(sDDMvt0,layout=c(3,5))

# Can also see on the natural scale, here using the posterior scale to focus on
# estimated values rather than dominance. Parameters are fairly well localized,
# and this can also be seen in tabular form.
plot_pars(sDDMvt0,layout=c(3,5),map=TRUE,use_prior_lim=FALSE)
posterior_summary(sDDMvt0,map=TRUE)

# Thresholds conform the expectation of a lower value for speed, with neutral 
# and accuracy quite similar. 

# Rates have a smaller magnitude for speed, with neutral and accuracy again 
# fairly similar. 

# Non-decision time follows the same pattern, less for speed and similar for
# neutral and accuracy.


# We can also look at variances. Estimates tend to be less well localized
# reflecting the fact that individual differences are not well estimated with
# only 19 participants.
plot_pars(sDDMvt0,layout=c(3,5),use_prior_lim=FALSE,selection="sigma2")
posterior_summary(sDDMvt0,selection="sigma2")

# NB: The contraction statistics is not very useful with variances as it is 
#     heavily influenced by tail values. In future we plan to investigate a 
#     robust version of this statistic.

# Note that map is not implemented for variances and correlations/covariances
# as they are most relevant to the unbounded scales on which estimation was
# performed. 

# For ease of reading correlation plots are grouped by the first the first part
# of each correlation (so plots are redundant, 9*8/2 =  36 unique).
plot_pars(sDDMvt0,layout=c(3,5),use_prior_lim=FALSE,selection="correlation")

# We can use flatten to look at them without redundancy
posterior_summary(sDDMvt0,selection="correlation",flatten=TRUE)

# These functions also have Sigma and covariance selections.

# We can look at individual participant posteriors, which are grouped
# by parameters, so each panel shows results for all subjects together.
plot_pars(sDDMvt0,layout=c(4,5),selection="alpha")
posterior_summary(sDDMvt0,selection="alpha")
# Typically these are not as dominated as population means, but remember priors 
# are now very informative, as they are based on population distributions) and 
# so are more influential, particularly evident with sv and SZ.

# We can plot the individual parameters mapped.
plot_pars(sDDMvt0,selection="alpha",layout=c(4,5),map=TRUE)
posterior_summary(sDDMvt0,selection="alpha",map=TRUE)

# We can plot each participant's density overlaid on the hierarchical 
# prior (i.e., the prior implied by the estimated population distribution). This 
# reveals that the participant sv posteriors are strongly pulled to the prior.
plot_pars(sDDMvt0,layout=c(3,5),all_subjects = TRUE, map=TRUE)

# Finally, we could select an individual subject, either by name or with an
# integer, in which case all parameters are grouped on one page. Again the
# prior implied by the population model is shown in red.
plot_pars(sDDMvt0,layout=c(3,5),selection="alpha",subject="as1t")

# NB1: If the subject has more than one entry plotting reverts to grouping by 
#      parameter.

# Parameters from all chains can be combined into a single data frame for 
# further processing, e.g., 
head(round(parameters(sDDMvt0,selection="mu"),3))

# This function also allows you to get a subset by random sampling, either
# without resampling (in which case N must be less than the total number of
# samples)
parameters(sDDMvt0,selection="sigma2",N=2,resample=FALSE)

# or with resampling, in which case it can be more
dim(parameters(sDDMvt0,selection="correlation",N=4000,resample=TRUE))

# For alphas a column indexing subjects is added
head(parameters(sDDMvt0,selection="alpha"))

### Goodness-of-fit ----

# To examine goodness we first calculate posterior predictives

# # We use multi-core to speed computation (note that in this case there is no 
# # reporting of progress).
# ppDDM <- predict(sDDMvt0,n_cores=12)
# save(ppDDM,file="Hierarchical/fits_DDM/ppDDM.RData")

# We can focus in on one subject Fit for subject 1. By default fits are
# separated into the full set of design cells.
plot_fit(forstmann,ppDDM,layout=c(2,3),subject=1)

# But most commonly we look at the the average. Both data and fits are
# are aggregated in the same way (i.e.,treated as a single subject), which means
# that any averaging artifacts apply equally to both.

# We make the plot over the truncation range and again the default is to show 
# the full breakdown by design cells. Although accuracy is fairly well
# captured there is clear misfit to RT distribution (i.e., the grey dots do not
# align very well with the black dots). 
plot_fit(forstmann,ppDDM,layout=c(2,3),lpos="right",xlim=c(.25,1.5))

# Note that posterior predictive uncertainty is low because results are 
# aggregated over the > 2500 data points in each cell.

# Lets focus in on key statistics to evaluate fit in more details. We use an 
# invisible return to save the results as a table as well.

# Accuracy is over-predicted and outside the 95% CI in half of the cells.
pc <- function(d) 100*mean(d$S==d$R)
tab <- plot_fit(forstmann,ppDDM,factors=c("E","S"),layout=c(2,3),
                stat=pc,stat_name="Accuracy (%)",xlim=c(70,95))
round(tab,2)

# Mean RT for correct responses is also over-predicted, again outside the 95% CI 
# in half of the cells.
tab <- plot_fit(forstmann,ppDDM,factors=c("E","S"),layout=c(2,3),
  stat=function(d){mean(d$rt[d$R==d$S])},stat_name="Mean Correct RT (s)",xlim=c(0.375,.6))
round(tab,3)

# Fast correct responses (10th percentile) are strongly under-predicted
tab <- plot_fit(forstmann,ppDDM,xlim=c(0.275,.4),factors=c("E","S"),layout=c(2,3),
  stat=function(d){quantile(d$rt[d$R==d$S],.1)},stat_name="10th Percentile Correct (s)")
round(tab,2)

# Slow correct responses (90th percentile) are strongly over-predicted
tab <- plot_fit(forstmann,ppDDM,xlim=c(0.525,.975),factors=c("E","S"),layout=c(2,3),
  stat=function(d){quantile(d$rt[d$R==d$S],.9)},stat_name="90th Correct Percentile (s)")
round(tab,2)

# Hence, correct RT variability (SD) is also strongly over-predicted.
tab <- plot_fit(forstmann,ppDDM,factors=c("E","S"),layout=c(2,3),
  stat=function(d){sd(d$rt[d$R==d$S])},stat_name="SD Correct (s)",xlim=c(0.1,.325))
round(tab,3)

# As there are fewer errors we look at only mean error RT: like mean correct RT
# it is is over-predicted.
tab <- plot_fit(forstmann,ppDDM,xlim=c(0.375,.725),factors=c("E","S"),layout=c(2,3),
  stat=function(d){mean(d$rt[d$R!=d$S])},stat_name="Mean Error RT (s)")
round(tab,2)

# The DDM does not seen to be an adequate model, so we will not further examine
# its parameters but instead test if the LBA can do better.

#### LBA Models ----

# We will now create a set of 6 models to illustrate various aspects of the LBA.
# We then load in the fits (stored in hierarchical_fits along with scripts to
# run each), check them, and compare to see which gives the best account of the
# data.

### Revising race model contrasts ----

# To reiterate points made in the Basic lesson with individual fitting, for 
# race models the Rlevels argument automatically makes available a "latent
# response" (lR) factor that indexes each accumulator by its corresponding
# response. We discuss below how lR can be used to allow for response bias
# as race models  do not have a separate response-bias parameter as does the
# DDM (i.e., Z). Note that unlike the DDM the LBA and other such race models
# can accommodate more than two responses by having more than two Rlevels.

# For race models like the LBA a new argument to design, matchfun, which is
# supplied with a function that indicates which accumulator (indexed by the lR
# factor) corresponds to the a correct response got different stimuli,
# automatically making available a "latent match" factor (lM)
# with levels that can be used to set different rates for matching and
# mismatching accumulators.

# First lets create some contrast matrices to capture rate and threshold effects
# in theoretically interesting ways.

# RATES
# Average rate = intercept, and rate d = difference (match-mismatch) contrast,
# The intercept is of interest as an index of decision urgency and stimulus
# magnitude effects. d indexes how easy the choices are to discriminate and also
# how effective a subject's selective attention is in making this discrimination.
ADmat <- matrix(c(-1/2,1/2),ncol=1,dimnames=list(NULL,"d")); ADmat

# THRESHOLDS
# The standard assumption is that Emphasis affects thresholds. As with the DDM
# we use the following contrasts to compare accuracy - neutral (a-n) and 
# accuracy - speed (a-s)
Emat <- matrix(c(0,-1,0,-1,0,0),nrow=3)
dimnames(Emat) <- list(NULL,c("a-n","a-s")); Emat


### Old-style conventional model ----

# We make the conventional assumption that only B affected by E. To fix scale 
# we fix sv=1 (by default), as assumption commonly made in early applications
# of the LBA
design_B <- design(matchfun=function(d)d$S==d$lR,
                   data=forstmann,model=LBA,
                   contrasts=list(lM=ADmat,E=Emat),
                   formula=list(v~lM,B~E,A~1,t0~1))


# We put slightly informed priors on the parameters, again similar to the
# prior used in the Basic lesson.
pmean <-  c(v=3,v_lMd=2,B=log(2),'B_Ea-n'=0,'B_Ea-s'=0,A=log(1),t0=log(.2))
psd <- c(2,2,1,.25,.25,0.5,.5)
priorB <- prior(design_B,mu_mean = pmean, mu_sd = psd)
plot_prior(priorB,design_B,layout=c(2,4))

samplers <- make_emc(forstmann,design_B,prior=priorB)
# # Run previously so don't run it unless you want to overwrite the saved samples.
# save(samplers,file="Hierarchical/fits_LBA/B.RData")

# NB: All LBA models were run by a single (runLBA.R) script.

### Adding response bias ----

# Analogous to the DDM, we used an additive lR effect on the threshold to 
# accommodate an overall response bias that does not vary over emphasis 
# conditions (B_lRright), just as in the DDM where we estimated only on response
# bias (Z) value. All LBA models below make this assumption as well.
design_BBa <- design(data=forstmann,model=LBA,
                    matchfun=function(d)d$S==d$lR,
                     contrasts=list(lM=ADmat,E=Emat),
                     formula=list(v~lM,B~E+lR,A~1,t0~1))

# Again we use a population prior similar to the individual prior used in the
# Basic lesson.
pmean <- c(v=3,v_lMd=2,B=log(2),'B_Ea-n'=0,'B_Ea-s'=0,
        'B_lRright'=0,A=log(1),t0=log(.2))
psd <- c(2,2,1,.25,.25,.25,0.5,.5)
priorBBa <- prior(design_BBa,mu_mean = pmean, mu_sd = psd)
plot_prior(priorBBa,design_BBa,layout=c(2,5))

samplers <- make_emc(forstmann,design_BBa,prior=priorBBa)
# save(samplers,file="Hierarchical/fits_LBA/BBa.RData")

### Conventional new style B emphasis effect ----

# More recent applications of the LBA have allowed for sv differing between 
# matching and mismatching accumulators. We allow this extra flexibility in this
# "baseline" model, fixing the intercept sv = 1 to set the scale. 

# All following models will make this assumption.

# Note that as sv is on a log scale the the sv_lMd contrast is multiplicative. 
design_BBasv <- design(data=forstmann,model=LBA,
  matchfun=function(d)d$S==d$lR,
  contrasts=list(lM=ADmat,E=Emat),
  formula=list(v~lM,B~E+lR,A~1,t0~1,sv~lM),
  constants=c(sv=log(1)))

pmean = c(v=3,v_lMd=2,B=log(2),'B_Ea-n'=0,'B_Ea-s'=0,
          'B_lRright'=0,A=log(1),t0=log(.2),sv_lMd=0)
psd = c(2,2,1,.25,.25,.25,0.5,.5,2)
priorBBasv <- prior(design_BBasv,mu_mean = pmean, mu_sd = psd)
plot_prior(priorBBasv,design_BBasv,layout=c(2,6))

samplers <- make_emc(forstmann,design_BBasv,prior=priorBBasv)
# save(samplers,file="Hierarchical/fits_LBA/BBasv.RData")
# Lets call the last one the baseline model and add to it ...


###  B + t0 emphasis effects ----

# As for the DDM we use the Emat contrast for t0.
design_Eat0 <- design(data=forstmann,model=LBA,
  matchfun=function(d)d$S==d$lR,
  contrasts=list(lM=ADmat,E=Emat),
  formula=list(v~lM,B~E+lR,A~1,t0~E,sv~lM),
  constants=c(sv=log(1)))

pmean = c(v=3,v_lMd=2,B=log(2),'B_Ea-n'=0,'B_Ea-s'=0,
          'B_lRright'=0,A=log(1),t0=log(.2),'t0_Ea-n'=0,'t0_Ea-s'=0,sv_lMd=0)
psd = c(2,2,1,.25,.25,.25,0.5,.5,1,1,1)
priorEat0 <- prior(design_Eat0,mu_mean = pmean, mu_sd = psd)
plot_prior(priorEat0,design=design_Eat0,layout=c(3,5))

samplers <- make_emc(forstmann,design_Eat0,prior=priorEat0)
# save(samplers,file="Hierarchical/fits_LBA/Eat0.RData")

###  B + v emphasis effects ----

# NB: Following our DDM models, in the next two model we use a different 
#     contrast matrix for E when applied to v (E2mat) and B (Emat). 
E2mat <- cbind(neutral=c(0,1,0),speed=c(1,0,0)); E2mat

design_Eav <- design(data=forstmann,model=LBA,
  matchfun=function(d)d$S==d$lR,
  contrasts=list(v=list(E=E2mat,lM=ADmat),B=list(E=Emat),
                 sv=list(lM=ADmat)),
  formula=list(v~E/lM,B~E+lR,A~1,t0~1,sv~lM),
  constants=c(sv=log(1)))

pmean = c(
  v=3,'v_Eneutral'=0,'v_Espeed'=0,
  'v_Eaccuracy:lMd'=2,'v_Eneutral:lMd'=2,'v_Espeed:lMd'=2,
  B=log(2),'B_Ea-n'=0,'B_Ea-s'=0,'B_lRright'=0,
  A=log(1),t0=log(.2),
  sv_lMd=0)
psd = c(2,2,2,2,2,2,1,.25,.25,.25,0.5,.5,1)
priorEav <- prior(design_Eav,mu_mean = pmean, mu_sd = psd)

plot_prior(priorEav,design=design_Eav,layout=c(3,6))

samplers <- make_emc(forstmann,design_Eav,prior=priorEav)
# save(samplers,file="Hierarchical/fits_LBA/Eav.RData")

### B + t0 + v emphasis effects ----

# Emphasis effect on both v and t0
design_Eavt0 <- design(data=forstmann,model=LBA,
  matchfun=function(d)d$S==d$lR,
  contrasts=list(v=list(E=E2mat,lM=ADmat),B=list(E=Emat),
                 sv=list(lM=ADmat)),
  formula=list(v~E/lM,B~E+lR,A~1,t0~E,sv~lM),
  constants=c(sv=log(1)))

pmean = c(
  v=3,'v_Eneutral'=0,'v_Espeed'=0,'v_Eaccuracy:lMd'=2,'v_Eneutral:lMd'=2,'v_Espeed:lMd'=2,
  B=log(2),'B_Ea-n'=0,'B_Ea-s'=0,'B_lRright'=0,
  A=log(1),t0=log(.2),'t0_Eneutral'=0,'t0_Espeed'=0,sv_lMd=0)
psd = c(2,2,2,2,2,2,1,.25,.25,.25,0.5,.5,1,1,1)
priorEavt0 <- prior(design_Eavt0,mu_mean = pmean, mu_sd = psd)
plot_prior(priorEavt0,design=design_Eavt0,layout=c(3,6))

samplers <- make_emc(forstmann,design_Eavt0,prior=priorEavt0)
# save(samplers,file="Hierarchical/fits_LBA/Eavt0.RData")

#### Check fits ----

# As with the DDM we will leave correlation checks as an exercise.

# Old style, no bias, all good
print(load("Hierarchical/fits_LBA/B.RData"))
check(sB)

# Old style, additive bias, all good
print(load("Hierarchical/fits_LBA/BBa.RData"))
check(sBBa)

# Baseline (additive bias, sv match_, all good
print(load("Hierarchical/fits_LBA/BBasv.RData"))
check(sBBasv)

# Baseline + t0 varies with E, all good
print(load("Hierarchical/fits_LBA/Eat0.RData"))
check(sEat0)

# Baseline + v varies with E, all good
print(load("Hierarchical/fits_LBA/Eav.RData"))
check(sEav)

# Baseline, t0 and v vary with E. This is the only case with an Rhat > 1.1, but
# only minimally so.
print(load("Hierarchical/fits_LBA/Eavt0.RData"))
check(sEavt0)

#### Model comparison ----

# Both allowing response bias and particularly allowing sv to vary with lM 
# produce big jumps in fit. 

compare(list(B0=sB, B0a=sBBa, B=sBBasv,Bt0=sEat0, Bv=sEav, Bvt0=sEavt0),
             cores_per_prop = 3)
#          MD wMD    DIC wDIC   BPIC wBPIC EffectiveN  meanD  Dmean   minD
# B0   -14154   0 -15196    0 -15031     0        165 -15361 -15488 -15526
# B0a  -14331   0 -15528    0 -15308     0        220 -15748 -15900 -15968
# B    -15381   1 -16758    0 -16575     0        183 -16941 -17078 -17124
# Bt0  -15160   0 -16843    0 -16587     0        256 -17099 -17243 -17356
# Bv   -15297   0 -17169    1 -16877     1        292 -17461 -17684 -17753
# Bvt0 -14879   0 -17062    0 -16669     0        393 -17456 -17717 -17849

# MD selects the conventional model with these two features

# DIC and BPIC select the model that adds only a v effect, although the most
# complex model is a close second.

# Compared to the DDM, all LBA models are hugely better. Even the worst LBA is 
# better than the best DDM, even though it has only 9 parameters compared to 
# the DDM's 15.
compare(list(B=sBBasv,Eav=sEav,Evt0=sEavt0,DDMvt0=sDDMvt0),
             cores_per_prop = 3)
#            MD wMD    DIC wDIC   BPIC wBPIC EffectiveN  meanD  Dmean   minD
# B      -15386   1 -16758    0 -16575     0        183 -16941 -17078 -17124
# Eav    -15293   0 -17169    1 -16877     1        292 -17461 -17684 -17753
# Evt0   -14875   0 -17062    0 -16669     0        393 -17456 -17717 -17849
# DDMvt0 -11648   0 -14482    0 -14053     0        428 -14910 -15223 -15338

# Following we focus on the conventional LBA model picked by MD and the most
# complex model. In exercises at the end of this lesson we consider the LBA
# model selected by the ICs.

#### Goodness of fit of the best LBA models ----

# # First obtain posterior predictives 
# ppLBAic <- predict(sEavt0,n_cores=12)
# ppLBAmd <- predict(sBBasv,n_cores=12)
# save(ppLBAic,ppLBAmd,file="Hierarchical/fits_LBA/ppLBA.RData")

# We first compare the DDM and two LBA selected models.

# The most complex DDM displays very little misfit.
plot_fit(forstmann,ppLBAic,layout=c(2,3),factors=c("E","S"),lpos="right",xlim=c(.25,1.5))

# The conventional model shows much more misfit, particularly in accuracy
plot_fit(forstmann,ppLBAmd,layout=c(2,3),factors=c("E","S"),lpos="right",xlim=c(.25,1.5))

# The DDM misfit is much more in RT distribution.
plot_fit(forstmann,ppDDM,layout=c(2,3),factors=c("E","S"),lpos="right",xlim=c(.25,1.5))


# Clearly LBAic is descriptively much better, so we focus on it from here to see 
# how it fares in detail.

# Accuracy: good but missing a slightly different bias effect in speed not than
# the other conditions.
pc <- function(d) 100*mean(d$S==d$R)
tab <- plot_fit(forstmann,ppLBAic,layout=c(2,3),factors=c("E","S"),
                stat=pc,stat_name="Accuracy (%)",xlim=c(70,95))
round(tab,2)

# Note that we could also average over stimulus if our main interest was on the
# emphasis factor. The good fit here is becasue the bias effect averages out.
tab <- plot_fit(forstmann,ppLBAic,layout=c(1,3),factors=c("E"),
                stat=pc,stat_name="Accuracy (%)",xlim=c(70,95))
round(tab,2)


# Speed for correct responses

# Mean RT excellent
tab <- plot_fit(forstmann,ppLBAic,layout=c(2,3),factors=c("E","S"),
  stat=function(d){mean(d$rt[d$R==d$S])},stat_name="Mean Correct RT (s)",xlim=c(0.375,.6))
round(tab,2)

# Fast responses (10th percentile): also good
tab <- plot_fit(forstmann,ppLBAic,layout=c(2,3),factors=c("E","S"),xlim=c(0.275,.4),
  stat=function(d){quantile(d$rt[d$R==d$S],.1)},stat_name="10th Percentile Correct (s)")
round(tab,2)

# Slow responses (90th percentile): also good
tab <- plot_fit(forstmann,ppLBAic,layout=c(2,3),factors=c("E","S"),xlim=c(0.525,.975),
  stat=function(d){quantile(d$rt[d$R==d$S],.9)},stat_name="90th Correct Percentile (s)")
round(tab,2)

# RT variability (SD): hence variability pretty good
tab <- plot_fit(forstmann,ppLBAic,layout=c(2,3),factors=c("E","S"),
  stat=function(d){sd(d$rt[d$R==d$S])},stat_name="SD Correct (s)",xlim=c(0.1,.325))
round(tab,3)

# Errors speed: Even error speed well accommodated.
tab <- plot_fit(forstmann,ppLBAic,layout=c(2,3),factors=c("E","S"),xlim=c(0.375,.725),
  stat=function(d){mean(d$rt[d$R!=d$S])},stat_name="Mean Error RT (s)")
round(tab,2)


#### Parameter estimates for the most complex model ----

### Population means ----

# Priors mostly well dominated for population means
plot_pars(sEavt0,layout=c(3,5))

# Effects are easiest to interpret on the natural scales
posterior_summary(sEavt0,map=TRUE)

# v) Accuracy emphasis mainly affects FALSE (mismatch) rates, which decrease 
#    markedly, while TRUE (match) rates don't differ much. 

# B) There is an overall left response bias, and a marked shift up in thresholds
#    from speed to neutral, but perhaps a decrease from neutral to accuracy.

# t0) Although model selection favored an effect none is evident in the 
#     posterior estimates.

# sv) FALSE is much more variable than TRUE

# A) About half of B.


## Rate quality (match - mismatch) ----

# The overall picture is rate quality fairly similar for accuracy and neutral and
# much less for speed.

# Credibly Accuracy > Neutral, but only positive Bayes Factor evidence
credible(sEavt0,c("v_Eaccuracy:lMd","v_Eneutral:lMd"))
hypothesis(sEavt0,fun=\(x) diff(x[c("v_Eaccuracy:lMd","v_Eneutral:lMd")]))

# Credibly Neutral > Speed, with very strong evidence.
credible(sEavt0,c("v_Eneutral:lMd","v_Espeed:lMd"))
hypothesis(sEavt0,fun=\(x) diff(x[c("v_Eneutral:lMd","v_Espeed:lMd")]))

# NB: Because hypothesis uses sampling it can sometimes produce unstable results
#     that are particularly evident when the result is large, as in the last 
#     test. More stable results can be obtained by increasing the number of 
#     samples from the default 1e4, e.g.:
hypothesis(sEavt0,fun=\(x) diff(x[c("v_Eneutral:lMd","v_Espeed:lMd")]),N=1e5)

## Rate quantity ----

# Credibly accuracy > neutral, and very strong Bayes Factor
credible(sEavt0,"v_Eneutral")
hypothesis(sEavt0,"v_Eneutral")

# Accuracy > speed not credible and equivocal evidence of a difference
credible(sEavt0,"v_Espeed")
hypothesis(sEavt0,"v_Espeed")

# We could also test if the two quantity effects differ. Although speed is
# credibly greater than accuracy there is no Bayes Factor support.
credible(sEavt0,c("v_Espeed","v_Eneutral"))
hypothesis(sEavt0,fun=\(x) diff(x[c("v_Espeed","v_Eneutral")]))

## Response bias ----

# The left bias is just credible (the same interval is reported in the overall
# posterior_summary(sEavt0) table, but the Bayes Factor is equivocal.
credible(sEavt0,c("B_lRright"))
hypothesis(sEavt0,"B_lRright")

## Caution ----

# We will perform the caution tests on the natural scale as that is easier to
# interpret, although it requires fun's that are a little more complicated as 
# they must average over left and right.

# Caution is credibly much bigger in neutral than speed, with strong Bayes 
# factor support

fun=\(x) {
  diff(c(mean(x[c("B_Espeed_lRleft","B_Espeed_lRright")]),
         mean(x[c("B_Eneutral_lRleft","B_Eneutral_lRright")])))
}
credible(sEavt0,map=TRUE,x_fun=fun)
hypothesis(sEavt0,map=TRUE,fun=fun)

# Caution is not quite credibly greater in neutral than accuracy, and this
# difference is equivocal according to the Bayes Factor.

fun=\(x) {
  diff(c(mean(x[c("B_Eaccuracy_lRleft","B_Eaccuracy_lRright")]),
         mean(x[c("B_Eneutral_lRleft","B_Eneutral_lRright")])))
}

credible(sEavt0,map=TRUE,x_fun=fun)
hypothesis(sEavt0,map=TRUE,fun=fun)

## Non-decision time ----

# t0 has a small increase between neutral and the other two conditions. This
# is best appreciated on the natural scale so we make tests there

# The increase for neutral over speed is just not credible at 15ms, and Bayes 
# Factor positively supports no difference
credible(sEavt0,c("t0_Espeed","t0_Eneutral"),map=TRUE,digits=3)
hypothesis(sEavt0,fun=\(x)diff(x[c("t0_Espeed","t0_Eneutral")]),map=TRUE)

# The increase for neutral over accuracy just fails credibility, at 15ms, and 
# Bayes Factor supports no difference
credible(sEavt0,c("t0_Eaccuracy","t0_Eneutral"),map=TRUE,digits=3)
hypothesis(sEavt0,fun=\(x)diff(x[c("t0_Eaccuracy","t0_Eneutral")]),map=TRUE)

## Rate variability ----

# Reflecting the strong model selection results, sv is credibly much 
# greater for FALSE than TRUE. There is essentially no chance of
# quality under the posterior so the Bayes Factor is infinite.
credible(sEavt0,c("sv_lMFALSE","sv_lMTRUE"),map=TRUE,digits=3)
hypothesis(sEavt0,fun=\(x)diff(x[c("sv_lMFALSE","sv_lMTRUE")]),map=TRUE)

## Individual tests ----

# Note that we can also look at credibility on an individual subject , e.g., for 
# subject 2 the accuracy to neutral "quantity" (average over match and mismatch) 
# rate increase (as coded by the E2mat contrasts) is less than the accuracy
# to speed quantity change, but not for subject 1 
credible(sEavt0,c("v_Eneutral","v_Espeed"),selection="alpha",x_subject=1)
credible(sEavt0,c("v_Eneutral","v_Espeed"),selection="alpha",x_subject=2)

# Similarly can look at the credibility of between-subjects tests differences, 
# e.g., the quantity effect for accuracy for subject 2 than subject 1.
credible(x=sEavt0,x_name="v_Eaccuracy:lMd",x_subject=2,
         y=sEav,y_name="v_Eaccuracy:lMd",y_subject=1,selection="alpha")

# Can also make mapped comparisons, e.g., t0 is greater for subject 2 than 1
credible(x=sEavt0,selection="alpha",
       x_name="t0",x_subject=2,y=sEavt0,y_name="t0",y_subject=1)

# NB: hypothesis cannot be used with individuals in a hierarchical model

### Population variability and correlations ----

# NB: No mapped tests for correlations and variances!

# The level of individual differences is estimated by the variance parameters
# at the population level.

# Here priors are not so dominated for variance as the level of individual 
# differences is harder to estimate with only 19 participants
plot_pars(sEavt0,selection="sigma2", use_prior_lim = FALSE,layout=c(3,5))
posterior_summary(sEav, selection = "sigma2")

# Can look at the credibility of variance differences, e.g., are individual 
# differences in the accuracy-speed emphasis effect on rate quantity greater 
# than the accuracy-neutral?

# The difference for neutral is credibly greater than for speed with strong 
# Bayes factor support.
credible(sEavt0,c("v_Espeed","v_Eneutral"),selection="sigma2")
hypothesis(sEavt0,fun=\(x)diff(x[c("v_Espeed","v_Eneutral")]),selection="sigma2")
       
# Finally, we can also test correlations. For example quality is highly 
# positively correlated, with the correlation of accuracy and neutral .979,
# whereas the correlation of accuracy and speed is .893. Lets test the 
# difference.
credible(sEavt0,c("v_Eaccuracy:lMd.v_Eneutral:lMd","v_Eaccuracy:lMd.v_Espeed:lMd"),
         selection="correlation")
hypothesis(sEavt0,selection="correlation",
  fun=\(x)diff(x[c("v_Eaccuracy:lMd.v_Eneutral:lMd","v_Eaccuracy:lMd.v_Espeed:lMd")]))

# Although just credibly different there is not much support from the test. 


#### Using custom factors ----

# Above we say that the accuracy - neutral effect is fairly weak for
# thresholds in the winning model, so you might want to test a version of this 
# model in which the effect of these two levels are not differentiated. One way 
# to do this is to make a new factor (lets call it E2) that splits the data into 
# "speed" and "notspeed".
levels(factor(forstmann$E=="speed",labels=c("notspeed","speed")))

# However, you might have also noticed that there seems to be a
# stronger accuracy vs. neutral difference in rates so you might want to keep
# that difference (i.e., keep using the E factor for rates). You can't have
# both E and E2 in factors because that would be redundant, so instead you
# can keep E in factors and use the functions argument of design to
# make E2 available. 

# In general, functions specified a list of functions that make factors with
# arbitrary combinations of the cells of the full design specified by 
# crossing of the elements of factors (i.e., the elements of factors "span"
# every possible difference that is relevant to your full model). 

# In the resulting model we use the following contrast for the (2-level) ES
# factor, so the intercept is the threshold for speed and 'ns-s' is  
# notspeed-speed (which we would expect to be positive) 
ESmat <- matrix(c(1,0),ncol=1)
dimnames(ESmat)[[2]] <- 'ns-s'

# As a result we drop one parameter
design_Eavt02 <- design(data=forstmann,model=LBA,
  matchfun=function(d)d$S==d$lR,
  contrasts=list(v=list(E=E2mat,lM=ADmat),B=list(E2=ESmat),
                 sv=list(lM=ADmat)),
  functions=list(E2=function(d){factor(d$E=="speed",labels=c("notspeed","speed"))}),
  formula=list(v~E/lM,B~E2+lR,A~1,t0~E,sv~lM),
  constants=c(sv=log(1)))

pmean = c(
  v=3,'v_Eneutral'=0,'v_Espeed'=0,'v_Eaccuracy:lMd'=2,'v_Eneutral:lMd'=2,'v_Espeed:lMd'=2,
  B=log(2),'B_Ens-s'=0,'B_lRright'=0,
  A=log(1),t0=log(.2),'t0_Eneutral'=0,'t0_Espeed'=0,sv_lMd=0)
psd = c(2,2,2,2,2,2,1,.25,.25,0.5,.5,1,1,1)
priorEavt02 <- prior(design_Eavt02,mu_mean = pmean, mu_sd = psd)
plot_prior(priorEavt02,design_Eavt02,layout=c(2,8))

samplers <- make_emc(forstmann,design_Eavt02,prior=priorEavt02)
# save(samplers,file="Hierarchical/fits_LBA/Eavt02.RData")

print(load("Hierarchical/fits_LBA/Eavt02.RData"))
# Good convergence
check(sEavt02)

# Simpler model wins 
compare(list(E=sEavt0,E2=sEavt02),cores_for_props = 3)
#        MD wMD    DIC wDIC   BPIC wBPIC EffectiveN  meanD  Dmean   minD
# E  -14882   0 -17062    0 -16669     0        393 -17456 -17717 -17849
# E2 -15098   1 -17273    1 -16953     1        320 -17593 -17782 -17913


# In this model t0 effects appear weak
posterior_summary(sEavt02)

# t0 effects very similar in neutral and speed
credible(sEavt02,c("t0_Eneutral","t0_Eaccuracy"))

# strong evidence against each effect
1/hypothesis(sEavt02,"t0_Eneutral")
1/hypothesis(sEavt02,"t0_Eaccuracy")

# even stronger evidence for no difference between effects
1/hypothesis(sEavt02,fun=\(x) diff(x[c("t0_Eneutral","t0_Eaccuracy")]))

#### Exercises ----

# 1. The t0 effects in the last LBA models are small, so fit a model in which
#    they are dropped (i.e., t0 ~ 1). Note that this model is the same as the
#    Eav model, but with the custom factor so thresholds only differ between
#    speed and non-speed, so lets call it Eav2
# a) Compare this model to the Eavt0, Eavt02 and Eav models: which win 
#    according to each model selection measure.
# b) How does it fit the data (use the same plots as for the most complicated
#    model).
#
# 2. Earlier we noted that there may be a difference in response bias between
#    accuracy and non-accuracy conditions. Allow this difference on top of the
#    Eav2 and Eavt02 models (call them Eav3 and Eavt03).
# a) Compare these models to the Eavt0, Eavt02, Eav and Eav2 models: which win 
#    according to each model selection measure.
# b) Do these models address any misfit seen with the earlier models?
#
# 3. The Eav3 model assumes three things explain differences among conditions.
#    For each use posterior testing to characterize the causes of the differences
#    in terms of model parameters on the natural scale. Report results in terms 
#    of median differences with 95% CIs (in []) and Bayes Factor (BF) tests 
#    using conventional positive (> 3, < 10) and strong (>10) classifications.
#    a) Neutral vs. accuracy is purely due to accumulation rate (quantity and 
#       quality).
#    b) Left vs. right is purely due to threshold bias.
#    c) Speed differs from the other two emphasis conditions due to both 
#       rates and thresholds. 

