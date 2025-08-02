rm(list = ls())
library(EMC2)
dat <- forstmann

# This script is heavily based on the LBA design of the EMC2 paper. 
# I've kept relevant descriptions of that tutorial in here, but will
# not go into detail on that. 


### Design Specification ---

# We will now try to incorporate our previous findings in an LBA model.
# While we would normally do full model development for the DDM and LBA
# separately, this tutorial will focus on only part of that process for 
# the LBA for illustrative purposes.
#
# As the LBA is a race model it has one accumulator representing each
# possible response, in this case separate accumulators for the left and
# right response. EMC2 uses the levels of the "R" factor found in the data
# to construct a factor which denotes the accumulators "lR" (the latent
# response) with level names taken from Rlevels.
#
# The lR factor allows for response bias, analogous to Z in the DDM. In
# race models, response bias can be modeled by allowing different
# thresholds (B) for different responses. Here we will use B~lR to allow
# for different thresholds for the accumulator corresponding to left and
# right stimuli (e.g., a bias to respond left occurs if the left threshold
# is lower than the right threshold).
#
# For race models, you often need to provide a "matchfun" function that
# marks which accumulator matches the correct response. Typically, this
# function takes the latent response (lR) and returns TRUE when it matches
# the stimulus (S):

matchfun=function(data)data$S==data$lR

# matchfun is used to automatically to create a latent match (lM) factor
# with levels "FALSE" (i.e., the stimulus does not match the accumulator)
# and "TRUE" (i.e., the stimulus does match the accumulator). This is
# added internally and can also be used in model formula, typically for
# parameters related to the rate of accumulation.
#
# When working with lM, it's often useful to design an "average and
# difference" contrast matrix. For binary responses, this has a simple
# form:

ADmat <- cbind(d = c(-1/2,1/2))
ADmat

# When applied to drift rates (v~lM), the contrast parameter v_lMd
# represents the difference between correct and incorrect response drift
# rates. This "rate quality" parameter is analogous to the DDM drift rate.
# When accuracy exceeds chance, rate quality is positive.
#
# The ADmat intercept (v) represents the average accumulation rate across
# accumulators, measuring both stimulus magnitude and decision urgency
# effects. Unlike the DDM, which represents only evidence differences
# between possible responses, this "rate quantity" parameter captures the
# overall tendency to respond. The DDM does not directly represent this
# tendency.
#
# For simplicity, we omit rate bias, although it could be included by
# allowing v~S*lM.
#
# We also allow the drift rate standard deviation to vary as a function of
# match (sv~lM), as has become standard in more recent applications of
# the LBA.
#
# Just as with the DDM, we fix one parameter to make the model
# identifiable; in this case, as is common, we chose to fix drift rate
# standard deviation intercept (sv).
#
# Given the findings of the hierarchical DDM that there's no evidence for
# a threshold or drift rate difference between neutral and accuracy
# emphasis, we construct our model to estimate one shared E parameter for
# the neutral and accuracy condition. To achieve this, we use the
# functions argument.

E2 <- function(data) factor(data$E!="speed",labels=c("speed","nonspeed"))
design_LBABvE <- design(data = dat,model=LBA,matchfun=matchfun,
                        formula=list(v~lM*E2,sv~lM,B~E2+lR,A~1,t0~1),
                        contrasts=list(lM = ADmat),constants=c(sv=log(1)),
                        functions=list(E2=E2))

# Here the E2 function creates a new factor in the data, with levels speed
# and nonspeed.
#
# Whereas the drift rate intercept in the DDM represents response bias,
# which is unlikely to be affected by emphasis instructions, in the LBA
# the intercept codes for urgency effects, which are more likely to be
# affected by emphasis instructions. We therefore allow the intercept to
# also vary by E, yielding a V_E2nonspeed parameter and a v_LMd:E2nonspeed
# parameter. The latter parameter represents the difference in
# accumulation rate between correct and incorrect responses under
# neutral/accuracy relative to the speed emphasis condition.
#
# Let's examine how our parameters map to the design:

mapped_pars(design_LBABvE)

### Prior Specification ----

# We now create a prior. Note that all parameters except v are estimated
# on a log scale (see ?LBA). For the LBA, t0 is typically smaller than for
# the DDM, so we set it to 0.2. We further set increasing thresholds and
# urgency for higher response caution conditions. As we have no response
# bias expectation, B_lRright is set to zero.

mu_mean <- c(v=1, v_E2nonspeed = -.2, v_lMd=1, "v_lMd:E2nonspeed"=.2,
             sv_lMd=log(1),B=log(1), B_E2nonspeed = log(1.5), B_lRright=0,
             A=log(0.25), t0=log(.2))
mu_sd <- c(v=1, v_lMd=0.5, "v_lMd:E2nonspeed"=0.5,
           sv_lMd=.5,B=0.3, B_E2nonspeed=0.3,B_lRright=0.3, A=0.4, t0=.5)

prior_LBABvE <- prior(design_LBABvE, type = 'standard',mu_mean=mu_mean,
                      mu_sd=mu_sd)
plot(prior_LBABvE, layout = c(2,3))

plot_design(prior_LBABvE, factors = list(v = c("E2", "lM"), B = "E2"),
            plot_factor = "E2")

### Model Estimation ----

# Now we can fit the model:

LBABvE <- make_emc(dat,design_LBABvE, prior=prior_LBABvE)
LBABvE <- fit(LBABvE, cores_per_chain = 4, fileName="tmp.RData")
save(LBABvE,file="samples/LBABvE.RData")

load("samples/LBABvE.RData")

# The model shows good convergence and generally good efficiency:

check(LBABvE)

### Model Inference ----

## Model Fit ----

# We generate posterior predictions for the LBA model to see how well the
# model fits the data:

pp_LBABvE <- predict(LBABvE, n_cores = 12)
acc_fun <- function(data) data$S == data$R
plot_cdf(LBABvE, pp_LBABvE, factors = "E2", functions = c(correct = acc_fun),
         layout= c(1,2), defective_factor = "correct",
         legendpos = c("right","topleft"))

# As with single-subject, we can also compare arbitrary descriptives 
# on the real data to the predictives of both models. 
# For example differences in response time and
# error rates between the emphasis conditions:

drt <- function(data){
  all <- tapply(data$rt,data$E,mean)
  out <- c(all['neutral'] - all['speed'],all['accuracy'] - all['speed'])
  names(out) <- c("NTR-SPD", "ACC-SPD")
  return(out)
}

derr <- function(data){
  data$correct <- data$S == data$R
  all <- tapply(data$correct,data$E,mean)*100
  out <- c(all['neutral'] - all['speed'],all['accuracy'] - all['speed'])
  names(out) <- c("NTR-SPD", "ACC-SPD")
  return(out)
}

par(mfrow = c(1,2))
plot_stat(LBABvE, list(LBA = pp_LBABvE), stat_fun = drt, 
  xlab = "RT (s) difference", layout = NULL,legendpos = c("topleft","topright"))
plot_stat(LBABvE, list(LBA = pp_LBABvE), stat_fun = derr, 
  xlab = "Accuracy (%) difference", layout = NULL,legendpos = c("topleft","topright"))

# So accuracy is quite similar between Accuracy and Neutral, 
# but response times are a little different, with people being slower in 
# the accuracy condition. This cannot be captured by the model, since it 
# forces these to be the same for these effects across all parameters. 

# One idea would be to allow urgency to vary for all levels of E. 
# But keep the other parameters the same (with only speed vs non-speed)
# Unfortunately this is non-trivial to achieve but we'll showcase step-by-step
# As a reminder this is our old design:
summary(design_LBABvE)

# Now we use a custom contrast matrix
E_incr <- contr.increasing(3)
colnames(E_incr) <- c("nonspeed", "acc")

# To make our design matrix for E:
# first illustrating with only E
design_LBABvE2 <- design(data = dat,model=LBA,matchfun=matchfun,
                        formula=list(v~E,sv~lM,B~E2+lR,A~1,t0~1),
                        contrasts=list(lM = ADmat, E = E_incr),
                        constants=c(sv=log(1)),
                        functions=list(E2=E2))

# Now with lM, we set v_Eacc:lMd to be a constant
design_LBABvE2 <- design(data = dat,model=LBA,matchfun=matchfun,
                        formula=list(v~E*lM,sv~lM,B~E2+lR,A~1,t0~1),
                        contrasts=list(lM = ADmat, E = E_incr),
                        constants=c(sv=log(1), `v_Eacc:lMd` = 0),
                        functions=list(E2=E2))

# This design satisfies our aim of having urgency vary by all levels of E
# whereas the drift rate difference only varies between speed non-speed
prior_LBABvE2 <- prior(design_LBABvE2, update = prior_LBABvE)
# # AH you had the wrong design
# prior_LBABvE2 <- prior(design_LBABvE, update = prior_LBABvE)

LBABvE2 <- make_emc(data = dat, design = design_LBABvE2, prior_list = prior_LBABvE2)

LBABvE2 <- fit(LBABvE2, cores_per_chain = 4)
save(LBABvE2, file = "samples/LBABvE2.RData")
load("samples/LBABvE2.RData")

# We'll first use the hypothesis function to test whether the new parameter 
# has support on the group-level
hypothesis(LBABvE2, parameter = "v_Eacc")
# The parameter has strong group-level support. 

# Just because there's a group-level difference doesn't necessarily mean that 
# there's also individual-level support and vice-versa. We can have only 
# fixed effects support, or only random effects support or both.
# 
# The hypothesis functions tests the group-level mean, whereas the compare 
# function tests whether there's a group-level mean AND individual
# differences support. 

# To illustrate, let's visualize the individual differences:
plot_pars(LBABvE2, all_subjects = TRUE)

# For the v_Eacc parameter we see large degree of individual differences, 
# but almost everyone has a negative effect, meaning lower urgency in the 
# accuracy condition compared to the neutral condition, similarly urgency
# is lower in the neutral condition compared to the speed condition. 

# Let's use bridge sampling to investigate the model support for these individual 
# differences
compare(list(old = LBABvE, new = LBABvE2), cores_per_prop = 3)

# The new version is overwhelmingly supported!

# We should also see better model fit in the posterior-predictives
# Create some new posterior predictives
pp_LBABvE2 <- predict(LBABvE2, n_cores = 12)

par(mfrow = c(1,2))
plot_stat(LBABvE, list(old = pp_LBABvE, new = pp_LBABvE2), stat_fun = drt, 
  xlab = "RT (s) difference", layout = NULL,legendpos = c("topleft","topright"))
plot_stat(LBABvE, list(old = pp_LBABvE, new = pp_LBABvE2), stat_fun = derr, 
  xlab = "Accuracy (%) difference", layout = NULL,legendpos = c("topleft","topright"))

# As expected the new model better accounts for these discrepancies


# Another interesting one is the response bias:
plot_pars(LBABvE2, all_subjects = TRUE, use_par = "B_lRright")

# On average response bias is only slight positive:
hypothesis(LBABvE2, parameter = "B_lRright")
# And there's inconclusive evidence whether it should be included in the model
# But potentially because of the large degree of individual differences
# A model without it would be weaker overall. This you can test in the exercises. 


# Hierarchical Shrinkage ----------------------------------------------------
# Now with our winning model we'll showcase what the effect was of fitting it 
# hierarchically, on average the degree of individual differences should be lower
# because the parameter estimates are shrunk towards the group-level mean. 
# To make this comparison we'll also fit this model non-hierarchically.
# To that end we first set up a new prior, the ideal model would have a super
# vague prior, but that wouldn't converge.
prior_LBABvE2_single <- prior(LBABvE2, update = prior_LBABvE2, type = "single",
                              psd = 2)

# The reason for this vague prior, is that a direct prior on the single-subject
# parameters is super influential compared to an indirect prior on the group-level
# The most appropriate comparison uses a super vague prior on the single subject, 
# and then  an informed prior on the two-step group-level model you would run on the 
# posterior medians.

summary(prior_LBABvE2_single)

LBABvE2_single <- make_emc(dat, design_LBABvE2, prior_list = prior_LBABvE2_single)
LBABvE2_single <- fit(LBABvE2_single, cores_per_chain = 4)

# Checking we see the first 500 iterations are still settling in, so we
# remove them with the subset command.
check(LBABvE2_single)
LBABvE2_single <- subset(LBABvE2_single, filter = 500)
save(LBABvE2_single, file = "samples/LBABvE2_single.RData")
load("samples/LBABvE2_single.RData")

# Let's first inspect the visual fit of an individual subject
pp_LBABvE2_single <- predict(LBABvE2_single, n_cores = 12)

# To find our most extreme subject, let's look at the rts a plot of mean RT
# vs. accuracy
plot(aggregate(rt ~ subjects, dat, mean)[,2],
     aggregate(dat$S==dat$R ~ subjects, dat, mean)[,2],pch=letters)
# subjects 11 (= k)/(15=o) have very fast/slow rts but high/low accuracy

# AH you choose subject 11 on speed alone, but subject 13 is even faster

# For subject 11 we see shrinkage towards lower accuracy/faster for 
# hierarchical (blue below green), especially in the speed condition
plot_cdf(dat, post_predict = list(single = pp_LBABvE2_single, hier = pp_LBABvE2), 
         subject = 11, factors = "E", functions = list(correct = acc_fun),
         defective_factor = "correct", layout = c(1,1))

# Whereas for 15 we see shrinkage towards faster RTs (green below blue)
plot_cdf(dat, post_predict = list(single = pp_LBABvE2_single, hier = pp_LBABvE2), 
         subject = 15, factors = "E", functions = list(correct = acc_fun),
         defective_factor = "correct", layout = c(1,1))

# We can see the hierarchical model fits the individual data slightly worse
# since it's shrunk the individual estimates a bit towards the group-level mean

# We can also illustrate this effect of shrinkage, by looking at the variance
# of the estimates graphically. 

# To that end we'll use the recovery function (more on this later)
# plotting the hierarchical estimates on the y axis and the single estimates on
# the x-axis. We see that a few outlying and uncertian single estimates mean 
# that the y-axis (hierarhcial) range is shrunk relative to the x-axis (single)
# range.

# AH great plot, I added a little more explanation. 
recovery(LBABvE2, true_pars = LBABvE2_single, selection = "alpha", xlab = "Single",
         ylab = "Hierarchical")

# We can also illustrate this effect of shrinkage, by looking at the variance
# of the estimates numerically. 

# Let's first group the single-subject estimates:
alpha <- do.call(cbind, credint(LBABvE2_single, probs = .5))
# And then calculate the variances
round(apply(alpha, 2, sd),3)

# The group-level variance estimates are much smaller
credint(LBABvE2, selection = "sigma2", probs = .5)

# Let's see which model is preferred. The DIC and BPIC are only loosely
# informed by the group-level, the Bayes Factor (MD) is the gold
# standard for comparing group-level models
compare(list(single = LBABvE2_single, hier = LBABvE2), cores_for_props = 4)
# Unsurprisingly, the group-level model is preferred by the MD and the single 
# model by the DIC/BPIC
# AH This is one place I get a difference with you, DIC likes hier more in my
#    samples (by 21) vs. 72 the other way in your samples (my single EffN is 
#    positive whereas yours is negative). Maybe phrase as their being little 
#    difference? 
#    In both sets of samples absolute fit (minD and Dmean) are BETTER for hier
#    which you may get questions about. I guess worse fit for outliers is 
#    compensated for by better fit for the others??

# Here the MD (Marginal Deviance; used for Bayes Factors), DIC (Deviance
# Information Criterion) and BPIC (Bayesian Predictive Information
# Criterion) are measures of model fit. The lower the better. The weights
# (wX) are the probabilities of the models, the higher the better.

# AH Above you say more on recovery later, but that seems to be missing.

