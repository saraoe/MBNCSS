rm(list = ls())
library(EMC2)
load("Joint/data_conflict4.RData")

#
# In this joint modelling example we'll look at four decision-making tasks
# That all employ a cognitive conflict manipulation using congruent and
# incongruent stimuli. The data is stored as a list of data sets
# An attention network task
plot_defective_density(data_list[[1]], factors = c("stim", "conflict_type"))
# A cued forgetting task
plot_defective_density(data_list[[2]], factors = c("stim", "conflict_type"))
# A location matching task
plot_defective_density(data_list[[3]], factors = c("stim", "conflict_type"))
# A simon task
plot_defective_density(data_list[[4]], factors = c("stim", "conflict_type"))

# Let's set up our design for the DDM we'll be using
ADmat <- matrix(c(-1/2,1/2),ncol=1,dimnames=list(NULL,"_d"))
Smat <- matrix(c(-1,1),ncol=1,dimnames=list(NULL,"_diff"))

# Since the design we can use across tasks is the same we can just create
# one design and use the same design across all tasks.
design_conflict <- design(data = data_list[[1]], contrasts = list(v = list(stim= Smat, conflict_type = ADmat)),
                          formula =list(v~stim*conflict_type,a~1, t0~1, s ~ 1, Z ~ 1),
                          constants=c(s=log(1), v_conflict_type_d = 0),
                          model = DDM)

# And then to make a joint design, simply replicate it across tasks
designs <- rep(list(design_conflict), 4)

# If the designs were different to each other we could have set it up as such:
# designs <- list(design1, design2, design3, design4)

# Let's set up a prior, here we use the same mean prior across tasks.
# We'll first fit every subject individually, to illustrate the effect of attenuation
prior_single <- prior(design = designs, type = 'single',
                      pmean = rep(c(1, 1, .5, log(1.5), log(.2), qnorm(.5)), 4))

# joint_single <- make_emc(data_list, designs, prior_list = prior_single, type = "single")
# joint_single <- fit(joint_single, cores_per_chain = 8, fileName = "joint_single.RData")
# save(joint_single, file = "joint_single.RData")
load("Joint/joint_single.RData")
# In a typical application, the parameter estimates of each person would be
# correlated afterwards. This is how we would do that
alpha_single <- posterior_summary(joint_single, selection = "alpha", prob = .5)
# Bind across parameters
alpha_single2 <- do.call(cbind, alpha_single)
# Take the correlation
cor_alpha <- cor(alpha_single2)
# Let's make some nice names we can keep on using throughout the code
nice_names <- rownames(cor_alpha)
nice_names[grepl("diff", nice_names) & !grepl("conflict", nice_names)] <- paste0(1:4, "|dS")
nice_names[grepl("conflict", nice_names)] <- paste0(1:4, "|dC")
rownames(cor_alpha) <- colnames(cor_alpha) <- nice_names

# We'll calculate p-values for these single subject correlations
# So that using the corrplot package we can only plot significant correlations
pval <- psych::corr.test(alpha_single2, adjust="none")$p
rownames(pval) <- colnames(pval) <- nice_names
par(mfrow = c(1,1))
corrplot::corrplot(cor_alpha, tl.col = "black", p.mat = pval, insig = "blank")

# Now we're going to be plotting these correlations the proper way
# using a joint hierarchical model.
# Let's first set up a prior for the hierarchical model
prior_standard <- prior(design = designs, type = 'standard',
                     mu_mean = rep(c(0, 2, .5, log(1.5), log(.2), qnorm(.5)), 4))

# This will plot a whole bunch of correlations, they're priors are all uniform
plot_prior(prior_standard, design = designs,
           selection = "correlation", N = 1e4)

# Construct and run our joint model
# joint_standard <- make_emc(data_list, designs,
#                            prior_list = prior_standard)
#
# joint_standard <- fit(joint_standard, cores_per_chain = 10, fileName = "joint_standard.RData")
# save(joint_standard, file = "joint_standard.RData")
load("Joint/joint_standard.RData")
# We can use the plot_relations function to replicate our previous plot
# But now for the joint model, determining credible correlations using the
# 95% posterior density.
plot_relations(joint_standard, only_cred = T, plot_cred = F, plot_means = F,
               nice_names = nice_names)

# First let's take a look at the attenuation factor, by looking at the variances.
# For this we'll use the recovery function in a slightly different way
recovery(joint_standard, diag(var(alpha_single2)), selection = "sigma2", ci_plot_args = list(col = "brown"),
         xlab = "Non-hierarchical", ylab = "Hierarchical", xlim = c(0, 0.35), ylim = c(0, .35))

# As we can see, almost all variances are larger for the non-hierarchical model
# indicating attenuation and shrinkage

# Despite that, we have high overlap between the joint model and the non-hierarchical model
# Sometimes correlations are even larger in the non-hierarchical model, this is probably
# because of the sloppy correlations getting in for the non-hierarchical model (i.e. z and v)
# Remember Michelle's talk?
recovery(joint_standard, cor_alpha, selection = "correlation", ci_plot_args = list(col = "brown"))

# See the strong correlations between v and Z in the single subject model,
# They might also affect the correlations between v and Z in the between subject correlations
pairs_posterior(joint_single, use_par = c("1|v", "1|Z"))

# Let's use our hierarchical model to perform some inference
# First let's plot the fit
# pps <- predict(joint_standard, n_cores = 25)
# save(pps, file = "pp_standard.RData")
load("Joint/pp_standard.RData")

plot_fit(data_list[[1]], pps[[1]], factors = c("stim", "conflict_type"), layout = c(2,2))
plot_fit(data_list[[2]], pps[[2]], factors = c("stim", "conflict_type"), layout = c(2,2))
plot_fit(data_list[[3]], pps[[3]], factors = c("stim", "conflict_type"), layout = c(2,2))
plot_fit(data_list[[4]], pps[[4]], factors = c("stim", "conflict_type"), layout = c(2,2))

## Pretty poor fit, but the aim here is not model comparison of the subject-level model

# However, with our full-hierarchical model, we've estimated [(4*6)*(4*6) - 4*6]/2 = 276 (!!!) correlations
# 4 being the number of tasks, 6 being the number of parameters.
# That's more correlations than we have participants, one could argue we're overfitting the structural
# relationships between participants.

# One way to massively reduce that dimensionality is to use a blocked covariance matrix
# To do so in EMC, we specify which parameters should be blocked together.
# Here we've chosen to block the parameters per `type`
par_groups <- rep(1:6, 4)
names(par_groups) <- nice_names
par_groups
prior_blocked <- prior(design = designs, type = 'blocked',
                       mu_mean = rep(c(1, 1, .5, log(1.5), log(.2), qnorm(.5)), 4),
                       par_groups = par_groups)
# As you can see only the correlations between types is estimated now
plot_prior(prior_blocked, design = designs,
           selection = "correlation", par_groups = par_groups, N = 1e4)

# Now let's make and estimated our blocked model again:
# joint_blocked <- make_emc(data_list, designs,
#                           prior_list = prior_blocked, type = "blocked",
#                           par_groups = par_groups)
#
# joint_blocked <- fit(joint_blocked, cores_per_chain = 10,fileName = "joint_blocked.RData")
# save(joint_blocked, file = "joint_blocked.RData")
load("Joint/joint_blocked.RData")

# First compare the correlations between the two types of models
par(mfrow = c(1,2))
plot_relations(joint_standard, only_cred = T, plot_cred = F,
               nice_names = nice_names)
plot_relations(joint_blocked, only_cred = T, plot_cred = F,
               nice_names = nice_names)
# Again no correlations between the dC parameters

# A more sophisticated way to reduce the dimensionality is using factor analysis

# Hierarchical factor analysis can also be done in EMC, by specifying a constraint matrix
# We use factor analysis to our advantage to make a more constraint model
# to test certain hypotheses.
p_vector <- sampled_p_vector(designs)
p_names <- names(p_vector)
Lambda_mat <- matrix(0, length(p_vector), 4)
Lambda_mat[grepl("diff", p_names) &! grepl("conflict", p_names),1] <- Inf
Lambda_mat[grepl("conflict", p_names),2] <- Inf
Lambda_mat[grepl("a", p_names),3] <- Inf
Lambda_mat[grepl("t0", p_names),4] <- Inf
#
Lambda_mat[2,1] <- Lambda_mat[3,2] <- Lambda_mat[4,3] <- Lambda_mat[5,4] <- 1
colnames(Lambda_mat) <- c("sC", "sD", "a", "t0")
rownames(Lambda_mat) <- p_names

# The Infs here represent freely estimated parameters
# The other items are constrained at that value
Lambda_mat

prior_factor <- prior(design = designs, type = 'factor',
                      mu_mean = rep(c(0, 2, .5, log(1.5), log(.2), qnorm(.5)), 4),
                      Lambda_mat = Lambda_mat)

plot_prior(prior_factor, design = designs,
           selection = "correlation", Lambda_mat = Lambda_mat)

plot_prior(prior_factor, design = designs,
           selection = "loadings", Lambda_mat = Lambda_mat)

plot_prior(prior_factor, design = designs,
           selection = "residuals", Lambda_mat = Lambda_mat)


# joint_factor <- make_emc(data_list, designs,
#                          prior_list = prior_factor, type = "factor",
#                          Lambda_mat = Lambda_mat)
#
# joint_factor <- fit(joint_factor, cores_per_chain = 10, fileName = "joint_factor_4F.RData")
# save(joint_factor, file = "joint_factor_4F.RData")
load("Joint/joint_factor_4F.RData")
# Let's compare the group-level mean estimated with the joint factor model to the group-level mean
# of the blocked model
# In green the factor model, in black the blocked model. They're in high agreement
plot_pars(joint_blocked, true_pars = joint_factor, selection = "mu",
          use_prior_lim = F)
# Unfortunately we cannot compare correlations using the plot_pars function for
# models with different correlations estimated, but that will be implemented soon


# We can use FA to test support for the existence of a latent construct using
# model comparison. To do so we set up a factor model that doesn't allow for a dC factor
#
p_vector <- sampled_p_vector(designs)
p_names <- names(p_vector)
Lambda_mat <- matrix(0, length(p_vector), 3)
Lambda_mat[grepl("diff", p_names) &! grepl("conflict", p_names),1] <- Inf
Lambda_mat[grepl("a", p_names),2] <- Inf
Lambda_mat[grepl("t0", p_names),3] <- Inf

Lambda_mat[2,1] <- Lambda_mat[4,2] <- Lambda_mat[5,3] <- 1
colnames(Lambda_mat) <- c("sC", "a", "t0")
rownames(Lambda_mat) <- p_names

prior_factor <- prior(design = designs, type = 'factor',
                      mu_mean = rep(c(0, 2, .5, log(1.5), log(.2), qnorm(.5)), 4),
                      Lambda_mat = Lambda_mat)

plot_prior(prior_factor, design = designs,
           selection = "correlation", Lambda_mat = Lambda_mat)

# joint_factor_no_sD <- make_emc(data_list, designs,
#                          prior_list = prior_factor, type = "factor",
#                          Lambda_mat = Lambda_mat)
#
# joint_factor_no_sD <- fit(joint_factor_no_sD, cores_per_chain = 10, fileName = "joint_factor_3F_no_sD.RData")
# save(joint_factor_no_sD, file = "joint_factor_3F_no_sD.RData")
load("Joint/joint_factor_3F_no_sD.RData")
# Now we can test if there's support for this additional factor through model comparison.
# Note that we can also run this with BayesFactor, but we need a lot more MCMC iterations to estimate BayesFactors
# for such massive models. Note that we could also compare the full model to the FA models,
# but for that we would need BayesFactors to more accurately incorporate the massive prior space of the full model
compare(list(F4 = joint_factor, F3 = joint_factor_no_sD), BayesFactor = F)
# So apparently there's no latent conflict processing factor, since the
# model that did not include a latent factor for conflict processing across tasks
# is not supported by the data.


# Now that we know there's no support for an sD (congruency drift rate) factor
# we can also see if there's support for a sC (general drift rate) factor
p_vector <- sampled_p_vector(designs)
p_names <- names(p_vector)
Lambda_mat <- matrix(0, length(p_vector), 2)
Lambda_mat[grepl("a", p_names),1] <- Inf
Lambda_mat[grepl("t0", p_names),2] <- Inf

Lambda_mat[4,1] <- Lambda_mat[5,2] <- 1
colnames(Lambda_mat) <- c("a", "t0")
rownames(Lambda_mat) <- p_names

#
prior_factor <- prior(design = designs, type = 'factor',
                            mu_mean = rep(c(0, 2, .5, log(1.5), log(.2), qnorm(.5)), 4),
                            Lambda_mat = Lambda_mat)

plot_prior(prior_factor, design = designs,
           selection = "correlation", Lambda_mat = Lambda_mat)

# joint_factor_no_sC <- make_emc(data_list, designs,
#                                prior_list = prior_factor, type = "factor",
#                                Lambda_mat = Lambda_mat)
#
# joint_factor_no_sC <- fit(joint_factor_no_sC, cores_per_chain = 10, fileName = "joint_factor_2F_no_sC.RData")
# save(joint_factor_no_sC, file = "joint_factor_2F_no_sC.RData")
load("Joint/joint_factor_2F_no_sC.RData")
compare(list(F4 = joint_factor, F3 = joint_factor_no_sD,
             F2 = joint_factor_no_sC), BayesFactor = F)

# And as we can see there is support for a general sD factor
# Now we can do some inference on our winning model
loadings <- EMC2:::standardize_loadings(joint_factor_no_sD)
plot_relations(loadings = loadings, only_cred = T, plot_cred = F,
               nice_names = nice_names)

plot_pars(joint_factor_no_sD, selection = "loadings", use_prior_lim = F)
plot_pars(joint_factor_no_sD, selection = "loadings", use_prior_lim = F, flatten = TRUE)
hypothesis(joint_factor_no_sD, parameter = "3|a.a", selection = "loadings")

