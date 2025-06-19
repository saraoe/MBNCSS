#### Exercises ----

# 1) No t0 effect, and E2 for B

design_Eav2 <- design(data=forstmann,model=LBA,
  matchfun=function(d)d$S==d$lR,
  contrasts=list(v=list(E=E2mat,lM=ADmat),B=list(E2=ESmat),sv=list(lM=ADmat)),
  functions=list(E2=function(d){factor(d$E=="speed",labels=c("notspeed","speed"))}),
  formula=list(v~E/lM,B~E2+lR,A~1,t0~1,sv~lM),
  constants=c(sv=log(1)))

pmean = c(
  v=3,'v_Eneutral'=0,'v_Espeed'=0,'v_Eaccuracy:lMd'=2,'v_Eneutral:lMd'=2,'v_Espeed:lMd'=2,
  B=log(2),'B_Ens-s'=0,'B_lRright'=0,
  A=log(1),t0=log(.2),sv_lMd=0)
psd = c(2,2,2,2,2,2,1,.25,.25,0.5,.5,1)
priorEav2 <- prior(design_Eav2,mu_mean = pmean, mu_sd = psd)
plot_prior(priorEav2,design_Eav2,layout=c(2,8))

samplers <- make_emc(forstmann,design_Eav2,prior=priorEav2)
# save(samplers,file="Hierarchical/fits_LBA/Eav2.RData")

print(load("Hierarchical/fits_LBA/Eav2.RData"))

# a) Eav2 wins by MD, and Eavt02 by DIC and BPIC
compare(list(Eav=sEav,Eav2=sEav2,Eavt0=sEavt0,Eavt02=sEavt02),cores_per_prop = 3)
#            MD wMD    DIC wDIC   BPIC wBPIC EffectiveN  meanD  Dmean   minD
# Eav    -15299   0 -17169    0 -16877     0        292 -17461 -17684 -17753
# Eav2   -15614   1 -17196    0 -16932     0        265 -17461 -17650 -17725
# Eavt0  -14881   0 -17062    0 -16669     0        393 -17456 -17717 -17849
# Eavt02 -15100   0 -17273    1 -16953     1        320 -17593 -17782 -17913

# b) Looking at misses of 95% CI (see bottom for code to run fit plots) 
# Accuracy: under-predicts for speed right and over for speed left
# Correct RT: 10th and 90th Percentile: slightly over predict for accuracy left
# Error mean RT: slightly under-predict accuracy left



# 2) Response bias 
design_Eavt03 <- design(data=forstmann,model=LBA,
  matchfun=function(d)d$S==d$lR,
  contrasts=list(v=list(E=E2mat,lM=ADmat),B=list(E=Emat,E2=ESmat),
                 sv=list(lM=ADmat)),
  functions=list(E2=function(d){factor(d$E=="speed",labels=c("notspeed","speed"))}),
  formula=list(v~E/lM,B~E2*lR,A~1,t0~E,sv~lM),
  constants=c(sv=log(1)))

pmean = c(
  v=3,'v_Eneutral'=0,'v_Espeed'=0,'v_Eaccuracy:lMd'=2,'v_Eneutral:lMd'=2,'v_Espeed:lMd'=2,
  B=log(2),'B_Ens-s'=0,'B_E2ns-s:lRright'=0,'B_lRright'=0,
  A=log(1),t0=log(.2),'t0_Eneutral'=0,'t0_Espeed'=0,sv_lMd=0)
psd = c(2,2,2,2,2,2,1,.25,.25,.25,0.5,.5,1,1,1)
priorEavt04 <- prior(design_Eavt03,mu_mean = pmean, mu_sd = psd)
plot_prior(priorEavt04,design=design_Eavt04,layout=c(3,6))

samplers <- make_emc(forstmann,design_Eavt03,prior=priorEavt03)
# save(samplers,file="Hierarchical/fits_LBA/Eavt03.RData")


design_Eav3 <- design(data=forstmann,model=LBA,
  matchfun=function(d)d$S==d$lR,
  contrasts=list(v=list(E=E2mat,lM=ADmat),B=list(E2=ESmat),sv=list(lM=ADmat)),
  functions=list(E2=function(d){factor(d$E=="speed",labels=c("notspeed","speed"))}),
  formula=list(v~E/lM,B~E2*lR,A~1,t0~1,sv~lM),
  constants=c(sv=log(1)))

pmean = c(
  v=3,'v_Eneutral'=0,'v_Espeed'=0,'v_Eaccuracy:lMd'=2,'v_Eneutral:lMd'=2,'v_Espeed:lMd'=2,
  B=log(2),'B_Ens-s'=0,'B_lRright'=0,'B_E2ns-s:lRright'=0,
  A=log(1),t0=log(.2),sv_lMd=0)
psd = c(2,2,2,2,2,2,1,.25,.25,.25,0.5,.5,1)
priorEav3 <- prior(design_Eav3,mu_mean = pmean, mu_sd = psd)
plot_prior(priorEav3,design_Eav3,layout=c(2,8))

samplers <- make_emc(forstmann,design_Eav3,prior=priorEav3)
# save(samplers,file="Hierarchical/fits_LBA/Eav3.RData")


print(load("Hierarchical/fits_LBA/Eavt03.RData"))
print(load("Hierarchical/fits_LBA/Eav3.RData"))

# a) Eav2 still wins on MD, Eavt03 wins on both ICs
compare(list(Eavt0=sEavt0,Eavt02=sEavt02,Eavt03=sEavt03,
             Eav=sEav,    Eav2=sEav2,     Eav3=sEav3))
#            MD wMD    DIC wDIC   BPIC wBPIC EffectiveN  meanD  Dmean   minD
# Eavt0  -14879   0 -17062    0 -16669     0        393 -17456 -17717 -17849
# Eavt02 -15095   0 -17273    0 -16953     0        320 -17593 -17782 -17913
# Eavt03 -14774   0 -17340    1 -16984     1        356 -17695 -17897 -18051
# Eav    -15292   0 -17169    0 -16877     0        292 -17461 -17684 -17753
# Eav2   -15614   1 -17196    0 -16932     0        265 -17461 -17650 -17725
# Eav3   -15320   0 -17213    0 -16867     0        346 -17560 -17757 -17906

# b) Same pattern of a few very small misses of the 95% CI, fixing a few of the
#    misses in the non-bias model.
# Eav3: very small under-prediction of accuracy in speed left and 10th percentile
#       in accuracy left
# Eavt03: very small under-prediction of accuracy in speed left and 10th percentile
#       in accuracy left

# 3. 
posterior_summary(sEav3,digits=2)
posterior_summary(sEav3,digits=2,map=TRUE)

# a) Rate quantity is greater for neutral than accuracy by .3 [.17-.44] with 
#    strong BF support.
hypothesis(sEav3,"v_Eneutral")
#    and quality is greater in accuracy by .31 [.11 to .53] has positive BF
#    support
credible(sEav3,c("v_Eaccuracy:lMd","v_Eneutral:lMd"))
hypothesis(sEav3,fun=\(x)diff(x[c("v_Eaccuracy:lMd","v_Eneutral:lMd")]))

# b) Thresholds are lower for left than right in the speed condition by .05
#    [.03 to.09] with positive BF support, and also in the non-speed conditions by 
#    .04 [.001 to .079] but with an equivocal BF. The larger bias in speed
#    .014 [-.011 to .043] has strong BF support for no difference.

credible(sEav3,c("B_E2speed_lRright","B_E2speed_lRleft"),map=TRUE)
hypothesis(sEav3,fun=\(x)diff(x[c("B_E2speed_lRright","B_E2speed_lRleft")]),map=TRUE)

credible(sEav3,c("B_E2notspeed_lRright","B_E2notspeed_lRleft"),map=TRUE,digits=3)
hypothesis(sEav3,fun=\(x)diff(x[c("B_E2notspeed_lRright","B_E2notspeed_lRleft")]),map=TRUE)

fun <- \(x){diff(x[c("B_E2notspeed_lRright","B_E2notspeed_lRleft")])-
            diff(x[c("B_E2speed_lRright","B_E2speed_lRleft")])}
credible(sEav3,x_fun=fun,map=TRUE,digits=3)
hypothesis(sEav3,fun=fun,map=TRUE)

# c) Rate quantity is greater for speed than accuracy by .78 [.41 to 1.11] with 
#    strong BF support and rate quality by 1.33 [.83 to 1.8] with strong BF
#    support. Caution for speed is less than caution for non-speed by .2 
#    [.15 to .25] with strong BF support.
credible(sEav3,"v_Espeed")
hypothesis(sEav3,"v_Espeed")

credible(sEav3,c("v_Eaccuracy:lMd","v_Espeed:lMd"))
hypothesis(sEav3,fun=\(x)diff(x[c("v_Eaccuracy:lMd","v_Espeed:lMd")]))

fun <- \(x){mean(x[c("B_E2notspeed_lRright","B_E2notspeed_lRleft")])-
            mean(x[c("B_E2speed_lRright","B_E2speed_lRleft")])}
credible(sEav3,x_fun=fun,map=TRUE)
hypothesis(sEav3,fun=fun,map=TRUE)



####  Fit plots ----
ppEav2 <- predict(sEav2,n_cores=12)
ppEav3 <- predict(sEav3,n_cores=12)
ppEavt03 <- predict(sEavt03,n_cores=12)
save(ppEav2,ppEav3,ppEavt03,file="Hierarchical/fits_LBA/ppExercises.RData")

pp <- ppEav2   
pp <- ppEav3
pp <- ppEavt03

# CDFs
plot_fit(forstmann,pp,layout=c(2,3))

# Accuracy
pc <- function(d) 100*mean(d$S==d$R)
tab <- plot_fit(forstmann,pp,layout=c(2,3),factors=c("E","S"),
                stat=pc,stat_name="Accuracy (%)",xlim=c(70,95))
round(tab,2)

# Note that we could also average over stimulus if our main interest was on the
# emphasis factor. The good fit here is because the bias effect averages out.
tab <- plot_fit(forstmann,pp,layout=c(1,3),factors=c("E"),
                stat=pc,stat_name="Accuracy (%)",xlim=c(70,95))
round(tab,3)


# Speed for correct responses

# Mean RT 
tab <- plot_fit(forstmann,pp,layout=c(2,3),factors=c("E","S"),
  stat=function(d){mean(d$rt[d$R==d$S])},stat_name="Mean Correct RT (s)",xlim=c(0.375,.6))
round(tab,3)

# Fast responses (10th percentile): also good
tab <- plot_fit(forstmann,pp,layout=c(2,3),factors=c("E","S"),xlim=c(0.275,.4),
  stat=function(d){quantile(d$rt[d$R==d$S],.1)},stat_name="10th Percentile Correct (s)")
round(tab,3)

# Slow responses (90th percentile)
tab <- plot_fit(forstmann,pp,layout=c(2,3),factors=c("E","S"),xlim=c(0.525,.975),
  stat=function(d){quantile(d$rt[d$R==d$S],.9)},stat_name="90th Correct Percentile (s)")
round(tab,3)

# RT variability (SD): hence variability pretty good
tab <- plot_fit(forstmann,pp,layout=c(2,3),factors=c("E","S"),
  stat=function(d){sd(d$rt[d$R==d$S])},stat_name="SD Correct (s)",xlim=c(0.1,.325))
round(tab,3)

# Errors speed: Even error speed well accommodated.
tab <- plot_fit(forstmann,pp,layout=c(2,3),factors=c("E","S"),xlim=c(0.375,.725),
  stat=function(d){mean(d$rt[d$R!=d$S])},stat_name="Mean Error RT (s)")
round(tab,3)
