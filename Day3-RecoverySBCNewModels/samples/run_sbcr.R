library(EMC2)

matchfun=function(data)data$S==data$lR
ADmat <- cbind(d = c(-1/2,1/2))
E2 <- function(data) factor(data$E!="speed",labels=c("speed","nonspeed"))
E_incr <- contr.increasing(3)
colnames(E_incr) <- c("nonspeed", "acc")
design_LBABvE2 <- design(
  factors=list(subjects=1,E=levels(forstmann$E),S=levels(forstmann$S)),
  Rlevels=levels(forstmann$R),model=LBA,matchfun=matchfun,
  formula=list(v~E*lM,sv~lM,B~E2+lR,A~1,t0~1),
  contrasts=list(lM = ADmat, E = E_incr),
  constants=c(sv=log(1), `v_Eacc:lMd` = 0),
  functions=list(E2=E2))

E2 <- get(load("E2.RData"))
pars <- parameters(E2,selection="alpha")
pmean <- apply(pars[,-1],2,mean)
psd <- apply(pars[,-1],2,sd)

priorE2 <- prior(design_LBABvE2,mu_mean=pmean,mu_sd=psd,type="single")
priorE2[[2]] <- cov(pars[,-1])

SBC_E2r <- run_sbc(design_LBABvE2, priorE2, replicates = 200, trials = 125, plot_data = FALSE,
                  iter = 1000, n_post = 1000, fileName = "SBC_E2r.RData",
                  cores_per_chain = 25)
save(SBC_E2r,file="SBC_E2r1.RData")


