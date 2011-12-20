###################################################
### code chunk number 1: Options
###################################################
#OPT <- options(prompt = "R> ", continue = "+  ", width = 70, useFancyQuotes = FALSE)
OPT <- options(prompt = "R> ", continue = "+  ", width = 75, useFancyQuotes = FALSE)

figSweave <- function(){par(mar = c(5, 4, 4, 1) + 0.1)}
options(SweaveHooks = list(fig = figSweave))


###################################################
### code chunk number 2: Load previously calculated results
###################################################
### This chunk loads results previously calculated and stored on the author's disk.
### The user should ignore this chunk since it is necessary to run all other commands first
### to get all the results.
if (!("mod" %in% ls())){ 
  print(load("/home/komarek/RESULT_OBJ/glmmClust-20111207/PBC_Mayo02-MCMC.RData"))
  ##print(load("/home/komarek/RESULT_OBJ/glmmClust-20111207/PBC_Mayo-mod_K2.RData"))
}

if (!("devs" %in% ls()) | !("PED" %in% ls())){ 
  Devs <- list()
  for (k in 1:4){
    print(load(paste("/home/komarek/RESULT_OBJ/glmmClust-20111207/PBC_Mayo-devs_K", k, ".RData", sep="")))
    Devs[[k]] <- devs[, "Dev1"]
    
    if (k == 1){
      PED <- as.data.frame(matrix(ped, nrow=1))
      colnames(PED) <- names(ped)
    }else PED <- rbind(PED, ped)        
  }  
  devs <- Devs
  rm(list="Devs")
}


###################################################
### code chunk number 3: load package and data
###################################################
library("mixAK")
data(PBCseq, package="mixAK")


###################################################
### code chunk number 4: take only people alive at 910 days
###################################################
idTake <- subset(PBCseq, day == 0 & alive >= 910)[, "id"]
pbc01 <- subset(PBCseq, id %in% idTake & day <= 910, 
   select=c("id", "day", "month", "fu.days", "delta.ltx.death", 
            "lbili", "platelet", "spiders"))
rownames(pbc01) <- 1:nrow(pbc01)
head(pbc01)
tail(pbc01)


###################################################
### code chunk number 5: summary
###################################################
length(unique(pbc01$id))
table(table(pbc01$id))
summary(pbc01)


###################################################
### code chunk number 6: extract longitudinal profiles
###################################################
ip <- getProfiles(t = "month", 
   y = c("lbili", "platelet", "spiders"), id = "id", data = pbc01)
print(ip[[1]])


###################################################
### code chunk number 7: Sweave options
###################################################
figSweave <- function(){par(mar = c(5, 4, 3, 0) + 0.1, bty = "n")}
options(SweaveHooks = list(fig = figSweave))


###################################################
### code chunk number 8: 01-long_prof
###################################################
getOption("SweaveHooks")[["fig"]]()
COL <- rainbow_hcl(3, start = 30, end = 210)
XLIM <- c(0, 910) / (365.25 / 12)
#
layout(autolayout(3))
plotProfiles(ip = ip, data = pbc01, var = "lbili", tvar = "month", 
   xlim = XLIM, xlab = "Time (months)", col = COL[1], 
   auto.layout = FALSE, main = "Log(bilirubin)")
plotProfiles(ip = ip, data = pbc01, var = "platelet", tvar = "month", 
   xlim = XLIM, xlab = "Time (months)", col = COL[2], 
   auto.layout = FALSE, main = "Platelet count")
plotProfiles(ip = ip, data = pbc01, var = "spiders",  tvar = "month", 
   xlim = XLIM, xlab = "Time (months)", col = COL[3], 
   auto.layout = FALSE, main = "Blood vessel malform.")


###################################################
### code chunk number 9: running MCMC (eval = FALSE)
###################################################
set.seed(20042007)
mod <- GLMM_MCMC(y = pbc01[, c("lbili", "platelet", "spiders")], 
    dist = c("gaussian", "poisson(log)", "binomial(logit)"), 
    id = pbc01[, "id"],
    x = list(lbili    = "empty", 
             platelet = "empty", 
             spiders  = pbc01[, "month"]),
    z = list(lbili    = pbc01[, "month"], 
             platelet = pbc01[, "month"], 
             spiders  = "empty"),
    random.intercept = rep(TRUE, 3),
    prior.b = list(Kmax = 2), 
    nMCMC = c(burn = 1000, keep = 10000, thin = 100, info = 1000),
    PED = TRUE, parallel = FALSE)


###################################################
### code chunk number 10: components of object of class GLMM_MCMClist
###################################################
class(mod)
names(mod)


###################################################
### code chunk number 11: components of object of class GLMM_MCMC
###################################################
class(mod[[1]])
names(mod[[1]])


###################################################
### code chunk number 12: shift vector and scale matrix
###################################################
print(mod[[1]]$scale.b)
print(mod[[2]]$scale.b)


###################################################
### code chunk number 13: prior for theta
###################################################
print(mod[[1]]$prior.b)


###################################################
### code chunk number 14: prior for theta - chain 2
###################################################
print(mod[[2]]$prior.b)


###################################################
### code chunk number 15: prior for alpha
###################################################
print(mod[[1]]$prior.alpha)


###################################################
### code chunk number 16: prior for alpha - chain 2
###################################################
print(mod[[2]]$prior.alpha)


###################################################
### code chunk number 17: prior for sigma2
###################################################
print(mod[[1]]$prior.eps)


###################################################
### code chunk number 18: prior for sigma2 - chain 2
###################################################
print(mod[[2]]$prior.eps)


###################################################
### code chunk number 19: initial values for random effects and mixture parameters
###################################################
print(mod[[1]]$init.b)


###################################################
### code chunk number 20: initial values for random effects and mixture parameters
###################################################
print(mod[[2]]$init.b)


###################################################
### code chunk number 21: initial values for fixed effects
###################################################
print(mod[[1]]$init.alpha)


###################################################
### code chunk number 22: initial values for dispersion parameters
###################################################
print(mod[[1]]$init.eps)


###################################################
### code chunk number 23: full specification of the prior hyperparameters by the user (eval = FALSE)
###################################################
set.seed(20042007)
mod02 <- GLMM_MCMC(y = pbc01[, c("lbili", "platelet", "spiders")], 
    dist = c("gaussian", "poisson(log)", "binomial(logit)"), 
    id = pbc01[, "id"],
    x = list(lbili    = "empty", 
             platelet = "empty", 
             spiders  = pbc01[, "month"]),
    z = list(lbili    = pbc01[, "month"], 
             platelet = pbc01[, "month"], 
             spiders  = "empty"),
    random.intercept = rep(TRUE, 3),
    scale.b = list(shift = c(0.315158108, 0.007654708, 5.526209768, 
                             -0.006634000, -2.749539170),
                   scale = c(0.86449212, 0.02007624, 0.34860152, 
                             0.01565365, 3.22849129)),                       
    prior.b = list(Kmax = 2, priormuQ = "independentC", 
                   delta = 1, xi = rep(0, 5), D = diag(rep(36, 5)), 
                   zeta = 6, g = rep(0.2, 5), h = rep(0.2777778, 5)), 
    prior.alpha = list(mean = 0, var = 10000),                       
    prior.eps = list(zeta = 2, g = 0.2, h = 2.755851),
    init.b  = mod[[1]]$state.last.b,
    init2.b = mod[[2]]$state.last.b,                   
    init.alpha  = mod[[1]]$state.last.alpha,
    init2.alpha = mod[[2]]$state.last.alpha,                   
    init.eps = mod[[1]]$state.last.eps,
    init2.eps = mod[[2]]$state.last.eps,                   
    nMCMC = c(burn = 0, keep = 1000, thin = 100, info = 1000),
    PED = TRUE, parallel = FALSE)


###################################################
### code chunk number 24: print part of w sample
###################################################
print(mod[[1]]$w_b[1:3,])


###################################################
### code chunk number 25: print part of mu sample
###################################################
print(mod[[1]]$mu_b[1:3,])


###################################################
### code chunk number 26: print part of D sample
###################################################
print(mod[[1]]$Sigma_b[1:3,])


###################################################
### code chunk number 27: print part of Q and Li samples
###################################################
print(mod[[1]]$Q_b[1:3,])
print(mod[[1]]$Li_b[1:3,])


###################################################
### code chunk number 28: print part of alpha sample
###################################################
print(mod[[1]]$alpha[1:3,])


###################################################
### code chunk number 29: print part of sigma sample
###################################################
print(mod[[1]]$sigma_eps[1:3,])


###################################################
### code chunk number 30: print part of mixture.b sample
###################################################
print(mod[[1]]$mixture_b[1:3,])


###################################################
### code chunk number 31: print part of Deviance
###################################################
print(mod[[1]]$Deviance[1:10])


###################################################
### code chunk number 32: running re-labelling algorithm and storing of sampled component probabilities (eval = FALSE)
###################################################
mod[[1]]  <- NMixRelabel(mod[[1]], type = "stephens", keep.comp.prob = TRUE)
mod[[2]]  <- NMixRelabel(mod[[2]], type = "stephens", keep.comp.prob = TRUE)


###################################################
### code chunk number 33: Sweave options
###################################################
figSweave <- function(){par(mar = c(4, 4, 0, 1) + 0.1)}
options(SweaveHooks = list(fig = figSweave))


###################################################
### code chunk number 34: 02-trace_deviance
###################################################
getOption("SweaveHooks")[["fig"]]()
tracePlots(mod[[1]], param = "Deviance")


###################################################
### code chunk number 35: traceplot of conditional deviance and fixed effects and residual variance
###################################################
tracePlots(mod[[1]], param = "Cond.Deviance")
tracePlots(mod[[1]], param = "alpha")
tracePlots(mod[[1]], param = "sigma_eps")


###################################################
### code chunk number 36: Sweave options
###################################################
figSweave <- function(){par(mar = c(4, 4, 4, 1) + 0.1)}
options(SweaveHooks = list(fig = figSweave))


###################################################
### code chunk number 37: 03-trace_Eb
###################################################
getOption("SweaveHooks")[["fig"]]()
COL <- rep(rainbow_hcl(3, start = 30, end = 210), c(2, 2, 1))
tracePlots(mod[[1]], param = "Eb", col = COL)


###################################################
### code chunk number 38: traceplots of SDb and Corb
###################################################
tracePlots(mod[[1]], param = "SDb")
tracePlots(mod[[1]], param = "Corb")


###################################################
### code chunk number 39: traceplots for mixture weighs and means and standard deviations
###################################################
tracePlots(mod[[1]], param = "w_b")
tracePlots(mod[[1]], param = "mu_b")
tracePlots(mod[[1]], param = "sd_b")


###################################################
### code chunk number 40: traceplots for mixture weighs and means and standard deviations after relabelling
###################################################
tracePlots(mod[[1]], param = "w_b", relabel = TRUE)
tracePlots(mod[[1]], param = "mu_b", relabel = TRUE)
tracePlots(mod[[1]], param = "sd_b", relabel = TRUE)


###################################################
### code chunk number 41: traceplots of variance hyperparameters
###################################################
tracePlots(mod[[1]], param = "gammaInv_b")
tracePlots(mod[[1]], param = "gammaInv_eps")


###################################################
### code chunk number 42: Sweave options
###################################################
figSweave <- function(){par(mar = c(4, 4, 1, 1) + 0.1, bty = "n")}
options(SweaveHooks = list(fig = figSweave))


###################################################
### code chunk number 43: 04-acor_deviance
###################################################
getOption("SweaveHooks")[["fig"]]()
autocorr.plot(mod[[1]]$Deviance, lag.max = 20, col = "blue4", 
              auto.layout = FALSE, lwd = 2)


###################################################
### code chunk number 44: autocorrelation plots of conditional deviance and fixed effects and residual standard deviation
###################################################
autocorr.plot(mod[[1]]$Cond.Deviance, col = "blue4", auto.layout = FALSE, lwd = 2)
autocorr.plot(mod[[1]]$alpha, col = "blue4", auto.layout = FALSE, lwd = 2)
autocorr.plot(mod[[1]]$sigma_eps, col = "blue4", auto.layout = FALSE, lwd = 2)


###################################################
### code chunk number 45: Sweave options
###################################################
figSweave <- function(){par(mar = c(4, 4, 4, 1) + 0.1, bty = "n")}
options(SweaveHooks = list(fig = figSweave))


###################################################
### code chunk number 46: 05-acor_Eb
###################################################
getOption("SweaveHooks")[["fig"]]()
layout(autolayout(5))
name.Eb <- paste("b.Mean.", 1:5, sep = "")
autocorr.plot(mod[[1]]$mixture_b[, name.Eb], lag.max = 20, col = "blue4", 
              auto.layout = FALSE, lwd = 2)


###################################################
### code chunk number 47: autocorrelation plots of SDb and Corb
###################################################
layout(autolayout(5))
name.SDb <- paste("b.SD.", 1:5, sep = "")
autocorr.plot(mod[[1]]$mixture_b[, name.SDb], lag.max = 20, col = "blue4", 
              auto.layout = FALSE, lwd = 2)
#
layout(autolayout(10))
name.Corb <- paste("b.Corr.", c(2:5, 3:5, 4:5, 5), ".", rep(1:4, 4:1), sep = "")
autocorr.plot(mod[[1]]$mixture_b[, name.Corb], lag.max = 20, col = "blue4", 
              auto.layout = FALSE, lwd = 2)


###################################################
### code chunk number 48: posterior summary statistics for basic model parameters
###################################################
print(mod)


###################################################
### code chunk number 49: coda posterior summary for regression parameters
###################################################
Regr <- cbind(mod[[1]]$mixture_b[, name.Eb], 
              mod[[1]]$alpha)
colnames(Regr) <- paste(rep(c("lbili", "platelet", "spiders"), each = 2), 
                        ":", rep(c("Intcpt", "Slope"), 3), sep="")
Regr <- mcmc(Regr)
summary(Regr)


###################################################
### code chunk number 50: coda posterior summary for residual standard deviation and SDb and Corb
###################################################
sigma1 <- mcmc(mod[[1]]$sigma_eps)
summary(sigma1)
#
SDb <- mcmc(mod[[1]]$mixture_b[, name.SDb])
summary(SDb)
#
Corb <- mcmc(mod[[1]]$mixture_b[, name.Corb])
summary(Corb)


###################################################
### code chunk number 51: HPD intervals for regression parameters
###################################################
HPDinterval(Regr)


###################################################
### code chunk number 52: HPD intervals for residual standard deviation and SDb and Corb
###################################################
HPDinterval(sigma1)
HPDinterval(SDb)
HPDinterval(Corb)


###################################################
### code chunk number 53: Sweave options
###################################################
figSweave <- function(){par(mar = c(4, 4, 4, 1) + 0.1, bty = "n")}
options(SweaveHooks = list(fig = figSweave))


###################################################
### code chunk number 54: 06-density_Eb
###################################################
getOption("SweaveHooks")[["fig"]]()
COL <- rep(rainbow_hcl(3, start = 30, end = 210), each = 2)
par(mfcol = c(2, 3))
for (i in 1:6){
  densplot(Regr[, i], show.obs = FALSE, col = COL[i], lwd = 2)
  title(main = colnames(Regr)[i])
}


###################################################
### code chunk number 55: calculation of fitted longitudinal profiles
###################################################
tpred <- seq(0, 30, by = 0.3)
fitMean <- fitted(mod[[1]], x = list("empty", "empty", tpred), 
                            z = list(tpred, tpred, "empty"), 
                            statistic = "mean", overall = TRUE)
fitMed <- fitted(mod[[1]], x = list("empty", "empty", tpred), 
                           z = list(tpred, tpred, "empty"), 
                           statistic = "median", overall = TRUE)


###################################################
### code chunk number 56: print parts of calculated longitudinal profiles for platelet counts
###################################################
print(fitMean[[2]][1:10,])


###################################################
### code chunk number 57: Sweave options
###################################################
figSweave <- function(){par(mar = c(4, 4, 4, 1) + 0.1, bty = "n")}
options(SweaveHooks = list(fig = figSweave))


###################################################
### code chunk number 58: 12-fitted_profiles
###################################################
getOption("SweaveHooks")[["fig"]]()
COL <- sequential_hcl(12, power = 2.2)[7]
COLmean <- "blue"
COLmed <- "red"
XLIM <- c(0, 910) / (365.25 / 12)
layout(autolayout(3))
#
plotProfiles(ip = ip, data = pbc01, var = "lbili", tvar = "month", 
   xlim = XLIM, xlab = "Time (months)", col = COL, 
   auto.layout = FALSE, main = "Log(bilirubin)")
lines(tpred, fitMean[[1]][,1], col = COLmean, lwd = 2)
lines(tpred, fitMed[[1]][,1], col = COLmed, lwd = 2)
#
plotProfiles(ip = ip, data = pbc01, var = "platelet", tvar = "month", 
   xlim = XLIM, xlab = "Time (months)", col = COL, 
   auto.layout = FALSE, main = "Platelet count")
lines(tpred, fitMean[[2]][,1], col = COLmean, lwd = 2)
lines(tpred, fitMed[[2]][,1], col = COLmed, lwd = 2)
#
plotProfiles(ip = ip, data = pbc01, var = "spiders",  tvar = "month", 
   xlim = XLIM, xlab = "Time (months)", col = COL, 
   auto.layout = FALSE, main = "Blood vessel malform.")
lines(tpred, fitMean[[3]][,1], col = COLmean, lwd = 2)
lines(tpred, fitMed[[3]][,1], col = COLmed, lwd = 2)


###################################################
### code chunk number 59: posterior means of individual random effects
###################################################
bhat <- mod[[1]]$poster.mean.profile[, 1:mod[[1]]$dimb]
print(bhat[1:10,])


###################################################
### code chunk number 60: print poster.mean.y component
###################################################
print(mod[[1]]$poster.mean.y$lbili[1:3,])
print(mod[[1]]$poster.mean.y$platelet[1:3,])
print(mod[[1]]$poster.mean.y$spiders[1:3,])


###################################################
### code chunk number 61: Sweave options
###################################################
figSweave <- function(){par(mar = c(4, 4, 4, 1) + 0.1, bty = "n")}
options(SweaveHooks = list(fig = figSweave))


###################################################
### code chunk number 62: 07-resplot
###################################################
getOption("SweaveHooks")[["fig"]]()
COL <- rainbow_hcl(3, start = 30, end = 210)
MAIN <- c("Log(bilirubin)", "Platelet count", "Blood vessel malform.")
layout(autolayout(3))
for (i in 1:3){
  plot(mod[[1]]$poster.mean.y[[i]]$fitted, mod[[1]]$poster.mean.y[[i]]$stres, 
       xlab = "Fitted", ylab = "Standard. residuals", col = COL[i])
  lines(lowess(mod[[1]]$poster.mean.y[[i]]$fitted, mod[[1]]$poster.mean.y[[i]]$stres),
        col = "red4")
  title(main = MAIN[i])
}


###################################################
### code chunk number 63: posterior means of mixture parameters after relabelling
###################################################
print(mod[[1]]$poster.mean.w_b)
print(mod[[1]]$poster.mean.mu_b)


###################################################
### code chunk number 64: posterior means of mixture covariance matrices
###################################################
print(mod[[1]]$poster.mean.Sigma_b)


###################################################
### code chunk number 65: posterior means of shifted and scaled mixture parameters
###################################################
NMixSummComp(mod[[1]])


###################################################
### code chunk number 66: estimated marginal densities of random effects (eval = FALSE)
###################################################
plugdm <- list()
plugdm[[1]] <- NMixPlugDensMarg(mod[[1]])
plugdm[[2]] <- NMixPlugDensMarg(mod[[2]])

pdm <- list()
pdm[[1]] <- NMixPredDensMarg(mod[[1]])
pdm[[2]] <- NMixPredDensMarg(mod[[2]])


###################################################
### code chunk number 67: basic plots of calculated estimates of marginal densities of random effects
###################################################
plot(plugdm[[1]])
plot(plugdm[[2]])

plot(pdm[[1]])
plot(pdm[[2]])


###################################################
### code chunk number 68: Sweave options
###################################################
figSweave <- function(){par(mar = c(4, 4, 4, 1) + 0.1, bty = "n")}
options(SweaveHooks = list(fig = figSweave))


###################################################
### code chunk number 69: 08-dens_b_marg
###################################################
getOption("SweaveHooks")[["fig"]]()
layout(autolayout(mod[[1]]$dimb))
blab <- c("Intcpt (lbili)", "Slope (lbili)", "Intcpt (platelet)", 
          "Slope (platelet)", "Intcpt (spiders)")
for (i in 1:mod[[1]]$dimb){
  plot(plugdm[[1]]$x[[i]], plugdm[[1]]$dens[[i]], type = "l", 
       xlab = "b", ylab = "Density", col = "red", main = blab[i], 
       ylim = c(0, max(pdm[[1]]$dens[[i]])))
  lines(pdm[[1]]$x[[i]], pdm[[1]]$dens[[i]], col = "darkblue")
}


###################################################
### code chunk number 70: estimated marginal densities of random effects with explicit specification of grid values (eval = FALSE)
###################################################
bgrid <- list(b1 = seq(-2.73, 3.37, length = 50),
              b2 = seq(-0.07, 0.08, length = 50),
              b3 = seq(4.28, 6.77, length = 50),
              b4 = seq(-0.06, 0.05, length = 50),
              b5 = seq(-16.05, 10.34, length = 50))

plugdm <- list()
plugdm[[1]] <- NMixPlugDensMarg(mod[[1]], grid = bgrid)
plugdm[[2]] <- NMixPlugDensMarg(mod[[2]], grid = bgrid)

pdm <- list()
pdm[[1]] <- NMixPredDensMarg(mod[[1]], grid = bgrid)
pdm[[2]] <- NMixPredDensMarg(mod[[2]], grid = bgrid)


###################################################
### code chunk number 71: estimated joint bivariate marginal densities of random effects (eval = FALSE)
###################################################
plugdj <- list()
plugdj[[1]] <- NMixPlugDensJoint2(mod[[1]])  
plugdj[[2]] <- NMixPlugDensJoint2(mod[[2]])  

pdj <- list()
pdj[[1]] <- NMixPredDensJoint2(mod[[1]])
pdj[[2]] <- NMixPredDensJoint2(mod[[2]])

plot(plugdj[[1]])
plot(plugdj[[2]])

plot(pdj[[1]])
plot(pdj[[2]])


###################################################
### code chunk number 72: Sweave options
###################################################
figSweave <- function(){par(mar = c(4, 4, 2, 1) + 0.1, bty = "n")}
options(SweaveHooks = list(fig = figSweave))


###################################################
### code chunk number 73: 09-dens_b_joint2
###################################################
getOption("SweaveHooks")[["fig"]]()
layout(matrix(c(1:9, 0, 10, 0), ncol = 3, byrow = TRUE))
for (i in 1:(mod[[1]]$dimb - 1)){
  for (j in (i+1):mod[[1]]$dimb){
    image(pdj[[1]]$x[[i]], pdj[[1]]$x[[j]], pdj[[1]]$dens[[paste(i, "-", j, sep = "")]], 
       col = rev(heat_hcl(33, c = c(80, 30), l = c(30, 90), power = c(1/5, 1.3))), 
       xlab = blab[i], ylab = blab[j])
    contour(pdj[[1]]$x[[i]], pdj[[1]]$x[[j]], pdj[[1]]$dens[[paste(i, "-", j, sep = "")]], 
         col = "brown", add = TRUE)
  }  
}  


###################################################
### code chunk number 74: estimated joint bivariate marginal densities of random effects with explicit specification of grid values (eval = FALSE)
###################################################
pdj[[1]]    <- NMixPredDensJoint2(mod[[1]], grid = bgrid)
pdj[[2]]    <- NMixPredDensJoint2(mod[[2]], grid = bgrid)

plugdj[[1]] <- NMixPlugDensJoint2(mod[[1]], grid = bgrid)
plugdj[[2]] <- NMixPlugDensJoint2(mod[[2]], grid = bgrid)


###################################################
### code chunk number 75: print part of poster.comp.prob3
###################################################
print(mod[[1]]$poster.comp.prob3[1:5,])


###################################################
### code chunk number 76: print part of quant.comp.prob3
###################################################
names(mod[[1]]$quant.comp.prob3)
print(mod[[1]]$quant.comp.prob3[["50%"]][1:5,])


###################################################
### code chunk number 77: print part of comp.prob3
###################################################
print(mod[[1]]$comp.prob3[1:5, 1:10])


###################################################
### code chunk number 78: Sweave options
###################################################
figSweave <- function(){par(mar = c(4, 4, 4, 1) + 0.1, bty = "n")}
options(SweaveHooks = list(fig = figSweave))


###################################################
### code chunk number 79: 11-post_dist_p
###################################################
getOption("SweaveHooks")[["fig"]]()
IDS <- unique(pbc01$id)
K <- mod[[1]]$prior.b$Kmax
N <- ncol(mod[[1]]$comp.prob3) / K
ID <- c(3, 7, 51)

par(mfrow = c(1, 3))
for (id in ID){
  i <- (1:N)[IDS == id]
  hist(mod[[1]]$comp.prob3[, (i - 1) * K + 1], xlim = c(0, 1), prob = TRUE, 
     xlab = expression(paste("P(u=1|", psi, ", ", theta, ", y)", sep = "")), 
     col = heat_hcl(12, c = c(80, 30), l = c(30, 90), power = c(1/5, 2))[12], 
     main = paste("ID", id))
}  


###################################################
### code chunk number 80: HPD intervals for component probabilities
###################################################
prob3HPD <- HPDinterval(mcmc(mod[[1]]$comp.prob3))
rownames(prob3HPD) <- paste("ID", rep(IDS, each = K), ", k = ", 1:K, ":", sep = "")
print(prob3HPD[1:6,])


###################################################
### code chunk number 81: posterior mean and median and HPD interval for component probabilities of selected patients
###################################################
Row <- (1:N)[IDS %in% ID]
Mean   <- mod[[1]]$poster.comp.prob3[Row, 1]
Median <- mod[[1]]$quant.comp.prob3[["50%"]][Row, 1]
HPD    <- prob3HPD[(Row - 1) * K + 1,] 
Pshow <- data.frame(Mean = Mean, Median = Median, 
                    HPD.lower = HPD[, 1], HPD.upper = HPD[, 2])
print(Pshow)


###################################################
### code chunk number 82: classification based on posterior means and medians of component probabilities
###################################################
groupMean <- apply(mod[[1]]$poster.comp.prob3, 1, which.max)
pMean <- apply(mod[[1]]$poster.comp.prob3, 1, max)

groupMed  <- apply(mod[[1]]$quant.comp.prob3[["50%"]], 1, which.max)
pMed <- apply(mod[[1]]$quant.comp.prob3[["50%"]], 1, max)

classif <- data.frame(id = IDS, groupMean = groupMean, pMean = pMean, 
                                groupMed = groupMed,   pMed = pMed)
print(classif[1:10,])


###################################################
### code chunk number 83: patients with different classicfication based on posterior means and medians of component probabilities
###################################################
classif[groupMean != groupMed, ]


###################################################
### code chunk number 84: proportions of patients classified into two groups
###################################################
table(groupMean)
round(prop.table(table(groupMean)) * 100, 2)

table(groupMed)
round(prop.table(table(groupMed)) * 100, 2)


###################################################
### code chunk number 85: lower limits of HPD intervals for component probabilities
###################################################
prob3HPDlower <- matrix(prob3HPD[, "lower"], ncol = 2, byrow = TRUE)
print(prob3HPDlower[1:5,])


###################################################
### code chunk number 86: classification which takes HPD interval into account
###################################################
groupHPD <- apply(prob3HPDlower, 1, which.max)
pHPD <- apply(prob3HPDlower, 1, max)
groupHPD[pHPD < 0.9] <- 3

classif$groupHPD <- groupHPD
classif$pHPD <- pHPD
print(classif[1:10,])


###################################################
### code chunk number 87: proportions of classified and unclassified patients
###################################################
table(groupHPD)
round(prop.table(table(groupHPD)) * 100, 2)


###################################################
### code chunk number 88: Sweave options
###################################################
figSweave <- function(){par(mar = c(4, 4, 2, 1) + 0.1, bty = "n")}
options(SweaveHooks = list(fig = figSweave))


###################################################
### code chunk number 89: 13-dens_b_joint2_bhat
###################################################
getOption("SweaveHooks")[["fig"]]()
COL <- c("darkgreen", "red4", "lightblue")
#
layout(matrix(c(1:9, 0, 10, 0), ncol = 3, byrow = TRUE))
for (i in 1:(mod[[1]]$dimb - 1)){
  for (j in (i+1):mod[[1]]$dimb){
    image(pdj[[1]]$x[[i]], pdj[[1]]$x[[j]], pdj[[1]]$dens[[paste(i, "-", j, sep = "")]], 
       col = rev(heat_hcl(33, c = c(80, 30), l = c(30, 90), power = c(1/5, 1.3))), 
       xlab = blab[i], ylab = blab[j])
    points(bhat[, i], bhat[, j], pch = 1, col = COL[groupHPD])
  }  
}  


###################################################
### code chunk number 90: calculation of group specific fitted longitudinal profiles
###################################################
tpred <- seq(0, 30, by = 0.3)
fitGroup <- fitted(mod[[1]], x = list("empty", "empty", tpred), 
                             z = list(tpred, tpred, "empty"), 
                             overall = FALSE)


###################################################
### code chunk number 91: print part of fitGroup for platelet counts
###################################################
print(fitGroup[[2]][1:10,])


###################################################
### code chunk number 92: add grouping variables to the original data
###################################################
TAB <- table(pbc01$id)
pbc01$groupMean <- factor(rep(groupMean, TAB))
pbc01$groupMed  <- factor(rep(groupMed, TAB))
pbc01$groupHPD  <- factor(rep(groupHPD, TAB))


###################################################
### code chunk number 93: extract observed longitudinal profiles
###################################################
ip <- getProfiles(t = "month", 
   y = c("lbili", "platelet", "spiders", "groupMean", "groupMed", "groupHPD"), 
   id = "id", data = pbc01)
print(ip[[1]])


###################################################
### code chunk number 94: Sweave options
###################################################
figSweave <- function(){par(mar = c(4, 4, 4, 1) + 0.1, bty = "n")}
options(SweaveHooks = list(fig = figSweave))


###################################################
### code chunk number 95: 14-fitted_profiles_group
###################################################
getOption("SweaveHooks")[["fig"]]()
#COL <- rainbow_hcl(2, start = 85, end = 40)
COL <- terrain_hcl(12, c = c(65, 0), l = c(45, 90), power = c(0.5, 1.5))[c(5, 9)]
names(COL) <- levels(pbc01$groupMed)

fitCOL <- c("darkgreen", "red4")
XLIM <- c(0, 910) / (365.25 / 12)
layout(autolayout(3))
#
plotProfiles(ip = ip, data = pbc01, var = "lbili", 
   gvar = "groupMed", tvar = "month", 
   xlim = XLIM, xlab = "Time (months)", col = COL, 
   auto.layout = FALSE, main = "Log(bilirubin)")
lines(tpred, fitGroup[[1]][,1], col = fitCOL[1], lwd = 2)
lines(tpred, fitGroup[[1]][,2], col = fitCOL[2], lwd = 2)
#
plotProfiles(ip = ip, data = pbc01, var = "platelet", 
   gvar = "groupMed", tvar = "month", 
   xlim = XLIM, xlab = "Time (months)", col = COL, 
   auto.layout = FALSE, main = "Platelet count")
lines(tpred, fitGroup[[2]][,1], col = fitCOL[1], lwd = 2)
lines(tpred, fitGroup[[2]][,2], col = fitCOL[2], lwd = 2)
#
plotProfiles(ip = ip, data = pbc01, var = "spiders",  
   gvar = "groupMed", tvar = "month", 
   xlim = XLIM, xlab = "Time (months)", col = COL, 
   auto.layout = FALSE, main = "Blood vessel malform.")
lines(tpred, fitGroup[[3]][,1], col = fitCOL[1], lwd = 2)
lines(tpred, fitGroup[[3]][,2], col = fitCOL[2], lwd = 2)


###################################################
### code chunk number 96: MCMC simulation for models with K equal to 1 and 2 and 3 and 4 (eval = FALSE)
###################################################
devs <- mods <- list()
for (K in 1:4){
  cat("Calculating K = ", K, "\n========================\n", sep="")

  if (K == 2){
    mods[[K]] <- mod
  }else{    
    set.seed(20042007)
    mods[[K]] <- GLMM_MCMC(y = pbc01[, c("lbili", "platelet", "spiders")], 
       dist = c("gaussian", "poisson(log)", "binomial(logit)"), 
       id = pbc01[, "id"],
       x = list(lbili    = "empty", 
               platelet = "empty", 
               spiders  = pbc01[, "month"]),
       z = list(lbili    = pbc01[, "month"], 
                platelet = pbc01[, "month"], 
                spiders  = "empty"),
       random.intercept = rep(TRUE, 3),
       prior.b = list(Kmax = K), 
       nMCMC = c(burn = 1000, keep = 10000, thin = 100, info = 1000),
       PED = TRUE, parallel = FALSE)
  }  

  devs[[K]] <- mods[[K]]$Deviance1    ### deviance from the first chain
  
  if (K == 1){
    PED <- as.data.frame(matrix(mods[[K]]$PED, nrow=1))
    colnames(PED) <- names(mods[[K]]$PED)
  }else PED <- rbind(PED, mods[[K]]$PED)            
}


###################################################
### code chunk number 97: PED
###################################################
print(PED)


###################################################
### code chunk number 98: posterior summary statistics for the difference between the observed data deviances in models with K 2 and K 1
###################################################
summaryDiff(devs[[2]], devs[[1]])


###################################################
### code chunk number 99: posterior summary statistics for the difference between the observed data deviances in models with K 3 and K 2
###################################################
summaryDiff(devs[[3]], devs[[2]])


###################################################
### code chunk number 100: posterior summary statistics for the difference between the observed data deviances in models with K 4 and K 3 or K 2
###################################################
summaryDiff(devs[[4]], devs[[3]])
summaryDiff(devs[[4]], devs[[2]])


###################################################
### code chunk number 101: Sweave options
###################################################
figSweave <- function(){par(mar = c(4, 4, 1, 1) + 0.1, bty = "n")}
options(SweaveHooks = list(fig = figSweave))


###################################################
### code chunk number 102: 10-cdf_deviance
###################################################
getOption("SweaveHooks")[["fig"]]()
COL <- terrain_hcl(4, c = c(65, 15), l = c(45, 80), power = c(0.5, 1.5))
plot(c(14000, 14275), c(0, 1), type="n", 
     xlab="Deviance", ylab="Posterior CDF")
for (K in 1:4){
  medDEV <- median(devs[[K]])
  ECDF <- ecdf(devs[[K]])
  plot(ECDF, col=COL[K], lwd=2, add=TRUE)
  text(medDEV+0.5, 0.5, labels=K)
}  


###################################################
### code chunk number 103: Options back to original
###################################################
options(prompt = OPT$prompt, continue = OPT$continue, width = OPT$width, useFancyQuotes = OPT$useFancyQuotes)
