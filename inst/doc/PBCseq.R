###################################################
### chunk number 1: Options
###################################################
#line 113 "mixAKclust01.Rnw"
#OPT <- options(prompt = "R> ", continue = "+  ", width = 70, useFancyQuotes = FALSE)
OPT <- options(prompt = "R> ", continue = "+  ", width = 75, useFancyQuotes = FALSE)

figSweave <- function(){par(mar = c(5, 4, 4, 1) + 0.1)}
options(SweaveHooks = list(fig = figSweave))


###################################################
### chunk number 2: Load previously calculated results
###################################################
#line 122 "mixAKclust01.Rnw"
### This chunk loads results previously calculated and stored on the author's disk.
### The user should ignore this chunk since it is necessary to run all other commands first
### to get all the results.
if (!("mod01" %in% ls())){ 
  load("/home/komarek/RESULT_OBJ/glmmClust-20110120/PBC_Mayo02-MCMC.RData")
  mod01 <- mod02
  pdm01 <- pdm02
  plugdm01 <- plugdm02
  pdj01 <- pdj02
  plugdj01 <- plugdj02
  rm(list=c("mod02", "pdm02", "plugdm02", "pdj02", "plugdj02"))
}

if (!("devs" %in% ls())){ 
  #load("/home/komarek/RESULT_OBJ/glmmClust-20110120/PBC_Mayo04-modelComparison-devs.RData")
  load("/home/komarek/RESULT_OBJ/glmmClust-20101223/PBC_Mayo04-modelComparison-devs.RData")  
}


###################################################
### chunk number 3: load package and data
###################################################
#line 191 "mixAKclust01.Rnw"
library("mixAK")
data(PBCseq, package="mixAK")


###################################################
### chunk number 4: take only people alive at 910 days
###################################################
#line 227 "mixAKclust01.Rnw"
idTake <- subset(PBCseq, day == 0 & alive >= 910)[, "id"]
pbc01 <- subset(PBCseq, id %in% idTake & day <= 910, 
   select=c("id", "day", "month", "fu.days", "delta.ltx.death", 
            "lbili", "platelet", "spiders"))
rownames(pbc01) <- 1:nrow(pbc01)
head(pbc01)
tail(pbc01)


###################################################
### chunk number 5: summary
###################################################
#line 237 "mixAKclust01.Rnw"
length(unique(pbc01$id))
table(table(pbc01$id))
summary(pbc01)


###################################################
### chunk number 6: extract longitudinal profiles
###################################################
#line 270 "mixAKclust01.Rnw"
ip <- getProfiles(t = "month", 
   y = c("lbili", "platelet", "spiders"), id = "id", data = pbc01)
print(ip[[1]])


###################################################
### chunk number 7: Sweave options
###################################################
#line 278 "mixAKclust01.Rnw"
figSweave <- function(){par(mar = c(5, 4, 3, 0) + 0.1, bty = "n")}
options(SweaveHooks = list(fig = figSweave))


###################################################
### chunk number 8: 01-long_prof
###################################################
#line 282 "mixAKclust01.Rnw"
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
### chunk number 9: running MCMC eval=FALSE
###################################################
## #line 422 "mixAKclust01.Rnw"
## set.seed(20042007)
## mod01 <- GLMM_MCMC(y = pbc01[, c("lbili", "platelet", "spiders")], 
##     dist = c("gaussian", "poisson(log)", "binomial(logit)"), 
##     id = pbc01[, "id"],
##     x = list(lbili    = "empty", 
##              platelet = "empty", 
##              spiders  = pbc01[, "month"]),
##     z = list(lbili    = pbc01[, "month"], 
##              platelet = pbc01[, "month"], 
##              spiders  = "empty"),
##     random.intercept = rep(TRUE, 3),
##     prior.b = list(Kmax = 2), 
##     nMCMC = c(burn = 1000, keep = 10000, thin = 100, info = 1000))


###################################################
### chunk number 10: components of object of class GLMM_MCMC
###################################################
#line 471 "mixAKclust01.Rnw"
names(mod01)


###################################################
### chunk number 11: shift vector and scale matrix
###################################################
#line 480 "mixAKclust01.Rnw"
print(mod01$scale.b)


###################################################
### chunk number 12: prior for theta
###################################################
#line 504 "mixAKclust01.Rnw"
print(mod01$prior.b)


###################################################
### chunk number 13: prior for alpha
###################################################
#line 584 "mixAKclust01.Rnw"
print(mod01$prior.alpha)


###################################################
### chunk number 14: prior for sigma2
###################################################
#line 615 "mixAKclust01.Rnw"
print(mod01$prior.eps)


###################################################
### chunk number 15: initial values for random effects and mixture parameters
###################################################
#line 645 "mixAKclust01.Rnw"
print(mod01$init.b)


###################################################
### chunk number 16: initial values for fixed effects
###################################################
#line 712 "mixAKclust01.Rnw"
print(mod01$init.alpha)


###################################################
### chunk number 17: initial values for dispersion parameters
###################################################
#line 722 "mixAKclust01.Rnw"
print(mod01$init.eps)


###################################################
### chunk number 18: full specification of the prior hyperparameters by the user eval=FALSE
###################################################
## #line 756 "mixAKclust01.Rnw"
## set.seed(20042007)
## mod02 <- GLMM_MCMC(y = pbc01[, c("lbili", "platelet", "spiders")], 
##     dist = c("gaussian", "poisson(log)", "binomial(logit)"), 
##     id = pbc01[, "id"],
##     x = list(lbili    = "empty", 
##              platelet = "empty", 
##              spiders  = pbc01[, "month"]),
##     z = list(lbili    = pbc01[, "month"], 
##              platelet = pbc01[, "month"], 
##              spiders  = "empty"),
##     random.intercept = rep(TRUE, 3),
##     scale.b = list(shift = c(0.315158108, 0.007654708, 5.526209768, 
##                              -0.006634000, -2.749539170),
##                    scale = c(0.86449212, 0.02007624, 0.34860152, 
##                              0.01565365, 3.22849129)),                       
##     prior.b = list(Kmax = 2, priormuQ = "independentC", 
##                    delta = 1, xi = rep(0, 5), D = diag(rep(36, 5)), 
##                    zeta = 6, g = rep(0.2, 5), h = rep(0.2777778, 5)), 
##     prior.alpha = list(mean = 0, var = 10000),                       
##     prior.eps = list(zeta = 2, g = 0.2, h = 2.755851),
##     init.b = mod01$state.last.b,
##     init.alpha = mod01$state.last.alpha,
##     init.eps = mod01$state.last.eps,
##     nMCMC = c(burn = 0, keep = 1000, thin = 100, info = 1000))


###################################################
### chunk number 19: print part of w sample
###################################################
#line 792 "mixAKclust01.Rnw"
print(mod01$w_b[1:3,])


###################################################
### chunk number 20: print part of mu sample
###################################################
#line 798 "mixAKclust01.Rnw"
print(mod01$mu_b[1:3,])


###################################################
### chunk number 21: print part of D sample
###################################################
#line 805 "mixAKclust01.Rnw"
print(mod01$Sigma_b[1:3,])


###################################################
### chunk number 22: print part of Q and Li samples
###################################################
#line 813 "mixAKclust01.Rnw"
print(mod01$Q_b[1:3,])
print(mod01$Li_b[1:3,])


###################################################
### chunk number 23: print part of alpha sample
###################################################
#line 830 "mixAKclust01.Rnw"
print(mod01$alpha[1:3,])


###################################################
### chunk number 24: print part of sigma sample
###################################################
#line 835 "mixAKclust01.Rnw"
print(mod01$sigma_eps[1:3,])


###################################################
### chunk number 25: print part of mixture.b sample
###################################################
#line 844 "mixAKclust01.Rnw"
print(mod01$mixture_b[1:3,])


###################################################
### chunk number 26: print part of Deviance
###################################################
#line 861 "mixAKclust01.Rnw"
print(mod01$Deviance[1:10])


###################################################
### chunk number 27: running re-labelling algorithm and storing of sampled component probabilities eval=FALSE
###################################################
## #line 928 "mixAKclust01.Rnw"
## mod01  <- NMixRelabel(mod01, type = "stephens", keep.comp.prob = TRUE)


###################################################
### chunk number 28: Sweave options
###################################################
#line 974 "mixAKclust01.Rnw"
figSweave <- function(){par(mar = c(4, 4, 0, 1) + 0.1)}
options(SweaveHooks = list(fig = figSweave))


###################################################
### chunk number 29: 02-trace_deviance
###################################################
#line 978 "mixAKclust01.Rnw"
tracePlots(mod01, param = "Deviance")


###################################################
### chunk number 30: traceplot of conditional deviance and fixed effects and residual variance
###################################################
#line 985 "mixAKclust01.Rnw"
tracePlots(mod01, param = "Cond.Deviance")
tracePlots(mod01, param = "alpha")
tracePlots(mod01, param = "sigma_eps")


###################################################
### chunk number 31: Sweave options
###################################################
#line 993 "mixAKclust01.Rnw"
figSweave <- function(){par(mar = c(4, 4, 4, 1) + 0.1)}
options(SweaveHooks = list(fig = figSweave))


###################################################
### chunk number 32: 03-trace_Eb
###################################################
#line 997 "mixAKclust01.Rnw"
COL <- rep(rainbow_hcl(3, start = 30, end = 210), c(2, 2, 1))
tracePlots(mod01, param = "Eb", col = COL)


###################################################
### chunk number 33: traceplots of SDb and Corb
###################################################
#line 1003 "mixAKclust01.Rnw"
tracePlots(mod01, param = "SDb")
tracePlots(mod01, param = "Corb")


###################################################
### chunk number 34: traceplots for mixture weighs and means and standard deviations
###################################################
#line 1012 "mixAKclust01.Rnw"
tracePlots(mod01, param = "w_b")
tracePlots(mod01, param = "mu_b")
tracePlots(mod01, param = "sd_b")


###################################################
### chunk number 35: traceplots for mixture weighs and means and standard deviations after relabelling
###################################################
#line 1020 "mixAKclust01.Rnw"
tracePlots(mod01, param = "w_b", relabel = TRUE)
tracePlots(mod01, param = "mu_b", relabel = TRUE)
tracePlots(mod01, param = "sd_b", relabel = TRUE)


###################################################
### chunk number 36: traceplots of variance hyperparameters
###################################################
#line 1027 "mixAKclust01.Rnw"
tracePlots(mod01, param = "gammaInv_b")
tracePlots(mod01, param = "gammaInv_eps")


###################################################
### chunk number 37: Sweave options
###################################################
#line 1039 "mixAKclust01.Rnw"
figSweave <- function(){par(mar = c(4, 4, 1, 1) + 0.1, bty = "n")}
options(SweaveHooks = list(fig = figSweave))


###################################################
### chunk number 38: 04-acor_deviance
###################################################
#line 1043 "mixAKclust01.Rnw"
autocorr.plot(mod01$Deviance, lag.max = 20, col = "blue4", 
              auto.layout = FALSE, lwd = 2)


###################################################
### chunk number 39: autocorrelation plots of conditional deviance and fixed effects and residual standard deviation
###################################################
#line 1057 "mixAKclust01.Rnw"
autocorr.plot(mod01$Cond.Deviance, col = "blue4", auto.layout = FALSE, lwd = 2)
autocorr.plot(mod01$alpha, col = "blue4", auto.layout = FALSE, lwd = 2)
autocorr.plot(mod01$sigma_eps, col = "blue4", auto.layout = FALSE, lwd = 2)


###################################################
### chunk number 40: Sweave options
###################################################
#line 1066 "mixAKclust01.Rnw"
figSweave <- function(){par(mar = c(4, 4, 4, 1) + 0.1, bty = "n")}
options(SweaveHooks = list(fig = figSweave))


###################################################
### chunk number 41: 05-acor_Eb
###################################################
#line 1070 "mixAKclust01.Rnw"
layout(autolayout(5))
name.Eb <- paste("b.Mean.", 1:5, sep = "")
autocorr.plot(mod01$mixture_b[, name.Eb], lag.max = 20, col = "blue4", 
              auto.layout = FALSE, lwd = 2)


###################################################
### chunk number 42: autocorrelation plots of SDb and Corb
###################################################
#line 1084 "mixAKclust01.Rnw"
layout(autolayout(5))
name.SDb <- paste("b.SD.", 1:5, sep = "")
autocorr.plot(mod01$mixture_b[, name.SDb], lag.max = 20, col = "blue4", 
              auto.layout = FALSE, lwd = 2)
#
layout(autolayout(10))
name.Corb <- paste("b.Corr.", c(2:5, 3:5, 4:5, 5), ".", rep(1:4, 4:1), sep = "")
autocorr.plot(mod01$mixture_b[, name.Corb], lag.max = 20, col = "blue4", 
              auto.layout = FALSE, lwd = 2)


###################################################
### chunk number 43: posterior summary statistics for basic model parameters
###################################################
#line 1134 "mixAKclust01.Rnw"
print(mod01)


###################################################
### chunk number 44: coda posterior summary for regression parameters
###################################################
#line 1146 "mixAKclust01.Rnw"
Regr <- cbind(mod01$mixture_b[, name.Eb], 
              mod01$alpha)
colnames(Regr) <- paste(rep(c("lbili", "platelet", "spiders"), each = 2), 
                        ":", rep(c("Intcpt", "Slope"), 3), sep="")
Regr <- mcmc(Regr)
summary(Regr)


###################################################
### chunk number 45: coda posterior summary for residual standard deviation and SDb and Corb
###################################################
#line 1159 "mixAKclust01.Rnw"
sigma1 <- mcmc(mod01$sigma_eps)
summary(sigma1)
#
SDb <- mcmc(mod01$mixture_b[, name.SDb])
summary(SDb)
#
Corb <- mcmc(mod01$mixture_b[, name.Corb])
summary(Corb)


###################################################
### chunk number 46: HPD intervals for regression parameters
###################################################
#line 1183 "mixAKclust01.Rnw"
HPDinterval(Regr)


###################################################
### chunk number 47: HPD intervals for residual standard deviation and SDb and Corb
###################################################
#line 1186 "mixAKclust01.Rnw"
HPDinterval(sigma1)
HPDinterval(SDb)
HPDinterval(Corb)


###################################################
### chunk number 48: Sweave options
###################################################
#line 1198 "mixAKclust01.Rnw"
figSweave <- function(){par(mar = c(4, 4, 4, 1) + 0.1, bty = "n")}
options(SweaveHooks = list(fig = figSweave))


###################################################
### chunk number 49: 06-density_Eb
###################################################
#line 1202 "mixAKclust01.Rnw"
COL <- rep(rainbow_hcl(3, start = 30, end = 210), each = 2)
par(mfcol = c(2, 3))
for (i in 1:6){
  densplot(Regr[, i], show.obs = FALSE, col = COL[i], lwd = 2)
  title(main = colnames(Regr)[i])
}


###################################################
### chunk number 50: calculation of fitted longitudinal profiles
###################################################
#line 1236 "mixAKclust01.Rnw"
tpred <- seq(0, 30, by = 0.3)
fitMean <- fitted(mod01, x = list("empty", "empty", tpred), 
                         z = list(tpred, tpred, "empty"), 
                         statistic = "mean", overall = TRUE)
fitMed <- fitted(mod01, x = list("empty", "empty", tpred), 
                        z = list(tpred, tpred, "empty"), 
                        statistic = "median", overall = TRUE)


###################################################
### chunk number 51: print parts of calculated longitudinal profiles for platelet counts
###################################################
#line 1249 "mixAKclust01.Rnw"
print(fitMean[[2]][1:10,])


###################################################
### chunk number 52: Sweave options
###################################################
#line 1254 "mixAKclust01.Rnw"
figSweave <- function(){par(mar = c(4, 4, 4, 1) + 0.1, bty = "n")}
options(SweaveHooks = list(fig = figSweave))


###################################################
### chunk number 53: 12-fitted_profiles
###################################################
#line 1258 "mixAKclust01.Rnw"
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
### chunk number 54: posterior means of individual random effects
###################################################
#line 1293 "mixAKclust01.Rnw"
bhat <- mod01$poster.mean.profile[, 1:mod01$dimb]
print(bhat[1:10,])


###################################################
### chunk number 55: print poster.mean.y component
###################################################
#line 1314 "mixAKclust01.Rnw"
print(mod01$poster.mean.y$lbili[1:3,])
print(mod01$poster.mean.y$platelet[1:3,])
print(mod01$poster.mean.y$spiders[1:3,])


###################################################
### chunk number 56: Sweave options
###################################################
#line 1333 "mixAKclust01.Rnw"
figSweave <- function(){par(mar = c(4, 4, 4, 1) + 0.1, bty = "n")}
options(SweaveHooks = list(fig = figSweave))


###################################################
### chunk number 57: 07-resplot
###################################################
#line 1337 "mixAKclust01.Rnw"
COL <- rainbow_hcl(3, start = 30, end = 210)
MAIN <- c("Log(bilirubin)", "Platelet count", "Blood vessel malform.")
layout(autolayout(3))
for (i in 1:3){
  plot(mod01$poster.mean.y[[i]]$fitted, mod01$poster.mean.y[[i]]$stres, 
       xlab = "Fitted", ylab = "Standard. residuals", col = COL[i])
  lines(lowess(mod01$poster.mean.y[[i]]$fitted, mod01$poster.mean.y[[i]]$stres),
        col = "red4")
  title(main = MAIN[i])
}


###################################################
### chunk number 58: posterior means of mixture parameters after relabelling
###################################################
#line 1410 "mixAKclust01.Rnw"
print(mod01$poster.mean.w_b)
print(mod01$poster.mean.mu_b)


###################################################
### chunk number 59: posterior means of mixture covariance matrices
###################################################
#line 1415 "mixAKclust01.Rnw"
print(mod01$poster.mean.Sigma_b)


###################################################
### chunk number 60: posterior means of shifted and scaled mixture parameters
###################################################
#line 1427 "mixAKclust01.Rnw"
NMixSummComp(mod01)


###################################################
### chunk number 61: estimated marginal densities of random effects eval=FALSE
###################################################
## #line 1447 "mixAKclust01.Rnw"
## plugdm01 <- NMixPlugDensMarg(mod01)
## pdm01    <- NMixPredDensMarg(mod01)


###################################################
### chunk number 62: basic plots of calculated estimates of marginal densities of random effects
###################################################
#line 1452 "mixAKclust01.Rnw"
plot(plugdm01)
plot(pdm01)


###################################################
### chunk number 63: Sweave options
###################################################
#line 1461 "mixAKclust01.Rnw"
figSweave <- function(){par(mar = c(4, 4, 4, 1) + 0.1, bty = "n")}
options(SweaveHooks = list(fig = figSweave))


###################################################
### chunk number 64: 08-dens_b_marg
###################################################
#line 1465 "mixAKclust01.Rnw"
layout(autolayout(mod01$dimb))
blab <- c("Intcpt (lbili)", "Slope (lbili)", "Intcpt (platelet)", 
          "Slope (platelet)", "Intcpt (spiders)")
for (i in 1:mod01$dimb){
  plot(plugdm01$x[[i]], plugdm01$dens[[i]], type = "l", 
       xlab = "b", ylab = "Density", col = "red", main = blab[i], 
       ylim = c(0, max(pdm01$dens[[i]])))
  lines(pdm01$x[[i]], pdm01$dens[[i]], col = "darkblue")
}


###################################################
### chunk number 65: estimated marginal densities of random effects with explicit specification of grid values eval=FALSE
###################################################
## #line 1479 "mixAKclust01.Rnw"
## bgrid <- list(b1 = seq(-2.73, 3.37, length = 50),
##               b2 = seq(-0.07, 0.08, length = 50),
##               b3 = seq(4.28, 6.77, length = 50),
##               b4 = seq(-0.06, 0.05, length = 50),
##               b5 = seq(-16.05, 10.34, length = 50))
## 
## plugdm01 <- NMixPlugDensMarg(mod01, grid = bgrid)
## pdm01    <- NMixPredDensMarg(mod01, grid = bgrid)


###################################################
### chunk number 66: estimated joint bivariate marginal densities of random effects eval=FALSE
###################################################
## #line 1502 "mixAKclust01.Rnw"
## plugdj01 <- NMixPlugDensJoint2(mod01)  
## pdj01    <- NMixPredDensJoint2(mod01)
## 
## plot(plugdj01)
## plot(pdj01)


###################################################
### chunk number 67: Sweave options
###################################################
#line 1513 "mixAKclust01.Rnw"
figSweave <- function(){par(mar = c(4, 4, 2, 1) + 0.1, bty = "n")}
options(SweaveHooks = list(fig = figSweave))


###################################################
### chunk number 68: 09-dens_b_joint2
###################################################
#line 1517 "mixAKclust01.Rnw"
layout(matrix(c(1:9, 0, 10, 0), ncol = 3, byrow = TRUE))
for (i in 1:(mod01$dimb - 1)){
  for (j in (i+1):mod01$dimb){
    image(pdj01$x[[i]], pdj01$x[[j]], pdj01$dens[[paste(i, "-", j, sep = "")]], 
       col = rev(heat_hcl(33, c = c(80, 30), l = c(30, 90), power = c(1/5, 1.3))), 
       xlab = blab[i], ylab = blab[j])
    contour(pdj01$x[[i]], pdj01$x[[j]], pdj01$dens[[paste(i, "-", j, sep = "")]], 
         col = "brown", add = TRUE)
  }  
}  


###################################################
### chunk number 69: estimated joint bivariate marginal densities of random effects with explicit specification of grid values eval=FALSE
###################################################
## #line 1531 "mixAKclust01.Rnw"
## plugdj01 <- NMixPlugDensJoint2(mod01, grid = bgrid)
## pdj01    <- NMixPredDensJoint2(mod01, grid = bgrid)


###################################################
### chunk number 70: print part of poster.comp.prob3
###################################################
#line 1584 "mixAKclust01.Rnw"
print(mod01$poster.comp.prob3[1:5,])


###################################################
### chunk number 71: print part of quant.comp.prob3
###################################################
#line 1589 "mixAKclust01.Rnw"
names(mod01$quant.comp.prob3)
print(mod01$quant.comp.prob3[["50%"]][1:5,])


###################################################
### chunk number 72: print part of comp.prob3
###################################################
#line 1606 "mixAKclust01.Rnw"
print(mod01$comp.prob3[1:5, 1:10])


###################################################
### chunk number 73: Sweave options
###################################################
#line 1627 "mixAKclust01.Rnw"
figSweave <- function(){par(mar = c(4, 4, 4, 1) + 0.1, bty = "n")}
options(SweaveHooks = list(fig = figSweave))


###################################################
### chunk number 74: 11-post_dist_p
###################################################
#line 1631 "mixAKclust01.Rnw"
IDS <- unique(pbc01$id)
K <- mod01$prior.b$Kmax
N <- ncol(mod01$comp.prob3) / K
ID <- c(3, 7, 51)

par(mfrow = c(1, 3))
for (id in ID){
  i <- (1:N)[IDS == id]
  hist(mod01$comp.prob3[, (i - 1) * K + 1], xlim = c(0, 1), prob = TRUE, 
     xlab = expression(paste("P(u=1|", psi, ", ", theta, ", y)", sep = "")), 
     col = heat_hcl(12, c = c(80, 30), l = c(30, 90), power = c(1/5, 2))[12], 
     main = paste("ID", id))
}  


###################################################
### chunk number 75: HPD intervals for component probabilities
###################################################
#line 1650 "mixAKclust01.Rnw"
prob3HPD <- HPDinterval(mcmc(mod01$comp.prob3))
rownames(prob3HPD) <- paste("ID", rep(IDS, each = K), ", k = ", 1:K, ":", sep = "")
print(prob3HPD[1:6,])


###################################################
### chunk number 76: posterior mean and median and HPD interval for component probabilities of selected patients
###################################################
#line 1658 "mixAKclust01.Rnw"
Row <- (1:N)[IDS %in% ID]
Mean   <- mod01$poster.comp.prob3[Row, 1]
Median <- mod01$quant.comp.prob3[["50%"]][Row, 1]
HPD    <- prob3HPD[(Row - 1) * K + 1,] 
Pshow <- data.frame(Mean = Mean, Median = Median, 
                    HPD.lower = HPD[, 1], HPD.upper = HPD[, 2])
print(Pshow)


###################################################
### chunk number 77: classification based on posterior means and medians of component probabilities
###################################################
#line 1697 "mixAKclust01.Rnw"
groupMean <- apply(mod01$poster.comp.prob3, 1, which.max)
pMean <- apply(mod01$poster.comp.prob3, 1, max)

groupMed  <- apply(mod01$quant.comp.prob3[["50%"]], 1, which.max)
pMed <- apply(mod01$quant.comp.prob3[["50%"]], 1, max)

classif <- data.frame(id = IDS, groupMean = groupMean, pMean = pMean, 
                                groupMed = groupMed,   pMed = pMed)
print(classif[1:10,])


###################################################
### chunk number 78: patients with different classicfication based on posterior means and medians of component probabilities
###################################################
#line 1710 "mixAKclust01.Rnw"
classif[groupMean != groupMed, ]


###################################################
### chunk number 79: proportions of patients classified into two groups
###################################################
#line 1715 "mixAKclust01.Rnw"
table(groupMean)
round(prop.table(table(groupMean)) * 100, 2)

table(groupMed)
round(prop.table(table(groupMed)) * 100, 2)


###################################################
### chunk number 80: lower limits of HPD intervals for component probabilities
###################################################
#line 1728 "mixAKclust01.Rnw"
prob3HPDlower <- matrix(prob3HPD[, "lower"], ncol = 2, byrow = TRUE)
print(prob3HPDlower[1:5,])


###################################################
### chunk number 81: classification which takes HPD interval into account
###################################################
#line 1734 "mixAKclust01.Rnw"
groupHPD <- apply(prob3HPDlower, 1, which.max)
pHPD <- apply(prob3HPDlower, 1, max)
groupHPD[pHPD < 0.9] <- 3

classif$groupHPD <- groupHPD
classif$pHPD <- pHPD
print(classif[1:10,])


###################################################
### chunk number 82: proportions of classified and unclassified patients
###################################################
#line 1744 "mixAKclust01.Rnw"
table(groupHPD)
round(prop.table(table(groupHPD)) * 100, 2)


###################################################
### chunk number 83: Sweave options
###################################################
#line 1767 "mixAKclust01.Rnw"
figSweave <- function(){par(mar = c(4, 4, 2, 1) + 0.1, bty = "n")}
options(SweaveHooks = list(fig = figSweave))


###################################################
### chunk number 84: 13-dens_b_joint2_bhat
###################################################
#line 1771 "mixAKclust01.Rnw"
COL <- c("darkgreen", "red4", "lightblue")
#
layout(matrix(c(1:9, 0, 10, 0), ncol = 3, byrow = TRUE))
for (i in 1:(mod01$dimb - 1)){
  for (j in (i+1):mod01$dimb){
    image(pdj01$x[[i]], pdj01$x[[j]], pdj01$dens[[paste(i, "-", j, sep = "")]], 
       col = rev(heat_hcl(33, c = c(80, 30), l = c(30, 90), power = c(1/5, 1.3))), 
       xlab = blab[i], ylab = blab[j])
    points(bhat[, i], bhat[, j], pch = 1, col = COL[groupHPD])
  }  
}  


###################################################
### chunk number 85: calculation of group specific fitted longitudinal profiles
###################################################
#line 1814 "mixAKclust01.Rnw"
tpred <- seq(0, 30, by = 0.3)
fitGroup <- fitted(mod01, x = list("empty", "empty", tpred), 
                          z = list(tpred, tpred, "empty"), 
                          overall = FALSE)


###################################################
### chunk number 86: print part of fitGroup for platelet counts
###################################################
#line 1823 "mixAKclust01.Rnw"
print(fitGroup[[2]][1:10,])


###################################################
### chunk number 87: add grouping variables to the original data
###################################################
#line 1839 "mixAKclust01.Rnw"
TAB <- table(pbc01$id)
pbc01$groupMean <- factor(rep(groupMean, TAB))
pbc01$groupMed  <- factor(rep(groupMed, TAB))
pbc01$groupHPD  <- factor(rep(groupHPD, TAB))


###################################################
### chunk number 88: extract observed longitudinal profiles
###################################################
#line 1848 "mixAKclust01.Rnw"
ip <- getProfiles(t = "month", 
   y = c("lbili", "platelet", "spiders", "groupMean", "groupMed", "groupHPD"), 
   id = "id", data = pbc01)
print(ip[[1]])


###################################################
### chunk number 89: Sweave options
###################################################
#line 1859 "mixAKclust01.Rnw"
figSweave <- function(){par(mar = c(4, 4, 4, 1) + 0.1, bty = "n")}
options(SweaveHooks = list(fig = figSweave))


###################################################
### chunk number 90: 14-fitted_profiles_group
###################################################
#line 1863 "mixAKclust01.Rnw"
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
### chunk number 91: MCMC simulation for models with K equal to 1 and 2 and 3 and 4 eval=FALSE
###################################################
## #line 1923 "mixAKclust01.Rnw"
## devs <- list()
## for (K in 1:4){
##   cat("Calculating K = ", K, "\n========================\n", sep="")
##   
##   set.seed(20042007)
##   modK <- GLMM_MCMC(y = pbc01[, c("lbili", "platelet", "spiders")], 
##      dist = c("gaussian", "poisson(log)", "binomial(logit)"), 
##      id = pbc01[, "id"],
##      x = list(lbili    = "empty", 
##               platelet = "empty", 
##               spiders  = pbc01[, "month"]),
##      z = list(lbili    = pbc01[, "month"], 
##               platelet = pbc01[, "month"], 
##               spiders  = "empty"),
##      random.intercept = rep(TRUE, 3),
##      prior.b = list(Kmax = K), 
##      nMCMC = c(burn = 1000, keep = 10000, thin = 100, info = 1000))
## 
##   devs[[K]] <- mods[[K]]$Deviance
##   rm(list = "modK")
## }


###################################################
### chunk number 92: posterior summary statistics for the difference between the observed data deviances in models with K 2 and K 1
###################################################
#line 1948 "mixAKclust01.Rnw"
summaryDiff(devs[[2]], devs[[1]])


###################################################
### chunk number 93: posterior summary statistics for the difference between the observed data deviances in models with K 3 and K 2
###################################################
#line 1961 "mixAKclust01.Rnw"
summaryDiff(devs[[3]], devs[[2]])


###################################################
### chunk number 94: posterior summary statistics for the difference between the observed data deviances in models with K 4 and K 3 or K 2
###################################################
#line 1966 "mixAKclust01.Rnw"
summaryDiff(devs[[4]], devs[[3]])
summaryDiff(devs[[4]], devs[[2]])


###################################################
### chunk number 95: Sweave options
###################################################
#line 1974 "mixAKclust01.Rnw"
figSweave <- function(){par(mar = c(4, 4, 1, 1) + 0.1, bty = "n")}
options(SweaveHooks = list(fig = figSweave))


###################################################
### chunk number 96: 10-cdf_deviance
###################################################
#line 1978 "mixAKclust01.Rnw"
COL <- terrain_hcl(4, c = c(65, 15), l = c(45, 80), power = c(0.5, 1.5))
plot(c(14050, 14500), c(0, 1), type = "n", 
     xlab = "Deviance", ylab = "Posterior CDF")
for (K in 1:4){
  medDEV <- median(devs[[K]])
  ECDF <- ecdf(devs[[K]])
  plot(ECDF, col = COL[K], lwd = 2, add = TRUE)
  if (K <= 3) text(medDEV + 0.5, 0.5, labels = K) 
  else        text(14250, 0.23, labels = 4)
}  


###################################################
### chunk number 97: Options back to original
###################################################
#line 2012 "mixAKclust01.Rnw"
options(prompt = OPT$prompt, continue = OPT$continue, width = OPT$width, useFancyQuotes = OPT$useFancyQuotes)


