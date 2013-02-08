###
###  Full pdf document describing the code included here is available at
###  http://msekce.karlin.mff.cuni.cz/~komarek/software/mixAK/PBCseq.pdf
###
### ==============================================================================


###################################################
### code chunk number 1: Options
###################################################
OPT <- options(width = 135, digits = 3)


###################################################
### code chunk number 5: load package and data
###################################################
library("mixAK")
data("PBCseq", package="mixAK")


###################################################
### code chunk number 6: take only people alive at 910 days
###################################################
idTake <- subset(PBCseq, day == 0 & alive >= 910)[, "id"]
PBC910 <- subset(PBCseq, id %in% idTake & day <= 910, 
   select=c("id", "day", "month", "fu.days", "delta.ltx.death", 
            "lbili", "platelet", "spiders"))
rownames(PBC910) <- 1:nrow(PBC910)
head(PBC910)
tail(PBC910)


###################################################
### code chunk number 7: summary
###################################################
length(unique(PBC910$id))
table(table(PBC910$id))
summary(PBC910)


###################################################
### code chunk number 8: create jittered spiders
###################################################
set.seed(20111229)
spider0 <- !is.na(PBC910[, "spiders"]) & PBC910[, "spiders"] == 0
spider1 <- !is.na(PBC910[, "spiders"]) & PBC910[, "spiders"] == 1
PBC910[, "jspiders"] <- NA
PBC910[spider0, "jspiders"] <- runif(sum(spider0), 0, 0.3)
PBC910[spider1, "jspiders"] <- runif(sum(spider1), 0.7, 1)


###################################################
### code chunk number 9: extract longitudinal profiles
###################################################
ip <- getProfiles(t = "month", 
   y = c("lbili", "platelet", "spiders", "jspiders"), id = "id", data = PBC910)
print(ip[[1]])


###################################################
### code chunk number 11: 01-long_prof
###################################################
COL <- rainbow_hcl(3, start = 30, end = 210)
XLIM <- c(0, 910) / (365.25 / 12)
#
par(mar = c(5, 4, 3, 0) + 0.1, bty = "n")
layout(autolayout(3))
plotProfiles(ip = ip, data = PBC910, var = "lbili", tvar = "month", 
   xlim = XLIM, col = COL[1], 
   xlab = "Time (months)", ylab = "Log(bilirubin)",
   auto.layout = FALSE, main = "Log(bilirubin)")
plotProfiles(ip = ip, data = PBC910, var = "platelet", tvar = "month", 
   xlim = XLIM, col = COL[2], 
   xlab = "Time (months)", ylab = "Platelet count",             
   auto.layout = FALSE, main = "Platelet count")
plotProfiles(ip = ip, data = PBC910, var = "jspiders",  tvar = "month", 
   xlim = XLIM, col = COL[3], 
   xlab = "Time (months)", ylab = "Blood vessel malform. (jittered)",
   auto.layout = FALSE, main = "Blood vessel malform.")


###################################################
### code chunk number 12: running MCMC (eval = FALSE)
###################################################
set.seed(20042007)
mod <- GLMM_MCMC(y = PBC910[, c("lbili", "platelet", "spiders")], 
    dist = c("gaussian", "poisson(log)", "binomial(logit)"), 
    id = PBC910[, "id"],
    x = list(lbili    = "empty", 
             platelet = "empty", 
             spiders  = PBC910[, "month"]),
    z = list(lbili    = PBC910[, "month"], 
             platelet = PBC910[, "month"], 
             spiders  = "empty"),
    random.intercept = rep(TRUE, 3),
    prior.b = list(Kmax = 2), 
    nMCMC = c(burn = 100, keep = 1000, thin = 10, info = 100),
    parallel = FALSE)


###################################################
### code chunk number 13: components of object of class GLMM_MCMClist
###################################################
class(mod)
names(mod)


###################################################
### code chunk number 14: components of object of class GLMM_MCMC
###################################################
class(mod[[1]])
names(mod[[1]])


###################################################
### code chunk number 15: shift vector and scale matrix for chains 1 and 2
###################################################
print(mod[[1]]$scale.b)
print(mod[[2]]$scale.b)


###################################################
### code chunk number 17: prior for theta for chains 1 and 2
###################################################
print(mod[[1]]$prior.b)
print(mod[[2]]$prior.b)


###################################################
### code chunk number 19: prior for alpha for chains 1 and 2
###################################################
print(mod[[1]]$prior.alpha)
print(mod[[2]]$prior.alpha)


###################################################
### code chunk number 21: prior for dispersion parameter for chains 1 and 2
###################################################
print(mod[[1]]$prior.eps)
print(mod[[2]]$prior.eps)


###################################################
### code chunk number 23: initial values for random effects and mixture parameters - chains 1 and 2
###################################################
print(mod[[1]]$init.b)
print(mod[[2]]$init.b)


###################################################
### code chunk number 25: initial values for random effects - chain 2
###################################################
print(mod[[2]]$init.b$b[1:5,])


###################################################
### code chunk number 26: initial values for shifted-scaled mixture means - chain 2
###################################################
print(mod[[2]]$init.b$mu)


###################################################
### code chunk number 27: initial values for shifted-scaled mixture covariance matrices - chain 2
###################################################
print(mod[[2]]$init.b$Sigma)


###################################################
### code chunk number 28: initial values for the hyperparameter gamma.b - chain 2
###################################################
print(mod[[2]]$init.b$gammaInv)


###################################################
### code chunk number 29: initial values for fixed effects - both chains
###################################################
print(mod[[1]]$init.alpha)
print(mod[[2]]$init.alpha)


###################################################
### code chunk number 30: initial values for dispersion parameters
###################################################
print(mod[[1]]$init.eps)
print(mod[[2]]$init.eps)


###################################################
### code chunk number 31: full specification of the prior hyperparameters by the user and restart of the MCMC simulation (eval = FALSE)
###################################################
set.seed(20072011)
modContinue <- GLMM_MCMC(y = PBC910[, c("lbili", "platelet", "spiders")], 
    dist = c("gaussian", "poisson(log)", "binomial(logit)"), 
    id = PBC910[, "id"],
    x = list(lbili    = "empty", 
             platelet = "empty", 
             spiders  = PBC910[, "month"]),
    z = list(lbili    = PBC910[, "month"], 
             platelet = PBC910[, "month"], 
             spiders  = "empty"),
    random.intercept = rep(TRUE, 3),
    scale.b = list(shift = c(0.31516, 0.00765, 5.56807, -0.00665, -2.74954),
                   scale = c(0.8645, 0.0201, 0.7667, 0.0156, 3.2285)),                       
    prior.b = list(Kmax = 2, priormuQ = "independentC", 
                   delta = 1, xi = rep(0, 5), D = diag(rep(36, 5)), 
                   zeta = 6, gD = rep(0.2, 5), hD = rep(0.278, 5)), 
    prior.alpha = list(mean = 0, var = 10000),                       
    prior.eps = list(zeta = 2, g = 0.2, h = 2.76),
    init.b  = mod[[1]]$state.last.b,
    init2.b = mod[[2]]$state.last.b,                   
    init.alpha  = mod[[1]]$state.last.alpha,
    init2.alpha = mod[[2]]$state.last.alpha,                   
    init.eps = mod[[1]]$state.last.eps,
    init2.eps = mod[[2]]$state.last.eps,                   
    nMCMC = c(burn = 0, keep = 1000, thin = 10, info = 100),
    parallel = TRUE)


###################################################
### code chunk number 32: print part of w sample
###################################################
print(mod[[1]]$w_b[1:3,])


###################################################
### code chunk number 33: print part of mu sample
###################################################
print(mod[[1]]$mu_b[1:3,])


###################################################
### code chunk number 34: print part of D sample
###################################################
print(mod[[1]]$Sigma_b[1:3,])


###################################################
### code chunk number 35: print part of Q and Li samples
###################################################
print(mod[[1]]$Q_b[1:3,])
print(mod[[1]]$Li_b[1:3,])


###################################################
### code chunk number 36: print parts of alpha and sigma sample
###################################################
print(mod[[1]]$alpha[1:3,])
print(mod[[1]]$sigma_eps[1:3,])


###################################################
### code chunk number 37: print parts of hyperparameter sample
###################################################
print(mod[[1]]$gammaInv_b[1:3,])
print(mod[[1]]$gammaInv_eps[1:3,])


###################################################
### code chunk number 38: print part of mixture.b sample
###################################################
print(mod[[1]]$mixture_b[1:3,])


###################################################
### code chunk number 39: print part of Deviance
###################################################
print(mod[[1]]$Deviance[1:10], digits=6)


###################################################
### code chunk number 40: order and rank objects
###################################################
mod[[1]]$order_b[1:3,]
mod[[1]]$rank_b[1:3,]


###################################################
### code chunk number 41: running relabelling algorithm and storing of sampled component probabilities (eval = FALSE)
###################################################
mod[[1]]  <- NMixRelabel(mod[[1]], type = "stephens", keep.comp.prob = TRUE)
mod[[2]]  <- NMixRelabel(mod[[2]], type = "stephens", keep.comp.prob = TRUE)


###################################################
### code chunk number 44: 02-trace_deviance
###################################################
par(mar = c(4, 4, 0, 1) + 0.1, bty = "n")
tracePlots(mod[[1]], param = "Deviance")


###################################################
### code chunk number 45: traceplot of fixed effects and residual standard deviation (eval = FALSE)
###################################################
par(mar = c(4, 4, 4, 1) + 0.1, bty = "n")
tracePlots(mod[[1]], param = "alpha")
tracePlots(mod[[1]], param = "sigma_eps")


###################################################
### code chunk number 46: traceplot of Eb and standard deviations and correlations derived from varb (eval = FALSE)
###################################################
COL <- rep(rainbow_hcl(3, start = 30, end = 210), c(2, 2, 1))
tracePlots(mod[[1]], param = "Eb",  col = COL)
tracePlots(mod[[1]], param = "SDb", col=COL)
tracePlots(mod[[1]], param = "Corb")


###################################################
### code chunk number 47: traceplots for mixture weighs and means and standard deviations (eval = FALSE)
###################################################
tracePlots(mod[[1]], param = "w_b")
tracePlots(mod[[1]], param = "mu_b")
tracePlots(mod[[1]], param = "sd_b")


###################################################
### code chunk number 48: traceplots for mixture weighs and means and standard deviations after relabelling (eval = FALSE)
###################################################
tracePlots(mod[[1]], param = "w_b",  relabel = TRUE)
tracePlots(mod[[1]], param = "mu_b", relabel = TRUE)
tracePlots(mod[[1]], param = "sd_b", relabel = TRUE)


###################################################
### code chunk number 49: traceplots of variance hyperparameters (eval = FALSE)
###################################################
tracePlots(mod[[1]], param = "gammaInv_b")
tracePlots(mod[[1]], param = "gammaInv_eps")


###################################################
### code chunk number 50: autocorrelation of deviances (eval = FALSE)
###################################################
autocorr.plot(mod[[1]]$Deviance, lag.max = 20, col = "blue4", 
              auto.layout = FALSE, lwd = 2)


###################################################
### code chunk number 51: posterior summary statistics for basic model parameters
###################################################
print(mod)


###################################################
### code chunk number 52: coda posterior summary for regression parameters
###################################################
name.Eb <- paste("b.Mean.", 1:5, sep = "")
Regr <- cbind(mod[[1]]$mixture_b[, name.Eb], 
              mod[[1]]$alpha, mod[[1]]$sigma_eps)
colnames(Regr) <- c(paste(rep(c("lbili", "platelet", "spiders"), each = 2), 
                          ":", rep(c("Intcpt", "Slope"), 3), sep=""),
                    "lbili:eps")
Regr <- mcmc(Regr, start = mod[[1]]$nMCMC["burn"] + 1, 
                   end = mod[[1]]$nMCMC["burn"] + mod[[1]]$nMCMC["keep"])
summary(Regr)


###################################################
### code chunk number 53: coda posterior summary for SDb and Corb (eval = FALSE)
###################################################
name.SDb <- paste("b.SD.", 1:5, sep = "")
SDb <- mcmc(mod[[1]]$mixture_b[, name.SDb], 
            start = mod[[1]]$nMCMC["burn"] + 1, 
            end = mod[[1]]$nMCMC["burn"] + mod[[1]]$nMCMC["keep"])            
summary(SDb)
#
name.Corb <- paste("b.Corr.", c(2:5, 3:5, 4:5, 5), ".", rep(1:4, 4:1), sep = "")
Corb <- mcmc(mod[[1]]$mixture_b[, name.Corb],
            start = mod[[1]]$nMCMC["burn"] + 1, 
            end = mod[[1]]$nMCMC["burn"] + mod[[1]]$nMCMC["keep"])                          
summary(Corb)


###################################################
### code chunk number 54: HPD intervals for regression parameters
###################################################
HPDinterval(Regr)


###################################################
### code chunk number 55: HPD intervals for SDb and Corb (eval = FALSE)
###################################################
HPDinterval(SDb)
HPDinterval(Corb)


###################################################
### code chunk number 56: posterior densities of beta (eval = FALSE)
###################################################
COL <- rep(rainbow_hcl(3, start = 30, end = 210), each = 2)
par(mfcol = c(2, 3))
for (i in 1:6){
  densplot(Regr[, i], show.obs = FALSE, col = COL[i], lwd = 2)
  title(main = colnames(Regr)[i])
}
par(mfcol = c(1, 1))


###################################################
### code chunk number 57: posterior means of shifted and scaled mixture parameters
###################################################
NMixSummComp(mod[[1]])


###################################################
### code chunk number 58: calculation of cluster specific fitted longitudinal profiles based on posterior means without Gaussian quadrature
###################################################
delta <- 0.3
tpred <- seq(0, 30, by = delta)
fit0 <- fitted(mod[[1]], x = list("empty", "empty", tpred), 
                         z = list(tpred, tpred, "empty"))
names(fit0) <- c("lbili", "platelet", "spiders")


###################################################
### code chunk number 59: calculation of cluster specific fitted longitudinal profiles based on posterior means using Gaussian quadrature (eval = FALSE)
###################################################
tpred1 <- tpred - delta/2
tpred2 <- tpred + delta/2
fit <- fitted(mod[[1]], x =list("empty", "empty", tpred1), 
                        z =list(tpred, tpred, "empty"),
                        x2=list("empty", "empty", tpred2),
                        z2=list(tpred, tpred, "empty"),
              glmer = TRUE)
names(fit) <- c("lbili", "platelet", "spiders")


###################################################
### code chunk number 61: print parts of calculated longitudinal profiles for platelet counts
###################################################
print(fit[["platelet"]][1:10,], digits=5)


###################################################
### code chunk number 63: 03-fitted_prof
###################################################
par(mar = c(5, 4, 3, 0) + 0.1, bty = "n")
allCOL <- rainbow_hcl(1, c=30, l=85, start=200)
kCOL <- c("darkgreen", "red3")
K <- mod[[1]]$prior.b$Kmax
#
layout(autolayout(3))
plotProfiles(ip = ip, data = PBC910, var = "lbili", tvar = "month", 
   xlim = XLIM, col = allCOL, 
   xlab = "Time (months)", ylab = "Log(bilirubin)",
   auto.layout = FALSE, main = "Log(bilirubin)")
for (k in 1:K) lines(tpred, fit[["lbili"]][, k], col = kCOL[k], lwd = 2)
#
plotProfiles(ip = ip, data = PBC910, var = "platelet", tvar = "month", 
   xlim = XLIM, col = allCOL, 
   xlab = "Time (months)", ylab = "Platelet count",
   auto.layout = FALSE, main = "Platelet count")
for (k in 1:K) lines(tpred, fit[["platelet"]][, k], col = kCOL[k], lwd = 2)
#
plotProfiles(ip = ip, data = PBC910, var = "jspiders",  tvar = "month", 
   xlim = XLIM, col = allCOL, 
   xlab = "Time (months)", ylab = "Blood vessel malform. (jittered)",
   auto.layout = FALSE, main = "Blood vessel malform.")
for (k in 1:K) lines(tpred, fit[["spiders"]][, k], col = kCOL[k], lwd = 2)


###################################################
### code chunk number 64: posterior sample of the component probabilities
###################################################
print(mod[[1]]$comp.prob3[1:10, 1:8])


###################################################
### code chunk number 65: posterior component probabilities
###################################################
print(mod[[1]]$poster.comp.prob3[1:4,])


###################################################
### code chunk number 66: posterior medians of the component probabilities
###################################################
print(mod[[1]]$quant.comp.prob3[["50%"]][1:4,])


###################################################
### code chunk number 67: classification based on posterior component probabilities
###################################################
groupMean <- apply(mod[[1]]$poster.comp.prob3, 1, which.max)
pMean <- apply(mod[[1]]$poster.comp.prob3, 1, max)
table(groupMean)


###################################################
### code chunk number 68: classification based on posterior medians of component probabilities
###################################################
groupMed <- apply(mod[[1]]$quant.comp.prob3[["50%"]], 1, which.max)
pMed <- apply(mod[[1]]$quant.comp.prob3[["50%"]], 1, max)
table(groupMed)
table(groupMean, groupMed)


###################################################
### code chunk number 69: subjects with different classification using posterior mean and median
###################################################
pMeanMed <- data.frame(Mean = pMean, Median = pMed)
rownames(pMeanMed) <- unique(PBC910$id)
print(pMeanMed[groupMean != groupMed, ])


###################################################
### code chunk number 71: 04-post_dist_p
###################################################
par(mar = c(4, 4, 4, 1) + 0.1, bty = "n")
IDS <- unique(PBC910$id)
rpMean <- round(pMean, 3)
rpMed  <- round(pMed, 3)

N <- ncol(mod[[1]]$comp.prob3) / K
ID <- c(2, 7, 11)

par(mfrow = c(1, 3))
for (id in ID){
  i <- (1:N)[IDS == id]
  hist(mod[[1]]$comp.prob3[, (i - 1) * K + 1], xlim = c(0, 1), prob = TRUE, 
     xlab = expression(paste("P(u=1|", psi, ", ", theta, ", y)", sep = "")), 
     col = rainbow_hcl(1, start=60), 
     main = paste("ID ", id, "   (", rpMean[i], ",  ", rpMed[i], ")", sep=""))  
}  


###################################################
### code chunk number 72: classification using the HPD credible intervals of the component probabilities
###################################################
pHPD <- HPDinterval(mcmc(mod[[1]]$comp.prob3))
print(pHPD[1:8,])
pHPDlower <- matrix(pHPD[, "lower"], ncol=2, byrow=TRUE)
pHPDupper <- matrix(pHPD[, "upper"], ncol=2, byrow=TRUE)
rownames(pHPDlower) <- rownames(pHPDupper) <- IDS
#
groupHPD <- rep(NA, nrow(pHPDlower))
groupHPD[pHPDlower[, 1] > 0.5] <- 1
groupHPD[pHPDlower[, 2] > 0.5] <- 2
c(table(groupHPD), NAs=sum(is.na(groupHPD)))


###################################################
### code chunk number 73: variable groupMeanHPD
###################################################
groupMeanHPD <- as.character(groupMean)
groupMeanHPD[is.na(groupHPD) & groupMean == 1] <- "1_NA"
groupMeanHPD[is.na(groupHPD) & groupMean == 2] <- "2_NA"
groupMeanHPD <- factor(groupMeanHPD, levels = c("1", "2", "1_NA", "2_NA"))
table(groupMeanHPD)


###################################################
### code chunk number 74: add group indicators to data
###################################################
TAB <- table(PBC910$id)
#
PBC910$groupMean    <- factor(rep(groupMean, TAB))
PBC910$groupMed     <- factor(rep(groupMed, TAB))
PBC910$groupMeanHPD <- factor(rep(groupMeanHPD, TAB))


###################################################
### code chunk number 75: extract longitudinal profiles including the group indicators
###################################################
ip <- getProfiles(t = "month", 
   y = c("lbili", "platelet", "spiders", "jspiders",
         "groupMean", "groupMed", "groupMeanHPD"), 
   id = "id", data = PBC910)
print(ip[[1]])


###################################################
### code chunk number 77: 05-long_prof_by_group
###################################################
par(mar = c(5, 4, 3, 0) + 0.1, bty = "n")
GCOL <- rainbow_hcl(3, start=220, end=40, c=50, l=60)[c(2, 3, 1, 1)]
names(GCOL) <- levels(groupMeanHPD)
#
XLIM <- c(0, 910) / (365.25 / 12)
#
layout(autolayout(4))
#
# Log(bilirubin):
plotProfiles(ip = ip, data = PBC910, var = "lbili", tvar = "month", 
   gvar = "groupMeanHPD", xlim = XLIM, col = GCOL, 
   xlab = "Time (months)", ylab = "Log(bilirubin)",
   auto.layout = FALSE, main = "Log(bilirubin)")
#
# Legend:
plot(c(0, 100), c(0, 100), type = "n", 
     xaxt = "n", yaxt = "n", xlab = "", ylab = "")
legend(0, 90, legend = c("Group 1", "Group 2", "Uncertain classification"),
       lty = 1, lwd = 5, col = GCOL, y.intersp = 1.5, cex=1.1)
#
# Platelet count:
plotProfiles(ip = ip, data = PBC910, var = "platelet", tvar = "month", 
   gvar = "groupMeanHPD", xlim = XLIM, col = GCOL, 
   xlab = "Time (months)", ylab = "Platelet count",
   auto.layout = FALSE, main = "Platelet count")
#
# Blood vessel malformations
plotProfiles(ip = ip, data = PBC910, var = "jspiders",  tvar = "month", 
   gvar = "groupMeanHPD", xlim = XLIM, col = GCOL, 
   xlab = "Time (months)", ylab = "Blood vessel malform. (jittered)",
   auto.layout = FALSE, main = "Blood vessel malform.")


###################################################
### code chunk number 79: 06-long_prof_by_group
###################################################
par(mar = c(4, 4, 0, 0) + 0.1, bty = "n")
ips <- list()
ips[[1]] <- getProfiles(t = "month", 
   y = c("lbili", "platelet", "spiders", "jspiders", "groupMeanHPD"), 
   id = "id", data = subset(PBC910, groupMeanHPD %in% c("1", "1_NA")))
ips[[2]] <- getProfiles(t = "month", 
   y = c("lbili", "platelet", "spiders", "jspiders", "groupMeanHPD"), 
   id = "id", data = subset(PBC910, groupMeanHPD %in% c("2", "2_NA")))
#
yvars <- c("lbili", "platelet", "jspiders")
fit.yvars <- c("lbili", "platelet", "spiders")
ylabs <- c("Log(bilirubin)", "Platelet count", 
           "Blood vessel malform. (jittered)")
#
par(mfrow = c(3, 2))
for (v in 1:length(yvars)){
  for (k in 1:2){
    YLIM <- range(PBC910[, yvars[v]], na.rm = TRUE)
    plotProfiles(ip = ips[[k]], data = subset(PBC910, groupMeanHPD == k),
       var = yvars[v], tvar = "month", 
       gvar = "groupMeanHPD", xlim = XLIM, ylim = YLIM, col = GCOL, 
       xlab = ifelse(v == 3, "Time (months)", ""), 
       ylab = ifelse(k == 1, ylabs[v], ""),
       xaxt = ifelse(v == 3, "s", "n"),                 
       yaxt = ifelse(k == 1, "s", "n"),
       auto.layout = FALSE, main = "")
    lines(tpred, fit[[fit.yvars[v]]][, k], col = kCOL[k], lwd = 3)    
  }  
}  
par(mfrow = c(1, 1))


###################################################
### code chunk number 80: deviance values
###################################################
print(mod[[1]]$Deviance[1:12], digits=9)
print(mod$Deviance1[1:12], digits=9)


###################################################
### code chunk number 81: PED for model with two mixture components
###################################################
print(mod$PED)


###################################################
### code chunk number 82: running models for different numbers of mixture components (eval = FALSE)
###################################################
Devs1 <- Devs2 <- list()
for (K in 1:4){
  cat("Calculating K = ", K, "\n========================\n", sep="")
  if (K == 2){
    modK <- mod
  }else{    
    set.seed(20042008)
    modK <- GLMM_MCMC(y = PBC910[, c("lbili", "platelet", "spiders")], 
       dist = c("gaussian", "poisson(log)", "binomial(logit)"), 
       id = PBC910[, "id"],
       x = list(lbili    = "empty", 
               platelet = "empty", 
               spiders  = PBC910[, "month"]),
       z = list(lbili    = PBC910[, "month"], 
                platelet = PBC910[, "month"], 
                spiders  = "empty"),
       random.intercept = rep(TRUE, 3),
       prior.b = list(Kmax = K), 
       nMCMC = c(burn = 100, keep = 1000, thin = 10, info = 100),
       parallel = TRUE)
  }  
  Devs1[[K]] <- modK[["Deviance1"]]
  Devs2[[K]] <- modK[["Deviance2"]]
  
  if (K == 1){
    PED <- as.data.frame(matrix(modK[["PED"]], nrow=1))
    colnames(PED) <- names(modK[["PED"]])
  }else PED <- rbind(PED, modK[["PED"]])  
  
  rm(list = "modK")  
}


###################################################
### code chunk number 83: penalized expected deviances for models with different values of K
###################################################
print(PED, digits=6)


###################################################
### code chunk number 84: options
###################################################
options(digits = 5)


###################################################
### code chunk number 85: posterior summary statistics for the difference between the observed data deviances in models with K 2 and K 1
###################################################
summaryDiff(Devs1[[2]], Devs1[[1]])


###################################################
### code chunk number 86: posterior summary statistics for the difference between the observed data deviances in models with K 3 and K 2
###################################################
summaryDiff(Devs1[[3]], Devs1[[2]])


###################################################
### code chunk number 87: posterior summary statistics for the difference between the observed data deviances in models with K 4 and K 3
###################################################
summaryDiff(Devs1[[4]], Devs1[[3]])


###################################################
### code chunk number 89: 07-cdf_deviance
###################################################
par(mar = c(4, 4, 1, 1) + 0.1, bty = "n")
COL <- terrain_hcl(4, c = c(65, 15), l = c(45, 80), power = c(0.5, 1.5))
plot(c(14000, 14275), c(0, 1), type="n", 
     xlab="Deviance", ylab="Posterior CDF")
for (K in 1:4){
  medDEV <- median(Devs1[[K]])
  ECDF <- ecdf(Devs1[[K]])
  plot(ECDF, col = COL[K], lwd = 2, add = TRUE)
  text(medDEV + 0.5, 0.5, labels = K)
}  
