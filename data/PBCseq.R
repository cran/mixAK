PBCseq <- read.table("PBCseq.dat", skip=86, header=FALSE, na.string=".")
colnames(PBCseq) <- c("id", "fu.days", "status", "drug", "age", "sex", "day", "ascites", "hepatom", "spiders", "edema", "bili", "chol", "albumin", "alk.phos", "sgot", "platelet", "protime", "stage")


PBCseq <- transform(PBCseq, age   = age / 365.25,             ## age in days -> age in years
                            month = day / (365.25/12))      
PBCseq <- transform(PBCseq, fstatus  = factor(status, levels=0:2, labels=c("censored", "ltx.censored", "death")),
                            fdrug    = factor(drug, levels=0:1, labels=c("placebo", "D-penicillamine")),
                            fsex     = factor(sex, levels=0:1, labels=c("male", "female")),
                            fascites = factor(ascites, levels=0:1, labels=c("no", "yes")),
                            fhepatom = factor(hepatom, levels=0:1, labels=c("no", "yes")),
                            fspiders = factor(spiders, levels=0:1, labels=c("no", "yes")),
                            fedema   = factor(edema, levels=c(0, 0.5, 1)),                    
                            fstage   = factor(stage))


##### Log transformations
PBCseq <- transform(PBCseq, lbili     = log(bili),
                            lalbumin  = log(albumin),
                            lalk.phos = log(alk.phos),
                            lchol     = log(chol),
                            lsgot     = log(sgot),
                            lplatelet = log(platelet),
                            lprotime  = log(protime))


##### Censoring indicators
PBCseq[, "delta.death"]     <- 1 * (PBCseq[, "fstatus"] == "death")
PBCseq[, "delta.ltx.death"] <- 1 * (PBCseq[, "fstatus"] %in% c("death", "ltx.censored"))


##### Variable which says when was the last known moment that the patient lived without LTx
PBCseq[, "alive"] <- PBCseq[, "fu.days"]
PBCseq[PBCseq[, "delta.ltx.death"] == 1, "alive"] <- PBCseq[PBCseq[, "delta.ltx.death"] == 1, "alive"] - 1


#summary(PBCseq[, "alive"])
#summary(PBCseq)

##### Re-shuffle variable names
PBCseq <- PBCseq[, c("id", "sex", "fsex", "drug", "fdrug", "age", "fu.days", "alive", "status", "fstatus", "delta.death", "delta.ltx.death",
                     "day", "month",
                     "ascites", "fascites", "hepatom", "fhepatom", "spiders", "fspiders", "edema", "fedema", "stage", "fstage",
                     "bili", "lbili", "albumin", "lalbumin", "alk.phos", "lalk.phos", "chol", "lchol", "sgot", "lsgot",
                     "platelet", "lplatelet", "protime", "lprotime")]
