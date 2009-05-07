TandmobEmer <- read.table("Tandmob.dat", skip = 72, header = TRUE, as.is = c(TRUE, FALSE, TRUE, TRUE, rep(FALSE, 3), rep(TRUE, 149)))

### Take only needed variables
TandmobEmer <- subset(TandmobEmer, select=c("IDNR", "GENDER", "GENDERNum", "DOB", "PROVINCE", "EDUC", "STARTBR",
                                            paste("EBEG.", rep(c(10, 20, 30, 40), 7) + rep(1:7, each=4), sep=""),
                                            paste("EEND.", rep(c(10, 20, 30, 40), 7) + rep(1:7, each=4), sep="")))

### Change all left-censored emergence times into interval-censored with the lower limit equal to 5 years of age
TandmobEmer[, paste("EBEG.", rep(c(10, 20, 30, 40), 7) + rep(1:7, each=4), sep="")][is.na(TandmobEmer[, paste("EBEG.", rep(c(10, 20, 30, 40), 7) + rep(1:7, each=4), sep="")])] <- 5

### Create censoring indicator
TandmobEmer <- cbind(TandmobEmer, matrix(3, nrow=nrow(TandmobEmer), ncol=4*7))
colnames(TandmobEmer) <- c("IDNR", "GENDER", "GENDERNum", "DOB", "PROVINCE", "EDUC", "STARTBR",
                           paste("EBEG.", rep(c(10, 20, 30, 40), 7) + rep(1:7, each=4), sep=""),
                           paste("EEND.", rep(c(10, 20, 30, 40), 7) + rep(1:7, each=4), sep=""),
                           paste("CENSOR.", rep(c(10, 20, 30, 40), 7) + rep(1:7, each=4), sep=""))
TandmobEmer[, paste("CENSOR.", rep(c(10, 20, 30, 40), 7) + rep(1:7, each=4), sep="")][is.na(TandmobEmer[, paste("EEND.", rep(c(10, 20, 30, 40), 7) + rep(1:7, each=4), sep="")])] <- 0

