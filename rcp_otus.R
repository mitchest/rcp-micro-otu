#install.packages("RCPmod")
library(RCPmod)

source("rcp_otus_funs.R")

otu_matrix <- read.csv("otu_matrix.csv", header = T)#, stringsAsFactors = F)

names(otu_matrix)[3780:3786]
otu_covar <-otu_matrix[,3783:3786]
otu_abund <-otu_matrix[,-3783:-3786]

# just exploring the data a little
boxplot(P ~ sampleCode, data=otu_covar)
boxplot(Ca ~ sampleCode, data=otu_covar)
boxplot(P ~ size, data=otu_covar)
boxplot(Ca ~ size, data=otu_covar)

occurence <- sort(colSums(ifelse(otu_abund>0,1,0)))

sum(occurence > 3)


# data prep for models ----------------------------------------------------

# i like covars to come first...
otu_matrix <- data.frame(otu_covar, otu_abund)
# get rid of the data without covar data...
otu_matrix <- otu_matrix[apply(X = otu_matrix, MARGIN = 1, FUN = function(X){!any(is.na(X))}),]
# define where abundance data starts in the data
n.abund <- 5

# choose otus to model with
# create a list of otus with occurance > n to use for modelling
otus.n <- 10 # change this to desired minimum occurance
otu.count <- data.frame(count=sort(colSums(ifelse(otu_matrix>0,1,0)[,n.abund:ncol(otu_matrix)]), decreasing=T))
otu.count$otu <- row.names(otu.count)
model.otus.vector <- otu.count$otu[otu.count$count>otus.n]
model.otus.string <- paste0("cbind(", paste(model.otus.vector, collapse=","),")")

# choose which covariates to use
model.covariates.string <- "P + P.2 + Ca + Ca.2"
model.covariates.vector <- c("P","Ca")

covar.data <- otu_matrix[,model.covariates.vector]
# calculate quadratic polynomial cols
covar.data <- data.frame(poly(covar.data$P, 2),
                        poly(covar.data$Ca, 2))
names(covar.data) <- c("P","P.2",
                      "Ca","Ca.2")


## convert categorical variables to factors if you have them

# generate model data
model.data = data.frame(otu_matrix[,model.otus.vector], covar.data)

# define model form
RCP.form = paste0(model.otus.string,"~","1","+",model.covariates.string)
my.cont = list(maxit=3000, penalty=0.0001, penalty.tau=10, penalty.gamma=10)

# fit models (with n clusters hard set)
# I would paralellise this when it gets serious
fm_3rcp <- regimix(form.RCP = RCP.form, data = model.data, nRCP = 3, 
                   dist="Poisson", control=my.cont, inits="noPreClust", titbits=TRUE)
fm_10rcp <- regimix(form.RCP = RCP.form, data = model.data, nRCP = 10, 
                   dist="Poisson", control=my.cont, inits="noPreClust", titbits=TRUE)

regifits <- list()
for (i in 2:20) {
  fm <- regimix.multifit(form.RCP = RCP.form, data = model.data, nRCP = i, 
          dist="Poisson", control=my.cont, inits="noPreClust", titbits=TRUE, mc.cores = 4)
  regifits[[paste0("fm_",i,"RCP")]] <- fm
}
save(regifits, file="regifits_2-20RCP_mult.RData")

# plot the initial results
load("regifits_2-20RCP_mult.RData")
plot(1~1, type = "n", xlim = c(1,41), ylim = c(2000000,5000000))
for (i in 1:19) {
  points(unlist(lapply(regifits[[i]], function(x)x$BIC)) ~ rep(i+1, 10), pch=16)
}



# data prep for partial effects -------------------------------------------

# make a new data matrix
new.dat <- matrix(NA, ncol=ncol(covar.data), nrow=1000)

# get the col means
X.means <- colMeans(covar.data)
# get the col medians
X.meds = colMedian(covar.data)
# get col mins
X.mins <- colMin(covar.data)

# add meds or means to the new data matrix
for (i in 1:length(X.meds)) {
  new.dat[,i] <- X.meds[i]
}

# give it the names used for fitting
new.dat <- data.frame(new.dat)
names(new.dat) <- names(covar.data)

## once you have categorical variables
# # hack the catgorical X's into shape - need a better solution for this
# # choose the group to be the prediciton level
# sort(table(covar.data$XXXXXX))
# new.dat$XXXXXX <- as.factor(covar.data$XXXXXX)[1]
# # ## or test the XXXXXX varaible
# # new.dat = new.dat[1:length(unique(covar.data$XXXXXX)),]
# # new.dat$XXXXXX = as.factor(unique(covar.data$XXXXXX))
# # new.dat.XXXXXX = new.dat

# make the polynomials - ensure names match original call
new.dat <- data.frame(predict(poly(covar.data$P, 2), newdata=new.dat$P),
                     predict(poly(covar.data$Ca, 2), newdata=new.dat$Ca))
names(new.dat) <- c("P","P.2",
                      "Ca","Ca.2")



# plot partial effects ----------------------------------------------------

par(mfcol=c(3,4))
for (i in names(covar.data)) {
  ParPlotIndiv.regimix(fm_3rcp, covar.data, new.dat, i, 1:3)
}

par(mfrow=c(1,1))
for (i in names(covar.data)) {
  ParPlotMany.regimix(fm_3rcp, covar.data, new.dat, i, 1:3)
}

par(mfrow=c(2,2))
for (i in names(covar.data)) {
  ParPlotMany.regimix(fm_3rcp, covar.data, new.dat, i, 1:3)
}

par(mfrow=c(1,2))
for (i in names(covar.data)) {
  ParPlotMany.regimix(fm_10rcp, covar.data, new.dat, i, 1:10)
}
