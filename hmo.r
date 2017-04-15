## Data analytics question and Analysis
######################################################
## A large HMO wants to know what patient and physician factors are most related to 
## whether a patient’s lung cancer goes into remission after 
## treatment as part of a larger study of treatment outcomes and quality of life in patients with lunger cancer.
## 
## Using CancerSurvival xl (coverting it to csv) data
## Name ## Mayukh Mukhopadhyay
###################################################################

install.packages("ggplot2")
install.packages("GGally")
install.packages("reshape2")
install.packages("lme4")
install.packages("compiler")
install.packages("parallel")
install.packages("boot")


require(ggplot2)
require(GGally)
require(reshape2)
require(lme4)
require(compiler)
require(parallel)
require(boot)

cancerdata <- read.csv("~/CancerSurvival.csv")
cancerdata <- within(cancerdata, {
  Married <- factor(Married, levels = 0:1, labels = c("no", "yes"))
  DID <- factor(DID)
  HID <- factor(HID)
})
ggpairs(cancerdata[, c("IL6", "CRP", "LengthofStay", "Experience")])
ggplot(cancerdata, aes(x = CancerStage, y = LengthofStay)) +
  stat_sum(aes(size = ..n.., group = 1)) +
  scale_size_area(max_size=10)

 # il6 crp cancerstage plots
tmp <- melt(cancerdata[, c("CancerStage", "IL6", "CRP")], id.vars="CancerStage")
ggplot(tmp, aes(x = CancerStage, y = value)) +
  geom_jitter(alpha = .1) +
  geom_violin(alpha = .75) +
  facet_grid(variable ~ .) +
  scale_y_sqrt()
  
  #boxplot
tmp <- melt(cancerdata[, c("remission", "IL6", "CRP", "LengthofStay", "Experience")],
  id.vars="remission")
ggplot(tmp, aes(factor(remission), y = value, fill=factor(remission))) +
  geom_boxplot() +
  facet_wrap(~variable, scales="free_y")
  
# estimate the model and store results in m for nAGQ = 10
m <- glmer(remission ~ IL6 + CRP + CancerStage + LengthofStay + Experience +
    (1 | DID), data = cancerdata, family = binomial, control = glmerControl(optimizer = "bobyqa"),
    nAGQ = 10)

# BOBYQA (Bound Optimization BY Quadratic Approximation)

# print the mod results without correlations among fixed effects
print(m, corr = FALSE)

# Generalized linear mixed model fit by maximum likelihood (Adaptive Gauss-Hermite Quadrature, nAGQ
#   = 10) [glmerMod]
#  Family: binomial  ( logit )
# Formula: remission ~ IL6 + CRP + CancerStage + LengthofStay + Experience +      (1 | DID)
#    Data: cancerdata
#       AIC       BIC    logLik  deviance  df.resid 
#  7397.276  7460.733 -3689.638  7379.276      8516 
# Random effects:
#  Groups Name        Std.Dev.
#  DID    (Intercept) 2.015   
# Number of obs: 8525, groups:  DID, 407
# Fixed Effects:
#    (Intercept)             IL6             CRP   CancerStageII  CancerStageIII   CancerStageIV  
#       -2.05271        -0.05677        -0.02148        -0.41393        -1.00346        -2.33704  
#   LengthofStay      Experience  
#       -0.12118         0.12009 

# estimate the model and store results in m for nAGQ = 50
m <- glmer(remission ~ IL6 + CRP + CancerStage + LengthofStay + Experience +
    (1 | DID), data = cancerdata, family = binomial, control = glmerControl(optimizer = "bobyqa"),
    nAGQ = 50)

# BOBYQA (Bound Optimization BY Quadratic Approximation)

# print the mod results without correlations among fixed effects
print(m, corr = FALSE)

#Generalized linear mixed model fit by maximum likelihood (Adaptive Gauss-Hermite Quadrature, nAGQ
#  = 50) [glmerMod]
# Family: binomial  ( logit )
#Formula: remission ~ IL6 + CRP + CancerStage + LengthofStay + Experience +      (1 | DID)
#   Data: cancerdata
#      AIC       BIC    logLik  deviance  df.resid 
# 7397.233  7460.690 -3689.617  7379.233      8516 
#Random effects:
# Groups Name        Std.Dev.
# DID    (Intercept) 2.015   
#Number of obs: 8525, groups:  DID, 407
#Fixed Effects:
#   (Intercept)             IL6             CRP   CancerStageII  CancerStageIII   CancerStageIV  
#      -2.05325        -0.05677        -0.02148        -0.41393        -1.00349        -2.33711  
#  LengthofStay      Experience  
#      -0.12118         0.12011


# getting confidence intervals (CIs)
se <- sqrt(diag(vcov(m)))
# table of estimates with 95% CI
(tab <- cbind(Est = fixef(m), LL = fixef(m) - 1.96 * se, UL = fixef(m) + 1.96 *
    se))
	
##                        Est          LL           UL
## (Intercept)    -2.05325039 -3.09521978 -1.011280996
## IL6            -0.05677292 -0.07934869 -0.034197141
## CRP            -0.02148286 -0.04151069 -0.001455033
## CancerStageII  -0.41393231 -0.56242696 -0.265437658
## CancerStageIII -1.00348621 -1.19612106 -0.810851355
## CancerStageIV  -2.33710519 -2.64690344 -2.027306943
## LengthofStay   -0.12118413 -0.18710508 -0.055263171
## Experience      0.12010779  0.06628398  0.173931601

# odds ratios instead of coefficients on the logit scale	
exp(tab)

##                       Est         LL        UL
## (Intercept)    0.12831714 0.04526506 0.3637527
## IL6            0.94480860 0.92371778 0.9663810
## CRP            0.97874625 0.95933908 0.9985460
## CancerStageII  0.66104570 0.56982444 0.7668702
## CancerStageIII 0.36659917 0.30236479 0.4444795
## CancerStageIV  0.09660689 0.07087033 0.1316897
## LengthofStay   0.88587083 0.82935658 0.9462361
## Experience     1.12761839 1.06853011 1.1899742




##############multilevel bootstrapping with 100 replicates

sampler <- function(dat, clustervar, replace = TRUE, reps = 1) {
    cid <- unique(dat[, clustervar[1]])
    ncid <- length(cid)
    recid <- sample(cid, size = ncid * reps, replace = TRUE)
    if (replace) {
        rid <- lapply(seq_along(recid), function(i) {
            cbind(NewID = i, RowID = sample(which(dat[, clustervar] == recid[i]),
                size = length(which(dat[, clustervar] == recid[i])), replace = TRUE))
        })
    } else {
        rid <- lapply(seq_along(recid), function(i) {
            cbind(NewID = i, RowID = which(dat[, clustervar] == recid[i]))
        })
    }
    dat <- as.data.frame(do.call(rbind, rid))
    dat$Replicate <- factor(cut(dat$NewID, breaks = c(1, ncid * 1:reps), include.lowest = TRUE,
        labels = FALSE))
    dat$NewID <- factor(dat$NewID)
    return(dat)
}

set.seed(20)
tmp <- sampler(cancerdata, "DID", reps = 100)
bigdata <- cbind(tmp, cancerdata[tmp$RowID, ])

## First we store the estimates from our original model, which we will use as start values for the bootstrap models. 
## Then we make a local cluster with 2 nodes (the number of processors on our machine). 
## Next, we export the data and load the lme4 package on the cluster. 
## Finally, we write a function to fit the model and return the estimates. 
## The call to glmer() is wrapped in try because not all models may converge on the resampled data. 
## This catches the error and returns it, rather than stopping processing.

f <- fixef(m)
r <- getME(m, "theta")

cl <- makeCluster(2)
clusterExport(cl, c("bigdata", "f", "r"))
clusterEvalQ(cl, require(lme4))
## [[1]]
## [1] TRUE
## 
## [[2]]
## [1] TRUE


myboot <- function(i) {
    object <- try(glmer(remission ~ IL6 + CRP + CancerStage + LengthofStay +
        Experience + (1 | NewID), data = bigdata, subset = Replicate == i, family = binomial,
        nAGQ = 1, start = list(fixef = f, theta = r)), silent = TRUE)
    if (class(object) == "try-error")
        return(object)
    c(fixef(object), getME(object, "theta"))
}


#### Execute with your own risk...My laptop with 2 cores and 4GB got blackedout  5 times and returned output after 2 hours
start <- proc.time()
res <- parLapplyLB(cl, X = levels(bigdata$Replicate), fun = myboot)
end <- proc.time()

# shut down the cluster
stopCluster(cl)

# calculate proportion of models that successfully converged
success <- sapply(res, is.numeric)
mean(success)

# combine successful results
bigres <- do.call(cbind, res[success])

# calculate 2.5th and 97.5th percentiles for 95% CI
(ci <- t(apply(bigres, 1, quantile, probs = c(0.025, 0.975))))

##                       2.5%     97.5%
## (Intercept)       -3.61982 -0.985404
## IL6               -0.08812 -0.029664
## CRP               -0.04897  0.006824
## CancerStageII     -0.60754 -0.228019
## CancerStageIII    -1.30217 -0.754609
## CancerStageIV     -2.91414 -2.002643
## LengthofStay      -0.21596 -0.046420
## Experience         0.06819  0.207223
## NewID.(Intercept)  2.03868  2.476366



# All results
finaltable <- cbind(Est = c(f, r), SE = c(se, NA), BootMean = rowMeans(bigres),
    ci)
# round and print
round(finaltable, 3)

##                    Est    SE BootMean   2.5%  97.5%
## (Intercept)     -2.053 0.531   -2.205 -3.620 -0.985
## IL6             -0.057 0.012   -0.059 -0.088 -0.030
## CRP             -0.021 0.010   -0.022 -0.049  0.007
## CancerStageII   -0.414 0.076   -0.417 -0.608 -0.228
## CancerStageIII  -1.003 0.098   -1.043 -1.302 -0.755
## CancerStageIV   -2.337 0.158   -2.460 -2.914 -2.003
## LengthofStay    -0.121 0.034   -0.142 -0.216 -0.046
## Experience       0.120 0.027    0.128  0.068  0.207
## DID.(Intercept)  2.015    NA    2.263  2.039  2.476