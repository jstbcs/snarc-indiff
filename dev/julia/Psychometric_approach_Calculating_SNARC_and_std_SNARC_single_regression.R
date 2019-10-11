# K Cipora et al., 2018
# Retrieved from https://osf.io/n7szg/

### Calculating unstandardized and standardized SNARC slopes. The scipt calculates confidence intervals for the individual SNARC slopes at selected level using the psychometric approach.

# Data should be saved as a tab delimited text. Decimals should be separated with ".". Note that each row of the data file must denote 1 experimental trial.Within each participant trials should be ordered as they were presented in the experiment. Relevant variables for the calculations below are as follows: 
# part.code - participant's ID
# number - number presented in a given trial
# resp.side - response side (possible values: right vs left)
# rt - reaction time
# filter.seq - result of sequential flitering of outlier reaction times (refer to "SNARC_calculations_sequential_filter" script)
# Note also that (as in the sample file) there may be some more variables in the file and it does not disturb anything.
# Note that R is case sensitive!
################################################

# Loading required packages
require(reshape2)
require(plyr)
require(QuantPsyc)
require(GeneNet)
require(Hmisc)

# Loading the data
# setwd("C:/My_SNARC_data")
raw_data <- read.table("Cipora_2014_raw_data.txt", sep="\t", dec=".", header = TRUE)

# Changing part.code variable to be a factor. In case one uses pure numerical codes
raw_data$part.code <- as.factor(raw_data$part.code)

# Creating new data frame comprising only filtered trials (i.e. in which the variable filter.seq == 0)
filtered_raw_data <- raw_data[raw_data$filter.seq == 1, ]


## Julia data inspection
head(raw_data)
# Only use accurate trials, like they did
sort(with(raw_data, tapply(correct.experimental, part.code, mean)))
sort(with(filtered_raw_data, tapply(correct.experimental, part.code, mean)))

#Use RTs between 300 and 1500 ms, maybe? RTs seem to be quite fast
hist(raw_data$rt)
hist(filtered_raw_data$rt)

nrow(raw_data)
nrow(filtered_raw_data)

## I want data with accuracy and >200ms filter but not the 3SD filter they used/.

dat <- subset(raw_data, correct.experimental == 1 & filter.correct.exp.antic == 1)

nrow(dat)

with(dat, tapply(rt, list(number, resp.side, part.code), length))

## Let's try some modeling

## Define Design Matrices
prep.models <- function(id, side, number){
  
  #define the direction of the effect
  l <- levels(side)[1]
  h <- levels(side)[2]
  
  id <- as.numeric(factor(id))
  
  #design matrices
  I <- length(unique(id))
  R <- length(id)
  X.full <- matrix(nrow = R, ncol = 4 * I + 3, 0)
  for (r in 1:R){
    # baseline RT
    X.full[r, id[r]] <- 1
    # effect of side
    if (side[r] == l) {
        X.full[r, I + 1] <- 1
        X.full[r, I + 1 + id[r]] <- 1}
    # else{
    #     X.full[r, I + 1] <- -1/2
    #     X.full[r, I + 1 + id[r]] <- -1/2
    #     }
    # Effect of number
    X.full[r, 2*I + 2] <- number[r]
    X.full[r, 2*I + 2 + id[r]] <- number[r]
    ## interaction between number and side, i.e. SNARC effect
    if (side[r] == l) {
      X.full[r, 3*I + 3] <- number[r]
      X.full[r, 3*I + 3 + id[r]] <- number[r]}
    # else{
    #     X.full[r, 3*I + 3] <- -number[r]/2
    #     X.full[r, 3*I + 3 + id[r]] <- -number[r]/2
    #   }
    }
  
  gMap.full <- c(rep(0, I), 1, rep(2, I), 3, rep(4, I), 5, rep(6, I))
  
  # X.one <- matrix(nrow = R, ncol = I + 1, 0)
  # for (r in 1:R){
  #   X.one[r, id[r]] <- 1
  #   if (fcond[r] == h) {
  #     X.one[r, I + 1] <- 1
  #   }}
  # 
  # gMap.one <- c(rep(0, I), 1)
  
  X.null <- matrix(nrow = R, ncol = 3*I + 2, 0)
  for (r in 1:R){
    # baseline RT
    X.null[r, id[r]] <- 1
    # effect of side
    if (side[r] == l) {
      X.null[r, I + 1] <- 1
      X.null[r, I + 1 + id[r]] <- 1}
    # else{
    #     X.null[r, I + 1] <- -1/2
    #     X.null[r, I + 1 + id[r]] <- -1/2
      # }
    # Effect of number
    X.null[r, 2*I + 2] <- number[r]
    X.null[r, 2*I + 2 + id[r]] <- number[r]
  }
  
  gMap.null <- c(rep(0, I), 1, rep(2, I), 3, rep(4, I))
  
  return(list(X.full = X.full
              , gMap.full = gMap.full
              # , X.one = X.one
              # , gMap.one = gMap.one
              , X.null = X.null
              , gMap.null = gMap.null
              , R = R
              , I= I
              , sub = id))
}

prep <- prep.models(dat$part.code, side = dat$resp.side, number = dat$number)
meanScale <- 1/10
effectScale <- 1/12

make.bf <- function(y, meanScale, effectScale, prep = prep.1, iter = 10000, burnin = 1000)
{
  keep <- (burnin + 1) : iter
  #posterior computation from the unconstrained model
  mcmc.full <- BayesFactor::nWayAOV(y
                                    , prep$X.full
                                    , prep$gMap.full
                                    , rscale = c(1, .5, .5, .5, .5, meanScale, effectScale)
                                    , posterior = T
                                    , iterations = iter)
  #Bayes factor estimation for the unconstrained model
  bf.full <- BayesFactor::nWayAOV(y
                                  , prep$X.full
                                  , prep$gMap.full
                                  , rscale = c(1, .5, .5, .5, .5, meanScale, effectScale)
                                  , posterior = F
                                  , iterations = iter)
  #Bayes factor computation for the common-effect model (but with unconstrained common effect, i.e. can be positive or negative)
  # bf.one <- BayesFactor::nWayAOV(y
  #                                , prep$X.one
  #                                , prep$gMap.one
  #                                , rscale = c(1, meanScale)
  #                                , posterior = F
  #                                , iterations = iter)
  # #posterior estimates from the common-effect model (but with unconstrained common effect, i.e. can be positive or negative)
  # mcmc.one <- BayesFactor::nWayAOV(y
  #                                  , prep$X.one
  #                                  , prep$gMap.one
  #                                  , rscale = c(1, meanScale)
  #                                  , posterior = T
  #                                  , iterations = iter)
  #Bayes factor computation for the null model (with individual baseline for each participant, but no effect)
  bf.null <- BayesFactor::nWayAOV(y
                                  , prep$X.null
                                  , prep$gMap.null
                                  , rscale = c(1, rep(.5, 4))
                                  , posterior = F
                                  , iterations = iter)
  
  #Extracting individual effect estimates from the output above
  i.theta0 <- 3 * prep$I + 4 #index for the common effect
  i.theta <- (3 * prep$I + 5):(4 * prep$I + 4) #indices for the individual effects
  
  myTheta <- mcmc.full[keep, i.theta] + mcmc.full[keep, i.theta0] #add oberall effect to individuals' deviations
  
  #Evaluate how often all individuals' theta_i are positive in the posterior
  good <- myTheta < 0 #evaluate their sign
  all.good <- apply(good, 1, mean) #evaluate how often all theta estimates are postitive
  PostCount <- mean(all.good == 1) #Posterior probability of all theta_i being positive
  
  #Evaluate how often all individuals' theta_i are positive in the prior
  #prior settings
  R <- iter * 10
  beta <- .5 * effectScale^2
  alpha <- .5
  PriorCount <- prior.p.greater(M = R, I = prep$I, a = alpha, b = beta, rscale = meanScale)
  
  #Calculate Bayes factors
  bf.FP <- PriorCount/PostCount
  bf.F0 <- exp(bf.full$bf - bf.null$bf)
  # bf.F1 <- exp(bf.full$bf - bf.one$bf)
  
  #Return posterior estimates from the unconstrained model
  m <- colMeans(myTheta)
  new.sd <- sd(m)
  new.mean <- mean(mcmc.full[keep, i.theta0])
  # i.g3 <- 2 * prep$I + 6
  # new.pop.sd <- mean(mcmc.full[keep, i.g3])
  return(list(ind.effects = m
              , posterior.mean = new.mean, posterior.sd = new.sd
              , bfs = c(bf.pu =  1/ bf.FP, bf.0u = 1 / bf.F0)
              , theta = myTheta
              , bf.unconstrained = bf.full, bf.null = bf.null
              , mcmc.unconstrained = mcmc.full[keep, ]
              , prior.prob = PriorCount, posterior.prob = PostCount
              , design.matrices = prep))
}

quid <- function(id #vector with participant ID, can be a factor or numeric
                 , side #vector with condition ID, can be a factor or numeric
                 , number
                 , rt #vector with response times per trial, has to be numeric
                 , prior = c(1/10, 1/12) #prior scale settings for the standard deviation of the overall effect and the individual effects
                 , iter = 10000 #number of iterations for the posterior sampling
                 , burnin = 1000 #number of to be discarded burn-in iterations
                 , messages = TRUE
){
  
  require(BayesFactor, quietly = messages)
  
  # if(!(is.numeric(condition) | is.factor(condition))) stop("Condition has to be numeric or a factor.")
  # if(!(is.numeric(id) | is.factor(id))) stop("Id has to be numeric or a factor.")
  if(!(is.numeric(rt))) stop("Rt has to be numeric.")
  
  # fcond <- factor(condition)
  
  # if(length(fcond) != length(id)) stop("Your condition vector does not match your id vector.")
  # if(length(levels(fcond)) != 2) stop("Your condition vector has more than two levels.")
  
  sub <- as.numeric(factor(id))
  
  # get design matrix
  prep.1 <- prep.models(id = id, side = side, number = number)
  
  # get bayes factors
  make.bf(y = rt, meanScale = prior[1], effectScale = prior[2], prep = prep.1, iter = iter, burnin = burnin)
}

out <- quid(id = dat$part.code, side = dat$resp.side, number = dat$number, rt = dat$rt)

## check
i.nu0 <- 1 * prep$I + 2 #index for the common effect
i.nu <- (1 * prep$I + 3):(2 * prep$I + 2) #indices for the individual effects

myNu <- mcmc.full[keep, i.nu] + mcmc.full[keep, i.nu0] #add oberall effect to individuals' deviations
plot(colMeans(myNu), testSNARC$side.slope)

i.delta0 <- 2 * prep$I + 3 #index for the common effect
i.delta <- (2 * prep$I + 4):(3 * prep$I + 3) #indices for the individual effects

myDelta <- mcmc.full[keep, i.delta] + mcmc.full[keep, i.delta0] #add oberall effect to individuals' deviations
colMeans(myDelta)[1]
plot(colMeans(myDelta), testSNARC$number.slope)

i.mu <- 1 #index for the common effect
i.alpha <- (2):(prep$I + 1) #indices for the individual effects

myAlpha <- mcmc.full[keep, i.mu] + mcmc.full[keep, i.alpha] #add oberall effect to individuals' deviations
plot(colMeans(myAlpha), testSNARC$p.intercept)

plot(out$ind.effects, testSNARC$unstdSNARC.slope)

## Very basic data inspection - building the table with number of trials per participant in total (i.e. including those where filter equal to 0) and number of valid trials per participant. It may be useful for checking whether there are some largely incomplete data

basic_inspection <- merge(plyr::rename(aggregate(rt ~ part.code, raw_data, length), c("rt"="total.no.of.trials")), plyr::rename(aggregate(filter.seq ~ part.code, raw_data, sum), c("filter.seq"="no.of.valid.trials")), by = "part.code")
basic_inspection

###############################################
## reorder participants
dat$sub <- prep$sub
dat <- dat[order(dat$sub), ]

# Averaging reaction times for each [participant x number x response side] configuration and forming the data frame in convenient manner, calculating dRTs
aggregated_data <- aggregate(rt~ resp.side + number + part.code, dat, mean)
aggregated_data_dRT <- dcast(aggregated_data, part.code + number ~ resp.side)
aggregated_data_dRT$dRT <- aggregated_data_dRT$right - aggregated_data_dRT$left

# Calculating unstandardized slopes with main effect of side and number
testSNARC <- ddply(aggregated_data, "part.code", function(x){
  regression <- lm(x$rt ~ x$number * x$resp.side, data = aggregated_data)
  ind.slope <- regression$coefficients[[4]]
  side.slope <- regression$coefficients[[2]]
  number.slope <- regression$coefficients[[3]]
  p.intercept <- regression$coefficients[[1]] 
  data.frame(unstdSNARC.slope=ind.slope, side.slope, number.slope, p.intercept)
})
testSNARC


# Calculating unstandardized slopes
unstandardizedSNARC <- ddply(aggregated_data_dRT, "part.code", function(x){
  regression <- lm(x$dRT ~ x$number, data = aggregated_data_dRT)
  ind.slope <- regression$coefficients[[2]]
  data.frame(unstdSNARC.slope=ind.slope)
})
unstandardizedSNARC

plot(sort(unstandardizedSNARC$unstdSNARC.slope), pch = 19)
abline(h = 0)

# testing for the snarc signficance at the sample level
t.test(unstandardizedSNARC$unstdSNARC.slope, mu=0)


# # Calculating standardized slopes (with Fisher-z transfofrmation)
# standardizedSNARC <- ddply(aggregated_data_dRT, "part.code", function(x){
#   regression <- lm(x$dRT ~ x$number, data = aggregated_data_dRT)
#   std.ind.slope <- lm.beta(regression)
#   data.frame(std.slopes=std.ind.slope)
# })
# 
# standardizedSNARC$standardized.SNARC <- z.transform(standardizedSNARC$std.slopes)
# 
# standardizedSNARC
# plot(sort(standardizedSNARC$standardized.SNARC), pch = 19)
# abline(h = 0)
# 
# # testing for the standardized SNARC signficance at the sample level
# t.test(standardizedSNARC$standardized.SNARC, mu=0)


# # Saving the data files with unstandardized and standardized SNARC slopes
# 
# write.table(unstandardizedSNARC, file="Unstandardized_SNARC.txt", sep = "\t", dec=".", quote = FALSE, row.names=FALSE)
# 
# write.table(standardizedSNARC, file="Standardized_SNARC.txt", sep = "\t", dec=".", quote = FALSE, row.names=FALSE)

#################################################
### Calculating confidence intervals
## Unstandardized SNARC
# Type in unstandardized SNARC reliability
# Note that reliability can be calculated using other scripts provided
# Currently the value is the one obtained from the sample data file

unstd_rel <- .78

# If you want to provide the SD slope value based some external data (e.g. from other experiment in which more people were tested), please type it below, otherwise SD slope from the file being analysed.

ext_SD_unstd_slope <- "place the numerical value without a quotation mark instead of this text"
if(is.numeric(ext_SD_unstd_slope)) SD_unstd_slope <- ext_SD_unstd_slope  else SD_unstd_slope <- sd(unstandardizedSNARC$unstdSNARC.slope)
ext_SD_unstd_slope

SEMunstd <- SD_unstd_slope * (sqrt((1-unstd_rel)))

## Calculating confidence intervals
# 80%
unstandardizedSNARC$ci.80.low <-unstandardizedSNARC$unstdSNARC.slope - 1.29*SEMunstd 
unstandardizedSNARC$ci.80.hi <-unstandardizedSNARC$unstdSNARC.slope + 1.29*SEMunstd

# 90%
unstandardizedSNARC$ci.90.low <-unstandardizedSNARC$unstdSNARC.slope - 1.65*SEMunstd 
unstandardizedSNARC$ci.90.hi <-unstandardizedSNARC$unstdSNARC.slope + 1.65*SEMunstd

# 95%
unstandardizedSNARC$ci.95.low <-unstandardizedSNARC$unstdSNARC.slope - 1.96*SEMunstd 
unstandardizedSNARC$ci.95.hi <-unstandardizedSNARC$unstdSNARC.slope + 1.96*SEMunstd

# 99%
unstandardizedSNARC$ci.99.low <-unstandardizedSNARC$unstdSNARC.slope - 2.58*SEMunstd 
unstandardizedSNARC$ci.99.hi <-unstandardizedSNARC$unstdSNARC.slope + 2.58*SEMunstd


# Saving the data file with unstandardized SNARC slopes and respective confidence intervals

# write.table(unstandardizedSNARC, file="Unstandardized_SNARC_with_c_i.txt", sep = "\t", dec=".", quote = FALSE, row.names=FALSE)
plot(sort(unstandardizedSNARC$unstdSNARC.slope), pch = 19, ylim = c(-40, 30))
abline(h = 0)
ind <- order(unstandardizedSNARC$unstdSNARC.slope)
plotCI(x = 1: nrow(unstandardizedSNARC)
       , y = unstandardizedSNARC$unstdSNARC.slope[ind]
       , ui= unstandardizedSNARC$ci.95.hi[ind]
       , li= unstandardizedSNARC$ci.95.low[ind]
       , add = TRUE
       , col = "gray40"
       , pch = 21
       , pt.bg= "gray40"
       , cex = 1.5)

cis <- apply(out$theta, 2, quantile, probs =  c(.025, .975))

plot(sort(-testSNARC$unstdSNARC.slope), pch = 1, ylim = c(-30, 30))
abline(h = 0)
ind <- order(-testSNARC$unstdSNARC.slope)
plotCI(x = 1: nrow(testSNARC)
       , y = out$ind.effects[ind]
       , ui= cis[2, ind]
       , li= cis[1, ind]
       , add = TRUE
       , col = "gray40"
       , pch = 21
       , pt.bg= "gray40"
       , cex = 1.5)

## Standardized SNARC
# Type in unstandardized SNARC reliability
# Note that reliability can be calculated using other scripts provided
# Currently the value is the one obtained from the sample data file

std_rel <- .65

# If you want to provide the SD slope value based some external data (e.g. from other experiment in which more people were tested), please type it below, otherwise SD slope from the file being analysed.

ext_SD_std_slope <- "place the numerical value without a quotation mark instead of this text"
if(is.numeric(ext_SD_std_slope)) SD_std_slope <- ext_SD_std_slope  else SD_std_slope <- sd(standardizedSNARC$standardized.SNARC)
SD_std_slope

SEMstd <- SD_std_slope * (sqrt((1-std_rel)))

## Calculating confidence intervals
# 80%
standardizedSNARC$ci.80.low <-standardizedSNARC$standardized.SNARC - 1.29*SEMstd 
standardizedSNARC$ci.80.hi <-standardizedSNARC$standardized.SNARC + 1.29*SEMstd

# 90%
standardizedSNARC$ci.90.low <-standardizedSNARC$standardized.SNARC - 1.65*SEMstd 
standardizedSNARC$ci.90.hi <-standardizedSNARC$standardized.SNARC + 1.65*SEMstd

# 95%
standardizedSNARC$ci.95.low <-standardizedSNARC$standardized.SNARC - 1.96*SEMstd 
standardizedSNARC$ci.95.hi <-standardizedSNARC$standardized.SNARC + 1.96*SEMstd

# 99%
standardizedSNARC$ci.99.low <-standardizedSNARC$standardized.SNARC - 2.58*SEMstd 
standardizedSNARC$ci.99.hi <-standardizedSNARC$standardized.SNARC + 2.58*SEMstd

# Saving the data file with standardized SNARC slopes and respective confidence intervals
# write.table(standardizedSNARC, file="Standardized_SNARC_with_c_i.txt", sep = "\t", dec=".", quote = FALSE, row.names=FALSE)
plot(sort(standardizedSNARC$standardized.SNARC), pch = 19, ylim = c(-3, 1.5))
abline(h = 0)
ind <- order(standardizedSNARC$standardized.SNARC)
plotCI(x = 1: nrow(standardizedSNARC)
       , y = standardizedSNARC$standardized.SNARC[ind]
       , ui= standardizedSNARC$ci.95.hi[ind]
       , li= standardizedSNARC$ci.95.low[ind]
       , add = TRUE
       , col = "gray40"
       , pch = 21
       , pt.bg= "gray40"
       , cex = 1.5)

### Calculating summary tables for prevalence of the SNARC under different confidence levels
## Unstandardized SNARC
unstandardizedSNARC$ci.80.consistent.neg <- ifelse(unstandardizedSNARC$ci.80.hi < 0, 1, 0)
unstandardizedSNARC$ci.80.consistent.pos <- ifelse(unstandardizedSNARC$ci.80.low > 0, 1, 0)
unstandardizedSNARC$ci.80.inconsistent <- ifelse(unstandardizedSNARC$ci.80.low < 0 & unstandardizedSNARC$ci.80.hi > 0, 1, 0)
unstandardizedSNARC$ci.80.inconsistent.less.than.0 <- ifelse(unstandardizedSNARC$unstdSNARC.slope < 0 & unstandardizedSNARC$ci.80.inconsistent == 1, 1, 0)
unstd_inconsistent_80ci_less_than_zero_prop <- mean(subset(unstandardizedSNARC, ci.80.inconsistent == 1)$ci.80.inconsistent.less.than.0)

unstandardizedSNARC$ci.90.consistent.neg <- ifelse(unstandardizedSNARC$ci.90.hi < 0, 1, 0)
unstandardizedSNARC$ci.90.consistent.pos <- ifelse(unstandardizedSNARC$ci.90.low > 0, 1, 0)
unstandardizedSNARC$ci.90.inconsistent<- ifelse(unstandardizedSNARC$ci.90.low < 0 & unstandardizedSNARC$ci.90.hi > 0, 1, 0)
unstandardizedSNARC$ci.90.inconsistent.less.than.0 <- ifelse(unstandardizedSNARC$unstdSNARC.slope < 0 & unstandardizedSNARC$ci.90.inconsistent == 1, 1, 0)
unstd_inconsistent_90ci_less_than_zero_prop <- mean(subset(unstandardizedSNARC, ci.90.inconsistent == 1)$ci.90.inconsistent.less.than.0)

unstandardizedSNARC$ci.95.consistent.neg <- ifelse(unstandardizedSNARC$ci.95.hi < 0, 1, 0)
unstandardizedSNARC$ci.95.consistent.pos <- ifelse(unstandardizedSNARC$ci.95.low > 0, 1, 0)
unstandardizedSNARC$ci.95.inconsistent<- ifelse(unstandardizedSNARC$ci.95.low < 0 & unstandardizedSNARC$ci.95.hi > 0, 1, 0)
unstandardizedSNARC$ci.95.inconsistent.less.than.0 <- ifelse(unstandardizedSNARC$unstdSNARC.slope < 0 & unstandardizedSNARC$ci.95.inconsistent == 1, 1, 0)
unstd_inconsistent_95ci_less_than_zero_prop <- mean(subset(unstandardizedSNARC, ci.95.inconsistent == 1)$ci.95.inconsistent.less.than.0)

unstandardizedSNARC$ci.99.consistent.neg <- ifelse(unstandardizedSNARC$ci.99.hi < 0, 1, 0)
unstandardizedSNARC$ci.99.consistent.pos <- ifelse(unstandardizedSNARC$ci.99.low > 0, 1, 0)
unstandardizedSNARC$ci.99.inconsistent<- ifelse(unstandardizedSNARC$ci.99.low < 0 & unstandardizedSNARC$ci.99.hi > 0, 1, 0)
unstandardizedSNARC$ci.99.inconsistent.less.than.0 <- ifelse(unstandardizedSNARC$unstdSNARC.slope < 0 & unstandardizedSNARC$ci.99.inconsistent == 1, 1, 0)
unstd_inconsistent_99ci_less_than_zero_prop <- mean(subset(unstandardizedSNARC, ci.99.inconsistent == 1)$ci.99.inconsistent.less.than.0)

summary_table_snarc_unstd <- as.data.frame(cbind((rbind("ci.80.consistent.neg", "ci.80.inconsistent", "ci.80.consistent.pos","80.inconsistent.less.than0","ci.90.consistent.neg", "ci.90.inconsistent", "ci.90.consistent.pos","90.inconsistent.less.than0", "ci.95.consistent.neg", "ci.95.inconsistent", "ci.95.consistent.pos", "95.inconsistent.less.than0","ci.99.consistent.neg", "ci.99.inconsistent", "ci.99.consistent.pos", "99.inconsistent.less.than0")), (rbind(mean(unstandardizedSNARC$ci.80.consistent.neg), mean(unstandardizedSNARC$ci.80.inconsistent), mean(unstandardizedSNARC$ci.80.consistent.pos),unstd_inconsistent_80ci_less_than_zero_prop, mean(unstandardizedSNARC$ci.90.consistent.neg), mean(unstandardizedSNARC$ci.90.inconsistent), mean(unstandardizedSNARC$ci.90.consistent.pos),unstd_inconsistent_90ci_less_than_zero_prop, mean(unstandardizedSNARC$ci.95.consistent.neg), mean(unstandardizedSNARC$ci.95.inconsistent), mean(unstandardizedSNARC$ci.95.consistent.pos),unstd_inconsistent_95ci_less_than_zero_prop,mean(unstandardizedSNARC$ci.99.consistent.neg), mean(unstandardizedSNARC$ci.99.inconsistent), mean(unstandardizedSNARC$ci.99.consistent.pos), unstd_inconsistent_99ci_less_than_zero_prop))))
names(summary_table_snarc_unstd) <- c("SNARC.type", "proportion.of.participants")

# unstandardized snarc
write.table(summary_table_snarc_unstd, file="unstandardized_snarc_proportions_of_participants_final.txt", sep = "\t", dec=".", quote = FALSE, row.names=FALSE)


## Standardized SNARC
standardizedSNARC$ci.80.consistent.neg <- ifelse(standardizedSNARC$ci.80.hi < 0, 1, 0)
standardizedSNARC$ci.80.consistent.pos <- ifelse(standardizedSNARC$ci.80.low > 0, 1, 0)
standardizedSNARC$ci.80.inconsistent <- ifelse(standardizedSNARC$ci.80.low < 0 & standardizedSNARC$ci.80.hi > 0, 1, 0)
standardizedSNARC$ci.80.inconsistent.less.than.0 <- ifelse(standardizedSNARC$standardized.SNARC < 0 & standardizedSNARC$ci.80.inconsistent == 1, 1, 0)
std_inconsistent_80ci_less_than_zero_prop <- mean(subset(standardizedSNARC, ci.80.inconsistent == 1)$ci.80.inconsistent.less.than.0)

standardizedSNARC$ci.90.consistent.neg <- ifelse(standardizedSNARC$ci.90.hi < 0, 1, 0)
standardizedSNARC$ci.90.consistent.pos <- ifelse(standardizedSNARC$ci.90.low > 0, 1, 0)
standardizedSNARC$ci.90.inconsistent<- ifelse(standardizedSNARC$ci.90.low < 0 & standardizedSNARC$ci.90.hi > 0, 1, 0)
standardizedSNARC$ci.90.inconsistent.less.than.0 <- ifelse(standardizedSNARC$standardized.SNARC < 0 & standardizedSNARC$ci.90.inconsistent == 1, 1, 0)
std_inconsistent_90ci_less_than_zero_prop <- mean(subset(standardizedSNARC, ci.90.inconsistent == 1)$ci.90.inconsistent.less.than.0)

standardizedSNARC$ci.95.consistent.neg <- ifelse(standardizedSNARC$ci.95.hi < 0, 1, 0)
standardizedSNARC$ci.95.consistent.pos <- ifelse(standardizedSNARC$ci.95.low > 0, 1, 0)
standardizedSNARC$ci.95.inconsistent<- ifelse(standardizedSNARC$ci.95.low < 0 & standardizedSNARC$ci.95.hi > 0, 1, 0)
standardizedSNARC$ci.95.inconsistent.less.than.0 <- ifelse(standardizedSNARC$standardized.SNARC < 0 & standardizedSNARC$ci.95.inconsistent == 1, 1, 0)
std_inconsistent_95ci_less_than_zero_prop <- mean(subset(standardizedSNARC, ci.95.inconsistent == 1)$ci.95.inconsistent.less.than.0)

standardizedSNARC$ci.99.consistent.neg <- ifelse(standardizedSNARC$ci.99.hi < 0, 1, 0)
standardizedSNARC$ci.99.consistent.pos <- ifelse(standardizedSNARC$ci.99.low > 0, 1, 0)
standardizedSNARC$ci.99.inconsistent<- ifelse(standardizedSNARC$ci.99.low < 0 & standardizedSNARC$ci.99.hi > 0, 1, 0)
standardizedSNARC$ci.99.inconsistent.less.than.0 <- ifelse(standardizedSNARC$standardized.SNARC < 0 & standardizedSNARC$ci.99.inconsistent == 1, 1, 0)
std_inconsistent_99ci_less_than_zero_prop <- mean(subset(standardizedSNARC, ci.99.inconsistent == 1)$ci.99.inconsistent.less.than.0)

summary_table_snarc_std <- as.data.frame(cbind((rbind("ci.80.consistent.neg", "ci.80.inconsistent", "ci.80.consistent.pos","80.inconsistent.less.than0","ci.90.consistent.neg", "ci.90.inconsistent", "ci.90.consistent.pos","90.inconsistent.less.than0", "ci.95.consistent.neg", "ci.95.inconsistent", "ci.95.consistent.pos", "95.inconsistent.less.than0","ci.99.consistent.neg", "ci.99.inconsistent", "ci.99.consistent.pos", "99.inconsistent.less.than0")), (rbind(mean(standardizedSNARC$ci.80.consistent.neg), mean(standardizedSNARC$ci.80.inconsistent), mean(standardizedSNARC$ci.80.consistent.pos),std_inconsistent_80ci_less_than_zero_prop,mean(standardizedSNARC$ci.90.consistent.neg), mean(standardizedSNARC$ci.90.inconsistent), mean(standardizedSNARC$ci.90.consistent.pos),std_inconsistent_90ci_less_than_zero_prop, mean(standardizedSNARC$ci.95.consistent.neg), mean(standardizedSNARC$ci.95.inconsistent), mean(standardizedSNARC$ci.95.consistent.pos),std_inconsistent_95ci_less_than_zero_prop,mean(standardizedSNARC$ci.99.consistent.neg), mean(standardizedSNARC$ci.99.inconsistent), mean(standardizedSNARC$ci.99.consistent.pos),std_inconsistent_99ci_less_than_zero_prop))))
names(summary_table_snarc_std) <- c("SNARC.type", "proportion.of.participants")

# standardized snarc
write.table(summary_table_snarc_std, file="standardized_snarc_proportions_of_participants_final.txt", sep = "\t", dec=".", quote = FALSE, row.names=FALSE)
