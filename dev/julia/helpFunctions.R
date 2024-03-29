prior.p.greater <- function(M, I, a = alpha, b = beta, rscale){
  
  s2 <- MCMCpack::rinvgamma(M, a, b)
  var_mu <- MCMCpack::rinvgamma(M, 1/2, 1/2 * rscale^2)
  mu <- rnorm(M, 0, sqrt(var_mu))
  res <- exp(pnorm(0, mu, sqrt(s2), lower.tail = F, log.p = T) * I)
  
  return(mean(res))
}

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
    if (side[r] == h) {
      X.full[r, I + 1] <- 1
      X.full[r, I + 1 + id[r]] <- 1}
    else{
        X.full[r, I + 1] <- -1/2
        X.full[r, I + 1 + id[r]] <- -1/2
        }
    # Effect of number
    X.full[r, 2*I + 2] <- number[r]
    X.full[r, 2*I + 2 + id[r]] <- number[r]
    ## interaction between number and side, i.e. SNARC effect
    if (side[r] == h) {
      X.full[r, 3*I + 3] <- number[r]
      X.full[r, 3*I + 3 + id[r]] <- number[r]}
    else{
        X.full[r, 3*I + 3] <- -number[r]/2
        X.full[r, 3*I + 3 + id[r]] <- -number[r]/2
      }
  }
  
  gMap.full <- c(rep(0, I), 1, rep(2, I), 3, rep(4, I), 5, rep(6, I))
  
  X.one <- matrix(nrow = R, ncol = 3 * I + 3, 0)
  for (r in 1:R){
    # baseline RT
    X.one[r, id[r]] <- 1
    # effect of side
    if (side[r] == h) {
      X.one[r, I + 1] <- 1/2
      X.one[r, I + 1 + id[r]] <- 1/2}
    else{
      X.one[r, I + 1] <- -1/2
      X.one[r, I + 1 + id[r]] <- -1/2
        }
    # Effect of number
    X.one[r, 2*I + 2] <- number[r]
    X.one[r, 2*I + 2 + id[r]] <- number[r]
    ## interaction between number and side, i.e. SNARC effect
    if (side[r] == h) {
      X.one[r, 3*I + 3] <- number[r]/2
    }
    else{
      X.one[r, 3*I + 3] <- -number[r]/2
      }
  }
  
  gMap.one <- c(rep(0, I), 1, rep(2, I), 3, rep(4, I), 5)
  
  
  X.null <- matrix(nrow = R, ncol = 3*I + 2, 0)
  for (r in 1:R){
    # baseline RT
    X.null[r, id[r]] <- 1
    # effect of side
    if (side[r] == h) {
      X.null[r, I + 1] <- 1/2
      X.null[r, I + 1 + id[r]] <- 1/2}
    else{
        X.null[r, I + 1] <- -1/2
        X.null[r, I + 1 + id[r]] <- -1/2
    }
    # Effect of number
    X.null[r, 2*I + 2] <- number[r]
    X.null[r, 2*I + 2 + id[r]] <- number[r]
  }
  
  gMap.null <- c(rep(0, I), 1, rep(2, I), 3, rep(4, I))
  
  return(list(X.full = X.full
              , gMap.full = gMap.full
              , X.one = X.one
              , gMap.one = gMap.one
              , X.null = X.null
              , gMap.null = gMap.null
              , R = R
              , I= I
              , sub = id))
}

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
  bf.one <- BayesFactor::nWayAOV(y
                                 , prep$X.one
                                 , prep$gMap.one
                                 , rscale = c(1, .5, .5, .5, .5, meanScale)
                                 , posterior = F
                                 , iterations = iter)
  #posterior estimates from the common-effect model (but with unconstrained common effect, i.e. can be positive or negative)
  mcmc.one <- BayesFactor::nWayAOV(y
                                   , prep$X.one
                                   , prep$gMap.one
                                   , rscale = c(1, .5, .5, .5, .5, meanScale)
                                   , posterior = T
                                   , iterations = iter)
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
  
  ## Evaluate positive common effect model
  good0 <- mean(mcmc.one[keep, i.theta0] < 0)
  
  #Calculate Bayes factors
  bf.FP <- PriorCount/PostCount
  bf.F0 <- exp(bf.full$bf - bf.null$bf)
  bf.F1 <- exp(bf.full$bf - bf.one$bf) / 2
  bf.10 <- exp(bf.one$bf - bf.null$bf) * 2
  
  #Return posterior estimates from the unconstrained model
  m <- colMeans(myTheta)
  new.sd <- sd(m)
  new.mean <- mean(mcmc.full[keep, i.theta0])
  # i.g3 <- 2 * prep$I + 6
  # new.pop.sd <- mean(mcmc.full[keep, i.g3])
  return(list(ind.effects = m
              , posterior.mean = new.mean, posterior.sd = new.sd
              , bfs = c(bf.pu =  1/ bf.FP, bf.1u = 1/ bf.F1, bf.0u = 1 / bf.F0)
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