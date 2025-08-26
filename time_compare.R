# code for computational time comparison,
# for the paper Bayesian inference of mixed Gaussian phylogenetic models (2025)
# by Brahmantio, Bartoszek, Yapar.

# ============== load required libraries
library(PCMBase)
library(PCMBaseCpp)
library(PCMFit)
library(bgphy)
library(bayou)


# ----------------------------------- bgphy ------------------------------------
OUMixed <- setModel(tree = lizardTree, regime_names = c("Ancestral", "New"),
                    modeltypes = c("OU", "OU"), 
                    startNodes = list(Ancestral = 101, New = 135))
# ------------- set priors
OUMixed$priors$X0 <- OUMixed$priors$theta_1 <- 
  OUMixed$priors$theta_2 <- prior_normal(mean = 0, sd = 10)
OUMixed$priors$alpha_1 <- OUMixed$priors$alpha_2 <- OUMixed$priors$sigma_1 <- 
  OUMixed$priors$sigma_2 <- prior_halfnormal(sigma = 10)

# define data
X <- matrix(XMixedOU[1,], nrow = 1, dimnames = list(NULL, colnames(XMixedOU)))

bgphy_times <- c()
# run 10 times

for (i in 1:10){
  timestamp()
  start_time <- Sys.time()
  bgphy(OUMixed, X, nsample = 10000)
  bgphy_times[i] <- Sys.time() - start_time
}

print(paste0("mean: ", round(mean(bgphy_times)*60, 4), " secs."))
print(paste0("sd: ", round(sd(bgphy_times)*60, 4), " secs."))

# ----------------------------------- bayou ------------------------------------
# set priors

# sanity check
priorsOUMixed <- make.prior(lizardTree,
                            dists = list(dalpha = "dhalfcauchy", 
                                         dsig2 = "dhalfcauchy",
                                         dk = "fixed", dsb="fixed", 
                                         dtheta = "dnorm"),
                            param = list(dalpha = list(scale=10), 
                                         dsig2 = list(scale=10),
                                         dk = "fixed", dsb="fixed",
                                         dtheta = list(mean=0, sd=10)),
                            fixed = list(k = 1, sb = c(131), 
                                        ntheta = 2, loc = c(0), t2 = 2))

X_bayou <- setNames(XMixedOU[1,], colnames(XMixedOU))
# define model

# run 10 times
bayou_times <- c()
for (i in 1:10){
  start_time <- Sys.time()
  bmcmc <- bayou.makeMCMC(tree = lizardTree, dat = X_bayou, 
                          prior = priorsOUMixed,
                          plot.freq = NULL, samp = 10,
                          ticker.freq = 10000)
  timestamp()
  bmcmc$run(100000)
  
  
  chain <- bmcmc$load()
  Bk <- qbeta(seq(0,1, length.out=5), 0.3,1)
  sstone <- bmcmc$steppingstone(10000, chain, Bk)
  bayou_times[i] <- Sys.time() - start_time
}

print(paste0("mean: ", round(mean(bayou_times)*60, 4), " secs."))
print(paste0("sd: ", round(sd(bayou_times)*60, 4), " secs."))


# ----------------------------- Maximum likelihood -----------------------------
# -------------------- define upper and lower limits --------------------
PCMParamLowerLimit.OU <- function(o, k, R, ...) {
  o <- NextMethod()
  k <- attr(o, "k", exact = TRUE)
  R <- length(attr(o, "regimes", exact = TRUE))
  
  # Lower limit for Theta
  o$Theta[1] = -100
  
  # Lower limit for H
  o$H[,,1] = 0
  
  # Lower limits for Sigma_x
  o$Sigma_x[,,1] = 0
  o
}
PCMParamUpperLimit.OU <- function(o, k, R, ...) {
  o <- NextMethod()
  k <- attr(o, "k", exact = TRUE)
  R <- length(attr(o, "regimes", exact = TRUE))
  
  # Upper limit for Theta
  o$Theta[1] = 100
  
  # Upper limit for H
  o$H[,,1] = 100
  
  # Upper limits for Sigma_x
  o$Sigma_x[,,1] = 100
  o
}


# ------------- prepare a list of results
ML_times <- c()


set.seed(1234)
# for each dataset, find the model with the highest likelihood
# --------> one iteration takes ~50 secs 
for (i in 1:10){
  start_time <- Sys.time()
  # define the data
  X <- matrix(XMixedOU[1,], nrow = 1, dimnames = list(NULL, colnames(XMixedOU)))
  timestamp()
  PCMFit(X, OUMixed$tree, OUMixed$model, metaI = PCMBaseCpp::PCMInfoCpp)
  # in seconds
  ML_times[i] <- Sys.time() - start_time
}

print(paste0("mean: ", round(mean(ML_times), 4), " secs."))
print(paste0("sd: ", round(sd(ML_times), 4), " secs."))
