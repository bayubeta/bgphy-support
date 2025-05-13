# code for simulation study part 1: parameter estimation,
# for the paper Bayesian inference of mixed Gaussian phylogenetic models (2025)
# by Brahmantio, Bartoszek, Yapar.

# ============== load required libraries
library(PCMBase)
library(PCMBaseCpp)
library(PCMFit)
library(bgphy)



# =========================== Bayesian inference using bgphy ===================
# ------------- set model
OUMixed <- setModel(tree = lizardTree, regime_names = c("Ancestral", "New"),
                    modeltypes = c("OU", "OU"), startNodes = list(Ancestral = 101, New = 135))
# ------------- set priors
OUMixed$priors$X0 <- OUMixed$priors$theta_1 <- OUMixed$priors$theta_2 <- prior_normal(mean = 0, sd = 10)
OUMixed$priors$alpha_1 <- OUMixed$priors$alpha_2 <- OUMixed$priors$sigma_1 <- OUMixed$priors$sigma_2 <- prior_halfnormal(sigma = 10)

# ------------- prepare a list of model
bgphy_models <- list()

# ------------- run posterior inference 100 times
# --------> one iteration takes ~4 mins on 
# AMD Ryzen 5 PRO 4650U 2100 Mhz, 6 Core(s), 12 Logical Processor(s)
# 16 GB RAM

set.seed(1234)
for (i in 1:100){
  # define X
  X <- matrix(XMixedOU[i,], nrow = 1, dimnames = list(NULL, colnames(XMixedOU)))
  timestamp()
  bgphy_models[[i]] <- bgphy(OUMixed, X)
  print(i)
}

# save(bgphy_models, file = "bgphy_models-bkup.RData")
# load(file = "bgphy_models.RData")


# =========================== ML inference using PCMFit ========================
# ------------- set model: use OUMixed

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
MLres <- list()


set.seed(1234)
# for each dataset, find the model with the highest likelihood
# --------> one iteration takes ~50 secs 
for (i in 1:100){
  # reset model list
  modelset <- list()
  
  # define the data
  X <- matrix(XMixedOU[i,], nrow = 1, dimnames = list(NULL, colnames(XMixedOU)))
  timestamp()
  for (j in 1:50){
    modelset[[j]] <- PCMFit(X, OUMixed$tree, OUMixed$model, metaI = PCMBaseCpp::PCMInfoCpp)
  }
  
  MLres[[i]] <- modelset[[which.max(sapply(modelset, function(x){x$logLikOptim[1]}))]]
  
  print(i)
}

# save(MLres, file = "MLres.RData")
# load(file = "MLres.RData")



# ====================== Boxplots: posterior median vs ML ======================

parnames <- names(attr(bgphy_models[[1]], "model")$priors)
bgphy_quantiles <- list()
for (name in parnames){
  bgphy_quantiles[[name]] <- t(sapply(bgphy_models, function(x){x$quantiles[name,]}))
}
truepars <- c(0, 2, 2, 1, 5, 0, 0.5)
positive <- c(0,1,0,1,1,0,1)

dwplot <- function(quantiles, main, truepar, positive = FALSE, ...){
  maxreach <- max(abs(min(quantiles)), abs(max(quantiles)))
  
  if (positive){
    xlim_min <- 0
  }else{
    xlim_min <- truepar - maxreach
  }
  
  plot(rep(0, 100), 1:100, xlim = c(xlim_min, truepar + maxreach),
       pch = NA, xlab = "", yaxt = "n", ylab = "", main = main)
  grid()
  lines(c(truepar,truepar), c(0,100), col = "red", lty = 2, lwd = 2)
  for (i in 1:100){
    lines(c(quantiles[i,1],quantiles[i,3]), col = rgb(0.1,0.1,0.7, 0.4), c(i,i),
          lwd = 2)
    points(quantiles[i,2], i, pch = 16, col = rgb(0.1,0.1,0.7, 0.4))
  }
}


ML_pars <- matrix(nrow = 100, ncol = length(parnames))
colnames(ML_pars) <- parnames




for (i in 1:100){
  ML_pars[i,] <- PCMBase::PCMParamGetShortVector(MLres[[i]]$modelOptim)
}
bgphy_medians <- sapply(bgphy_quantiles, function(x){x[,2]})


M <- matrix(c(1:7,rep(8,7)), nrow = 2, ncol = 7, byrow = TRUE)


{
  pdf(file = "Fig5.pdf", height = 8, width = 15)
  layout(mat = rbind(M, M + max(M)), heights = c(8, 2, 8, 1))
  i <- 1
  par(mar = c(1,4,4,1))
  par(cex.main = 2)
  par(cex.axis = 1.5)
  boxplot(as.data.frame(cbind(bgphy_medians[,"X0"], ML_pars[,"X0"])), xaxt = "n", main = bquote(X[0]), col = c("green2", "skyblue2"), boxwex = 0.5)
  lines(c(0.5, 2.5), c(truepars[i], truepars[i]), col = "red", lty = 2, lwd = 2); i <- i+1
  mtext("(a)", side = 3, adj = 0, line = 1, font = 2)
  boxplot(as.data.frame(cbind(bgphy_medians[,"alpha_1"], ML_pars[,"alpha_1"])), xaxt = "n", main = bquote(alpha[1]), col = c("green2", "skyblue2"), boxwex = 0.5)
  lines(c(0.5, 2.5), c(truepars[i], truepars[i]), col = "red", lty = 2, lwd = 2); i <- i+1
  boxplot(as.data.frame(cbind(bgphy_medians[,"theta_1"], ML_pars[,"theta_1"])), xaxt = "n", main = bquote(theta[1]), col = c("green2", "skyblue2"), boxwex = 0.5)
  lines(c(0.5, 2.5), c(truepars[i], truepars[i]), col = "red", lty = 2, lwd = 2); i <- i+1
  boxplot(as.data.frame(cbind(bgphy_medians[,"sigma_1"], ML_pars[,"sigma_1"])), xaxt = "n", main = bquote(sigma[1]), col = c("green2", "skyblue2"), boxwex = 0.5)
  lines(c(0.5, 2.5), c(truepars[i], truepars[i]), col = "red", lty = 2, lwd = 2); i <- i+1
  boxplot(as.data.frame(cbind(bgphy_medians[,"alpha_2"], ML_pars[,"alpha_2"])), xaxt = "n", main = bquote(alpha[2]), col = c("green2", "skyblue2"), boxwex = 0.5)
  lines(c(0.5, 2.5), c(truepars[i], truepars[i]), col = "red", lty = 2, lwd = 2); i <- i+1
  boxplot(as.data.frame(cbind(bgphy_medians[,"theta_2"], ML_pars[,"theta_2"])), xaxt = "n", main = bquote(theta[2]), col = c("green2", "skyblue2"), boxwex = 0.5)
  lines(c(0.5, 2.5), c(truepars[i], truepars[i]), col = "red", lty = 2, lwd = 2); i <- i+1
  boxplot(as.data.frame(cbind(bgphy_medians[,"sigma_2"], ML_pars[,"sigma_2"])), xaxt = "n", main = bquote(sigma[2]), col = c("green2", "skyblue2"), boxwex = 0.5)
  lines(c(0.5, 2.5), c(truepars[i], truepars[i]), col = "red", lty = 2, lwd = 2)
  par(mar = c(2,0,0,0))
  plot(NA, xlim = c(-1,1), ylim = c(-1,1), xaxt = 'n', yaxt = 'n', bty = 'n',
       xlab = '', ylab = '')
  legend("center", legend = c("True value", "Posterior median", "Maximum likelihood"),
         lty = c(2, NA, NA), lwd = c(3,3,3), pch = c(NA, 15, 15),
         col = c("red", "green2", "skyblue3"), bty = "n", horiz = TRUE, cex = 2,
         x.intersp = c(1,1,0), text.width = c(0.3,0.3,0))
  

# =========================== Confidence intervals plot ========================
# confidence interval plots

  par(mar = c(2,4,4,1))
  dwplot(bgphy_quantiles$X0, bquote(X[0]), 0, positive = 0)
  mtext("(b)", side = 3, adj = 0, line = 1, font = 2)
  dwplot(bgphy_quantiles$alpha_1, bquote(alpha[1]), 2, positive = 1)
  dwplot(bgphy_quantiles$theta_1, bquote(theta[1]), 2, positive = 0)
  dwplot(bgphy_quantiles$sigma_1, bquote(sigma[1]), 1, positive = 1)
  dwplot(bgphy_quantiles$alpha_2, bquote(alpha[2]), 5, positive = 1)
  dwplot(bgphy_quantiles$theta_2, bquote(theta[2]), 0, positive = 0)
  dwplot(bgphy_quantiles$sigma_2, bquote(sigma[2]), 0.5, positive = 1)
  par(mar = c(0,0,0,0))
  plot(NA, xlim = c(-1,1), ylim = c(-1,1), xaxt = 'n', yaxt = 'n', bty = 'n',
       xlab = '', ylab = '')
  legend("center", legend = c("True value", "95% Marginal CI", "Median"),
         lty = c(2, 1, NA), lwd = c(2, 2, NA), pch = c(NA, NA, 16),
         col = c("red", rgb(0.1,0.1,0.7, 0.7), rgb(0.1,0.1,0.7, 0.7)), bty = "n",
         horiz = TRUE, cex = 2, x.intersp = c(1,1,0), text.width = c(0.3,0.3,0))
  dev.off()
}

