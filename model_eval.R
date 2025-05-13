# code for simulation study part 1: parameter estimation,
# for the paper Bayesian inference of mixed Gaussian phylogenetic models (2025)
# by Brahmantio, Bartoszek, Yapar.

# ============== load required libraries
library(PCMBase)
library(PCMBaseCpp)
library(PCMFit)
library(bgphy)



# ============= define model and simulate data using PCMBase ===================
# ----------------------- Global model
# already included in the bgphy package, 
# see bgphy::modelOU for the model, bgphy::XOU for the simulated data


# ---------------- Mixed OU, 2 Regimes
# list of models in each regimes
modeltypes2 <- setNames(c("OU", "OU"), c("ancestral", "R1"))

Model2 <- setModel(tree = lizardTree, 
                   regime_names = c("ancestral", "R1"),
                   modeltypes = modeltypes2, 
                   startNodes = list(ancestral = 101, R1 = 134))

# change the PCMBase model's parameter for simulation
M2 <- Model2$model
M2$ancestral$H[,,1] <- 2
M2$ancestral$Theta[1,] <- 2
M2$ancestral$Sigma_x[,,1] <- 2
M2$R1$H[,,1] <- 1
M2$R1$Theta[1,] <- -2
M2$R1$Sigma_x[,,1] <- 1

# simulate
set.seed(1234)
X2 <- t(replicate(100, PCMBase::PCMSim(Model2$tree, M2, M2$X0)[1,1:100]))
save(X2, file = "X2.RData")


# ---------------- Mixed OU, 4 Regimes
# list of models in each regimes
regime_names <- c("ancestral", "R1", "R2", "R3")
modeltypes <- setNames(c("OU", "OU", "OU", "OU"), regime_names)
startNodes <- setNames(as.list(c(101, 134, 160, 103)), regime_names)


Model3 <- setModel(tree = lizardTree, regime_names = regime_names, 
                   modeltypes = modeltypes, startNodes = startNodes)

# change the PCMBase model's parameter for simulation
M3 <- Model3$model
M3$ancestral$H[,,1] <- 3
M3$ancestral$Theta[1,] <- 2
M3$ancestral$Sigma_x[,,1] <- 2
# ----------
M3$R1$H[,,1] <- 2
M3$R1$Theta[1,] <- 1
M3$R1$Sigma_x[,,1] <- 1
# ----------
M3$R2$H[,,1] <- 2
M3$R2$Theta[1,] <- -1
M3$R2$Sigma_x[,,1] <- 1
# ----------
M3$R3$H[,,1] <- 3
M3$R3$Theta[1,] <- -2
M3$R3$Sigma_x[,,1] <- 2

# simulate
set.seed(1234)
X3 <- t(replicate(100, PCMBase::PCMSim(Model3$tree, M3, M3$X0)[1,1:100]))
save(X3, file = "X3.RData")



# ========================= Run inference using bgphy ==========================
# ------------------- set priors for all models
# more informative, but still weak:
# - unbounded = Normal(0,5)
# - bounded = Half-normal(5)

# Global
Model1 <- setModel(tree = lizardTree, 
                   regime_names = "ancestral", modeltypes = "OU")
Model1$priors$X0 <- Model1$priors$theta <- prior_normal(mean = 0, sd = 5)
Model1$priors$alpha <- Model1$priors$sigma <- prior_halfnormal(sigma = 5)
save(Model1, file = "Model1.RData")

# 2 Regimes
Model2$priors$X0 <- Model2$priors$theta_1 <- Model2$priors$theta_2 <- prior_normal(mean = 0, sd = 5)
Model2$priors$alpha_1 <- Model2$priors$alpha_2 <- prior_halfnormal(sigma = 5)
Model2$priors$sigma_1 <- Model2$priors$sigma_2 <- prior_halfnormal(sigma = 5)
save(Model2, file = "Model2.RData")

# 4 regimes
Model3$priors$X0 <- Model3$priors$theta_1 <- Model3$priors$theta_2 <- prior_normal(mean = 0, sd = 5)
Model3$priors$theta_3 <- Model3$priors$theta_4 <- prior_normal(mean = 0, sd = 5)
Model3$priors$alpha_1 <- Model3$priors$alpha_2 <- prior_halfnormal(sigma = 5)
Model3$priors$sigma_1 <- Model3$priors$sigma_2 <- prior_halfnormal(sigma = 5)
Model3$priors$alpha_3 <- Model3$priors$alpha_4 <- prior_halfnormal(sigma = 5)
Model3$priors$sigma_3 <- Model3$priors$sigma_4 <- prior_halfnormal(sigma = 5)
save(Model3, file = "Model3.RData")


# ------------------ run posterior inferences
# load models
load(file = "Model1.RData")
load(file = "Model2.RData")
load(file = "Model3.RData")

set.seed(1234)

# -------------- model 1 --------------
modelevX1 <- list(mod1 = list(), mod2 = list(), mod3 = list())
# set data
XData <- XOU

for (i in 1:100){
  # prepare data
  X <- matrix(XData[i,], nrow = 1, dimnames = list(NULL, colnames(XData)))
  
  timestamp()
  # model1
  modelevX1$mod1[[i]] <- bgphy(model = Model1, X = X)
  # model2
  modelevX1$mod2[[i]] <- bgphy(model = Model2, X = X)
  # model3
  modelevX1$mod3[[i]] <- bgphy(model = Model3, X = X)
  
  print(i)
}

save(modelevX1, file = "modelevX1.RData")




set.seed(1234)
# -------------- model 2 --------------
modelevX2 <- list(mod1 = list(), mod2 = list(), mod3 = list())
# set data
load(file = "X2.RData")
XData <- X2

for (i in 1:100){
  timestamp()
  # prepare data
  X <- matrix(XData[i,], nrow = 1, dimnames = list(NULL, colnames(XData)))
  
  # model1
  modelevX2$mod1[[i]] <- bgphy(model = Model1, X = X)
  # model2
  modelevX2$mod2[[i]] <- bgphy(model = Model2, X = X)
  # model3
  modelevX2$mod3[[i]] <- bgphy(model = Model3, X = X)
  print(i)
}

save(modelevX2, file = "modelevX2.RData")
timestamp()





set.seed(1234)
# -------------- model 3 --------------
modelevX3 <- list(mod1 = list(), mod2 = list(), mod3 = list())
# set data
load(file = "X3.RData")
XData <- X3

for (i in 1:100){
  timestamp()
  # prepare data
  X <- matrix(XData[i,], nrow = 1, dimnames = list(NULL, colnames(XData)))
  
  # model1
  modelevX3$mod1[[i]] <- bgphy(model = Model1, X = X)
  # model2
  modelevX3$mod2[[i]] <- bgphy(model = Model2, X = X)
  # model3
  modelevX3$mod3[[i]] <- bgphy(model = Model3, X = X)
  print(i)
}

save(modelevX3, file = "modelevX3.RData")
timestamp()





# ========================= plot the results ==========================

# load the posterior inference results
load(file = "modelevX1.RData")
load(file = "modelevX2.RData")
load(file = "modelevX3.RData")

# retrieve posterior predictive loss for each data
loss1 <- sapply(modelevX1, function(mod){sapply(mod, function(x){x$loss})})
loss2 <- sapply(modelevX2, function(mod){sapply(mod, function(x){x$loss})})
loss3 <- sapply(modelevX3, function(mod){sapply(mod, function(x){x$loss})})



{
  pdf(file = "Fig6.pdf", height = 5, width = 12)
  layout(matrix(c(1:3, rep(4,3)), nrow = 2, ncol = 3, byrow = TRUE), heights = c(8,1))
  par(mar = c(1,5,3,1))
  par(cex.main = 2)
  par(cex.lab = 1.5)
  par(cex.axis = 1.5)
  boxplot(as.data.frame(loss1), xaxt = "n", boxwex = 0.5, main = expression(bold(x[tips]^("1"))), ylab = "loss")
  grid()
  boxplot(as.data.frame(loss1), xaxt = "n", add = TRUE, boxwex = 0.5,
          col = c("darkolivegreen1", "chartreuse3", "darkolivegreen"))
  boxplot(as.data.frame(loss2), xaxt = "n", boxwex = 0.5, main = expression(bold(x[tips]^("2"))), ylab = "loss")
  grid()
  boxplot(as.data.frame(loss2), xaxt = "n", add = TRUE, boxwex = 0.5,
          col = c("darkolivegreen1", "chartreuse3", "darkolivegreen"))
  
  boxplot(as.data.frame(loss3), xaxt = "n", boxwex = 0.5, main = expression(bold(x[tips]^("3"))), ylab = "loss")
  grid()
  boxplot(as.data.frame(loss3), xaxt = "n", add = TRUE, boxwex = 0.5,
          col = c("darkolivegreen1", "chartreuse3", "darkolivegreen"))
  
  par(mar = c(0,0,0,0))
  plot.new()
  legend("center", legend = c(expression(bold(M)[1]), expression(bold(M)[2]), expression(bold(M)[3])),
         col = c("darkolivegreen1", "chartreuse3", "darkolivegreen"), pch = c(15,15,15),
         bty = "n", pt.cex = 4, inset = -0.1, y.intersp = 1.5, x.intersp = 2, horiz = TRUE,
         text.width = c(0.1), cex = 2)
  dev.off()
}
