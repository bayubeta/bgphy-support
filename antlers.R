library(bgphy)
library(ape)

# code for application to the antlers dataset,
# for the paper Bayesian inference of mixed Gaussian phylogenetic models (2025)
# by Brahmantio, Bartoszek, Yapar.

# read data
dat <- readxl::read_xlsx("antler_data.xlsx")
# read tree
tree <- read.nexus(file = "cervidae.tre")
# =========== the required files above are available upon request from
# =========== doi.org/10.1007/s11692-023-09624-1


# active taxa
active <- intersect(dat$Taxa, tree$tip.label)

# active data
dat2 <- dat
rownames(dat2) <- dat$Taxa
dat2 <- dat2[active,]

# new tree from active taxa
antler_tree <- drop.tip(tree, setdiff(tree$tip.label, active))

# active taxa id
id <- match(active, dat$Taxa)

# prepare the continuous measurements
# posterior skull length
logPSL <- matrix(log(dat$mean_psl[id]), nrow = 1)
colnames(logPSL) <- active

# ------------------------------------------------------------------------------
# prepare prior hyperparameters
# for X0 and theta
x1 <- 1
x2 <- 50
q1 <- 0.1
q2 <- 0.9
sigma_normal <- (log(x2) - log(x1))/(qnorm(q2) - qnorm(q1))
mu_normal <- log(x1) - qnorm(q1)*sigma_normal

# for alpha
# tree height
tH <- max(phytools::nodeHeights(antler_tree))
# high alpha value ~ half-life = 1% of tH => 95% percentile
alpha_high <- log(2)/(0.01*tH)
p <- 0.95
s_alpha <- alpha_high/qnorm((p+1)/2)

# for sigma
sigma_high <- sigma_normal*sqrt(2/tH)
p <- 0.95
s_sigma <- sigma_high/qnorm((p+1)/2)

# ================== Regime config 1: Neutral OU regime ========================

OU <- setModel(tree = antler_tree, regime_names = "R1", modeltypes = "OU")

# set priors
# for X0 and theta \in R
OU$priors$X0 <- prior_normal(mean = mu_normal, sd = sigma_normal)
OU$priors$theta <- prior_normal(mean = mu_normal, sd = sigma_normal)

# for alpha and sigma in R+
OU$priors$alpha <- prior_halfnormal(sigma = s_alpha)
OU$priors$sigma <- prior_halfnormal(sigma = s_sigma)

# run posterior inference
set.seed(1234)
timestamp()
post_OU <- bgphy(model = OU, X = logPSL, nsample = 100000)
save(post_OU, file = "post_OU.RData")

# run posterior predictive simulation
post_pred_OU <- post_simulate(post_OU, nsample = 100000)



# ============= Regime config 2: Old-world vs New-world regimes ================

OldNew <- setModel(tree = antler_tree,
                   regime_names = c("Old_world", "New_world"),
                   modeltypes = c("OU", "OU"),
                   startNodes = list(Old_world = c(48), New_world = c(49)))

# set priors
OldNew$priors$X0 <- OldNew$priors$theta_1 <- OldNew$priors$theta_2 <- 
  prior_normal(mean = mu_normal, sd = sigma_normal)
OldNew$priors$alpha_1 <- OldNew$priors$alpha_2 <- 
  prior_halfnormal(sigma = s_alpha)
OldNew$priors$sigma_1 <- OldNew$priors$sigma_2 <- 
  prior_halfnormal(sigma = s_sigma)

# run posterior inference
set.seed(1234)
timestamp()
post_OldNew <- bgphy(model = OldNew, X = logPSL, nsample = 100000)
save(post_OldNew, file = "post_OldNew.RData")

# run posterior predictive simulation
post_pred_OldNew <- post_simulate(post_OldNew, nsample = 100000)



# ============= Regime config 3: regimes based on antler's shape =====================

Shape <- setModel(tree = antler_tree,
                  regime_names = c("MB", "Palm", "Bifurcate"),
                  modeltypes = c("OU", "OU", "OU"),
                  startNodes = list(MB = c("48", dat2$Taxa[dat2$Type == "MB"]),
                                    Palm = c("87", "51", dat2$Taxa[dat2$Type == "Palm"]),
                                    Bifurcate = c("75", dat2$Taxa[dat2$Type == "Bifurcate"])))

# set priors
Shape$priors$X0 <- Shape$priors$theta_1 <- Shape$priors$theta_2 <- Shape$priors$theta_3 <- 
  prior_normal(mean = mu_normal, sd = sigma_normal)
Shape$priors$alpha_1 <- Shape$priors$alpha_2 <- Shape$priors$alpha_3 <- 
  prior_halfnormal(sigma = s_alpha)
Shape$priors$sigma_1 <- Shape$priors$sigma_2 <- Shape$priors$sigma_3 <- 
  prior_halfnormal(sigma = s_sigma)

# run posterior inference
set.seed(1234)
timestamp()
post_Shape <- bgphy(model = Shape, X = logPSL, nsample = 100000)
save(post_Shape, file = "post_Shape.RData")

# run posterior predictive simulation
post_pred_Shape <- post_simulate(post_Shape, nsample = 100000)



