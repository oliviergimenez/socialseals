### Count model

# Read in Marion Island data
nb_females_all <- read.table('females.txt',header=T)
nb_pups <- read.table('pups.txt',header=T)
N <- nb_pups$number_pups
# Nb of females per males
nb_females_per_males <- nb_females_all$nb_females[nb_females_all$component=='females_per_adult_male'] 
# Adult sex ratio (males per female)
adult_sexratio <- 1 / nb_females_per_males 
# Harem size
harem_size <- nb_females_all$nb_females[nb_females_all$component=='harem_size'] 
# R2 variable
R2 <- harem_size * adult_sexratio
# alpha = psiB_B
alpha <- 0.81
# phiPB0
phi0 <- 0.6
# birth sex ratio
bsr <- read.table('sexratio.txt',header=T)
bsr <- bsr / 2
# dataset
dat <- data.frame(year = 1973:2016, N, R2, rho = mean(as.matrix(bsr)))
dat

# build the function S to minimize, see Eq 7
dev_seals <- function(theta,dat) {
## theta is the vector of parameter (q,a) with 
# - q = rho * r * alpha 
# - F = (1+R2^a)^(1/a)
## dat is a data.frame containing the data
# - N the number of pups
# - R2 the social structure variable that is used to calculate the fertility function	
# - bsr birth sex ratio

## we use data from 1973 to 1977 to initialize
# we predict pop size from 1978-2016

# get parameters
t_start <- 1978
qq <- 1/(1+exp(-theta[1]))
sigma <- exp(theta[2])
pp <- 0.76
pi <- 0.3

dat$Npredict <- dat$N

# first year
Ft1 <- 1
Ft2 <- 1
# RP: modified to pi starting at 3 and 1-pi starting at 4
rho <- dat$rho[1]
dat$Npredict[dat$year==t_start] <- qq * rho * Ft1 * (pp*(1-pi)*dat$Npredict[dat$year==(t_start-4)] + pi*dat$Npredict[dat$year==(t_start-3)]) + pp * Ft1/Ft2 * dat$Npredict[dat$year==(t_start-1)]
# then subsequent years
for (t in (t_start+1):2016){
	Ft1 <- 1
	Ft2 <- 1
	dat$Npredict[dat$year==t] <- qq * rho * Ft1 * (pp*(1-pi)*dat$Npredict[dat$year==(t-4)] + pi*dat$Npredict[dat$year==(t-3)]) + pp * Ft1/Ft2 * dat$Npredict[dat$year==(t-1)]
}
sum((log(dat$N[dat$year>=t_start])-log(dat$Npredict[dat$year>=t_start]))^2)/(sigma*sigma) + 2 * log(sigma)
}

# initial values
set.seed(1)
theta_init <- c(runif(1,0,1),runif(1,0,1)) # qq, a, sigma

# evaluate deviance
dev_seals(theta_init,dat)

# minimization
tmpmin <- optim(par = theta_init, fn = dev_seals,gr = NULL, dat = dat, hessian=TRUE, method="BFGS", control=list(trace=1, REPORT=1, maxit = 10000))

# get parameter estimates
(qq <- 1/(1+exp(-tmpmin$par[1]))) # 0.4647273
(sigma <- exp(tmpmin$par[2])) # 1.680743

# standard errors on the link scale
SE_l <- sqrt(diag(solve(tmpmin$hessian)))

# confidence intervals
plogis(c(tmpmin$par[1] - 1.96 * SE_l[1],tmpmin$par[1] + 1.96 * SE_l[1])) # 0.4125385 0.5177002
exp(c(tmpmin$par[2] - 1.96 * SE_l[2],tmpmin$par[2] + 1.96 * SE_l[2])) # 0.6308016 4.4782638

# deviance is 
tmpmin$value # 2.038471

# AIC is 
tmpmin$value + 2*3 # 8.038471

