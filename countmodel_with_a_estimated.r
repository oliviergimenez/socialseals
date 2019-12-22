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
aa <- -exp(theta[2])
sigma <- exp(theta[3])
pp <- 0.76
pi <- 0.3

dat$Npredict <- dat$N

# first year
Ft1 <- (1 + dat$R2[dat$year==(t_start)]^aa)^(1/aa)
Ft2 <- (1 + dat$R2[dat$year==(t_start-1)]^aa)^(1/aa)
# RP: modified to pi starting at 3 and 1-pi starting at 4
rho <- dat$rho[1]
dat$Npredict[dat$year==t_start] <- qq * rho * Ft1 * (pp*(1-pi)*dat$Npredict[dat$year==(t_start-4)] + pi*dat$Npredict[dat$year==(t_start-3)]) + pp * Ft1/Ft2 * dat$Npredict[dat$year==(t_start-1)]
# then subsequent years
for (t in (t_start+1):2016){
	Ft1 <- (1 + dat$R2[dat$year==(t)]^aa)^(1/aa)
	Ft2 <- (1 + dat$R2[dat$year==(t-1)]^aa)^(1/aa)
	dat$Npredict[dat$year==t] <- qq * rho * Ft1 * (pp*(1-pi)*dat$Npredict[dat$year==(t-4)] + pi*dat$Npredict[dat$year==(t-3)]) + pp * Ft1/Ft2 * dat$Npredict[dat$year==(t-1)]
}
sum((log(dat$N[dat$year>=t_start])-log(dat$Npredict[dat$year>=t_start]))^2)/(sigma*sigma) + 2 * log(sigma)
}

# initial values
set.seed(1)
theta_init <- c(runif(1,0,1),runif(1,-5,0),runif(1,0,1)) # qq, a, sigma

# evaluate deviance
dev_seals(theta_init,dat)

# minimization
tmpmin <- optim(par = theta_init, fn = dev_seals,gr = NULL, dat = dat, hessian=TRUE, method="BFGS", control=list(trace=1, REPORT=1, maxit = 10000))

# get parameter estimates
(qq <- 1/(1+exp(-tmpmin$par[1]))) # 0.5333777
(aa <- -exp(tmpmin$par[2])) # -2.081313
(sigma <- exp(tmpmin$par[3])) # 1.620026

# standard errors on the link scale
SE_l <- sqrt(diag(solve(tmpmin$hessian)))

# confidence intervals
plogis(c(tmpmin$par[1] - 1.96 * SE_l[1],tmpmin$par[1] + 1.96 * SE_l[1])) # 0.01460742 0.98878173
-exp(c(tmpmin$par[2] - 1.96 * SE_l[2],tmpmin$par[2] + 1.96 * SE_l[2])) # -8.192536e-04 -5.287576e+03
exp(c(tmpmin$par[3] - 1.96 * SE_l[3],tmpmin$par[3] + 1.96 * SE_l[3])) # 0.6080048 4.3165502

# deviance is 
tmpmin$value # 1.964853

# AIC is 
tmpmin$value + 2*3 # 7.964853


# now predict pups abundance using estimated parameters
t_start <- 1978
dat$Npredict <- dat$N
pp <- 0.76
pi <- 0.3
rho <- dat$rho[1]
a <- -exp(tmpmin$par[2]) # -2.081313
Ft1 <- (1 + dat$R2[dat$year==(t_start)]^a)^(1/a)
Ft2 <- (1 + dat$R2[dat$year==(t_start-1)]^a)^(1/a)
dat$Npredict[dat$year==t_start] <- qq * rho * Ft1 * (pp*(1-pi)*dat$Npredict[dat$year==(t_start-4)] + pi*dat$Npredict[dat$year==(t_start-3)]) + pp * Ft1/Ft2 * dat$Npredict[dat$year==(t_start-1)]
# then subsequent years
for (t in (t_start+1):2016){
	Ft1 <- (1 + dat$R2[dat$year==(t)]^aa)^(1/aa)
	Ft2 <- (1 + dat$R2[dat$year==(t-1)]^aa)^(1/aa)
dat$Npredict[dat$year==t] <- qq * rho * Ft1 * (pp*(1-pi)*dat$Npredict[dat$year==(t-4)] + pi*dat$Npredict[dat$year==(t-3)]) + pp * Ft1/Ft2 * dat$Npredict[dat$year==(t-1)]
}
#dat

# compare predictions to observed values
plot(1973:2016,dat$N,xlab='year',ylab='Number of pups',ylim=c(50,800))
points(1973:2016,dat$Npredict,col='red')
legend(1990,800,c('Counts','Predicted values'), pch = c(1,1), col=c('black','red'))

