# model integrating capture-recapture data and pup counts 
# for elephant seals in R with a maximum likelihood approach
# october 2018

#------------------------------------------------------------------#
#---------- DEVIANCE OF THE INTEGRATED MODEL ----------------------#
#------------------------------------------------------------------#

dev_integrated_model <- function(b,dat,data,eff,e,garb,nh,km1){
	
#----------- 1. PARAMETERS
# b[1:14] = [phiPB0,phiPB1,phiPB2,phiPB3,phiPB4p,phi,psiPB_B2,psiPB_B3,psiPB_B4,psiB_B,psiNB_B,p1_PB,p2_B,p1_NB]
# b[15] = sigma observation error on counts

#----------- 2. DATA INPUTS
# dat is a data.frame containing N = number of pups, R2 = social structure variable used to calculate fertility function	
# data contains the encounter histories
# eff counts
# e vector of dates of first captures
# garb vector of initial states 
# km1 nb of recapture occasions (nb of capture occ - 1)
# nh nb ind

#----------- 3. CAPTURE-RECAPTURE LIKELIHOOD

# OBSERVATIONS (+1)
# 0 : not seen
# 1 : seen outside of breeding
# 2 : seen breeding
  
# STATES
# PB : pre-breeder
# B : breeder
# NB : non-breeder
# D : dead

# logit link for all parameters
lb <- 1/(1+exp(-b))
phiPB0 <- lb[1]
phiPB1 <- lb[2]
phiPB2 <- lb[3]
phiPB3 <- lb[4]
phiPB4p <- lb[5]
phi <- lb[6]
psiPB_B2 <- lb[7]
psiPB_B3 <- lb[8]
psiPB_B4 <- lb[9]
psiB_B <- lb[10]
psiNB_B <- lb[11]
p1_PB <- lb[12]
p2_B <- lb[13]
p1_NB <- lb[14]

#-- prob of obs (rows) cond on states (col)

# capture
B <- matrix(c(
1-p1_PB,p1_PB,0,
1-p2_B,0,p2_B,
1-p1_NB,p1_NB,0,
1,0,0),nrow=4,ncol=3,byrow=T)
B <- t(B)

#-- first encounter

BE <- matrix(c(
0,1,0,
0,0,1,
0,1,0,
1,0,0),nrow=4,ncol=3,byrow=T)
BE <- t(BE) 

#-- prob of states at t+1 given states at t

# survival
A1_age0 <- matrix(c(
phiPB0,0,0,1-phiPB0,
0,0,0,1,
0,0,0,1,
0,0,0,1),nrow=4,ncol=4,byrow=T)

A1_age1 <- matrix(c(
phiPB1,0,0,1-phiPB1,
0,0,0,1,
0,0,0,1,
0,0,0,1),nrow=4,ncol=4,byrow=T)

A1_age2 <- matrix(c(
phiPB2,0,0,1-phiPB2,
0,0,0,1,
0,0,0,1,
0,0,0,1),nrow=4,ncol=4,byrow=T)

A1_age3 <- matrix(c(
phiPB3,0,0,1-phiPB3,
0,phi,0,1-phi,
0,0,phi,1-phi,
0,0,0,1),nrow=4,ncol=4,byrow=T)

A1_age4p <- matrix(c(
phiPB4p,0,0,1-phiPB4p,
0,phi,0,1-phi,
0,0,phi,1-phi,
0,0,0,1),nrow=4,ncol=4,byrow=T)

# breeding
A2_age0 <- matrix(c(
1,0,0,0,
0,0,1,0,
0,0,1,0,
0,0,0,1),nrow=4,ncol=4,byrow=T)

A2_age1 <- matrix(c(
1,0,0,0,
0,0,1,0,
0,0,1,0,
0,0,0,1),nrow=4,ncol=4,byrow=T)

A2_age2 <- matrix(c(
1-psiPB_B2,psiPB_B2,0,0,
0,0,1,0,
0,0,1,0,
0,0,0,1),nrow=4,ncol=4,byrow=T)

A2_age3 <- matrix(c(
1-psiPB_B3,psiPB_B3,0,0,
0,psiB_B,1-psiB_B,0,
0,0,1,0,
0,0,0,1),nrow=4,ncol=4,byrow=T)

A2_age4p <- matrix(c(
1-psiPB_B4,psiPB_B4,0,0,
0,psiB_B,1-psiB_B,0,
0,psiNB_B,1-psiNB_B,0,
0,0,0,1),nrow=4,ncol=4,byrow=T)

A_age0 <- A1_age0 %*% A2_age0
A_age1 <- A1_age1 %*% A2_age1
A_age2 <- A1_age2 %*% A2_age2
A_age3 <- A1_age3 %*% A2_age3
A_age4p <- A1_age4p %*% A2_age4p

# init states
PI <- c(1,0,0,0)

# likelihood
   l <- 0
   for (i in 1:nh) # loop on ind
   {
      ei <- e[i] # date of first det
      oe <- garb[i] + 1 # init obs
      evennt <- data[,i] + 1 # add 1 to obs to avoid 0s in indexing
      ALPHA <- PI*BE[oe,]
      for (j in (ei+1):(km1+1)) # cond on first capture
      {
        if (j == (ei+1)) ALPHA <- (ALPHA %*% A_age0)*B[evennt[j],]
        if (j == (ei+2)) ALPHA <- (ALPHA %*% A_age1)*B[evennt[j],]
        if (j == (ei+3)) ALPHA <- (ALPHA %*% A_age2)*B[evennt[j],]
        if (j == (ei+4)) ALPHA <- (ALPHA %*% A_age3)*B[evennt[j],]
        if (j > (ei+4)) ALPHA <- (ALPHA %*% A_age4p)*B[evennt[j],]
     }
      l <- l + log(sum(ALPHA))*eff[i]
   }

dev_capturerecapture <- -2*l

#----------- 4. PUPS COUNTS LIKELIHOOD

# build the deviance function to minimize, see Eq 7 in Ferrari's paper
## we use data from 1973 to 1977 to initialize
# we predict pop size from 1978-2016

# get parameters
t_start <- 1978
sigma <- exp(b[15]) # observation error
pp <- phi
rho <- dat$rho[1]
qq <- 2 * rho * phiPB0 * (phiPB1 * phiPB2) # alpha is what Ferrari calls fertility constant 
																		# basically litter size, 1 in our case
dat$Npredict <- dat$N
Ft1 <- 1
Ft2 <- 1
dat$Npredict[dat$year==t_start] <- qq * Ft1 * (phiPB3 * (1-psiPB_B2) * dat$Npredict[dat$year==(t_start-4)] + psiPB_B2 * dat$Npredict[dat$year==(t_start-3)]) + pp * Ft1/Ft2 * dat$Npredict[dat$year==(t_start-1)]
for (t in (t_start+1):2016){
	Ft1 <- 1
	Ft2 <- 1
	dat$Npredict[dat$year==t] <- qq * Ft1 * (phiPB3 * (1-psiPB_B2) * dat$Npredict[dat$year==(t-4)] + psiPB_B2 * dat$Npredict[dat$year==(t-3)]) + pp * Ft1/Ft2 * dat$Npredict[dat$year==(t-1)]
}
dev_pupcounts <- sum((log(dat$N[dat$year>=t_start])-log(dat$Npredict[dat$year>=t_start]))^2)/(sigma*sigma) + 2 * log(sigma)

dev <- dev_capturerecapture + dev_pupcounts
dev

}

#------------------------------------------------------------------#
#--------------- MODEL FITTING / PARAMETER ESTIMATION -------------#
#------------------------------------------------------------------#

# read in capture-recapture data
data <- read.table('capturerecapturedata.txt')
head(data)
data <- R2ucare::group_data(data[,1:(ncol(data)-1)],rep(1,nrow(data)))
eff <- data[,ncol(data)]
data <- data[,1:(ncol(data)-1)]
# define various quantities
nh <- dim(data)[1]
k <- dim(data)[2]
km1 <- k-1
# compute the date of first capture fc, and state at initial capture init.state
fc <- NULL
init.state <- NULL
for (i in 1:nh){
temp <- 1:k
fc <- c(fc,min(which(data[i,]!=0)))
init.state <- c(init.state,data[i,fc[i]])
}
# transpose data
data <- t(data)

# read in pups counts data
nb_females_all <- read.table('females.txt',header=T)
# Npups
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
# birth sex ratio
bsr <- read.table('sexratio.txt',header=T)
bsr <- bsr / 2
# altogether
dat <- data.frame(year=1973:2016,N,R2,rho=mean(as.matrix(bsr)))

# init values
binit <- c(rep(0.5,14), 1.6)

# evaluate the integrated pop model deviance at the initial values, just to check
dev_integrated_model(binit,dat,data,eff,fc,init.state,nh,km1)

# fit model
deb=Sys.time()
tmpmin <- optim(binit,dev_integrated_model,NULL,hessian=TRUE,dat,data,eff,fc,init.state,nh,km1,method="BFGS",control=list(trace=1, REPORT=1, maxit = 10000))
fin=Sys.time()
fin-deb 

# get estimates and back-transform
b <- tmpmin$par

data.frame(param=c('phiPB0',
'phiPB1',
'phiPB2',
'phiPB3',
'phiPB4p',
'phi',
'psiPB_B2',
'psiPB_B3',
'psiPB_B4',
'psiB_B',
'psiNB_B',
'p1_PB',
'p2_B',
'p1_NB',
'sigma'), integrated = 
c(1/(1+exp(-b[1:14])),exp(b[15])))

#      param integrated
1    phiPB0  0.5982805
2    phiPB1  0.7645472
3    phiPB2  0.7804589
4    phiPB3  0.7909626
5   phiPB4p  0.7117037
6       phi  0.7609248
7  psiPB_B2  0.3076801
8  psiPB_B3  0.6752931
9  psiPB_B4  0.5642536
10   psiB_B  0.8180319
11  psiNB_B  0.7270515
12    p1_PB  0.7764135
13     p2_B  0.9092006
14    p1_NB  0.8721307
15    sigma 10.2934250

# get confidence intervals
SE_link <- sqrt(diag(solve(tmpmin$hessian)))
plogis(cbind(b[1:14]-1.96*SE_link[1:14],b[1:14]+1.96*SE_link[1:14]))
(sigma_CI <- exp(c(b[15]-1.96*SE_link[15],b[15]+1.96*SE_link[15]))) 


# deviance is 
tmpmin$value # 43170.68

# AIC is
tmpmin$value + 2*16 # 43202.68
