# model integrating capture-recapture data and pup counts 
# for elephant seals in R with a maximum likelihood approach
# october 2018

#------------------------------------------------------------------#
#---------- DEVIANCE OF THE INTEGRATED MODEL ----------------------#
#------------------------------------------------------------------#

dev_integrated_model <- function(b,dat,data,eff,e,garb,nh,km1){
	
#----------- 1. PARAMETERS
# b[1:14] = [phiPB0,phiPB1,phiPB2,phiPB3,phiPB4p,phi,psiPB_B2,psiPB_B3,psiPB_B4,psiB_B,psiNB_B,p1_PB,p2_B,p1_NB]
# b[15] = a in F = (1+R2^a)^(1/a)
# b[16] = sigma observation error on counts

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
aa <- -exp(b[15]) # a in F(t)
sigma <- exp(b[16]) # observation error
pp <- phi
rho <- dat$rho[1]
qq <- 2 * rho * phiPB0 * (phiPB1 * phiPB2) # alpha is what Ferrari calls fertility constant 
																		# basically litter size, 1 in our case
dat$Npredict <- dat$N
Ft1 <- (1 + dat$R2[dat$year==(t_start)]^aa)^(1/aa)
Ft2 <- (1 + dat$R2[dat$year==(t_start-1)]^aa)^(1/aa)
dat$Npredict[dat$year==t_start] <- qq * Ft1 * (phiPB3 * (1-psiPB_B2) * dat$Npredict[dat$year==(t_start-4)] + psiPB_B2 * dat$Npredict[dat$year==(t_start-3)]) + pp * Ft1/Ft2 * dat$Npredict[dat$year==(t_start-1)]
for (t in (t_start+1):2016){
	Ft1 <- (1 + dat$R2[dat$year==(t)]^aa)^(1/aa)
	Ft2 <- (1 + dat$R2[dat$year==(t-1)]^aa)^(1/aa)
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
binit <- c(rep(0.5,14), -1, 1.5)

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
'a',
'sigma'), integrated = 
c(1/(1+exp(-b[1:14])),-exp(b[15]),exp(b[16])))

#      param integrated
1    phiPB0  0.5984056
2    phiPB1  0.7647213
3    phiPB2  0.7806313
4    phiPB3  0.7910896
5   phiPB4p  0.7113818
6       phi  0.7611661
7  psiPB_B2  0.3077905
8  psiPB_B3  0.6751706
9  psiPB_B4  0.5641837
10   psiB_B  0.8181329
11  psiNB_B  0.7270977
12    p1_PB  0.7764373
13     p2_B  0.9090268
14    p1_NB  0.8731530
15        a -0.9874015
16    sigma  1.6316564

# conf intervals
SE_link <- sqrt(diag(solve(tmpmin$hessian)))
plogis(cbind(b[1:14]-1.96*SE_link[1:14],b[1:14]+1.96*SE_link[1:14]))


# deviance is 
tmpmin$value # 43166.99

# AIC is
tmpmin$value + 2*16 # 43198.99

#------------------------------------------------------------------#
#-------------------- PUPS COUNTS PREDICTION ----------------------#
#------------------------------------------------------------------#

t_start <- 1978

# now predict pups abundance using estimated parameters
pp <- 1/(1+exp(-b[6]))
a <- -exp(b[15])
phiPB0 <- 1/(1+exp(-b[1]))
phiPB1 <- 1/(1+exp(-b[2]))
phiPB2 <- 1/(1+exp(-b[3]))
phiPB3 <- 1/(1+exp(-b[4]))
psiPB_B2 <- 1/(1+exp(-b[7]))

qq <- 2 * mean(as.matrix(bsr)) * phiPB0 * (phiPB1 * phiPB2)
qq

dat$Npredict <- dat$N
Ft1 <- (1 + dat$R2[dat$year==(t_start)]^a)^(1/a)
Ft2 <- (1 + dat$R2[dat$year==(t_start-1)]^a)^(1/a)
dat$Npredict[dat$year==t_start] <- qq * Ft1 * (phiPB3 * (1-psiPB_B2) * dat$Npredict[dat$year==(t_start-4)] + psiPB_B2 * dat$Npredict[dat$year==(t_start-3)]) + pp * Ft1/Ft2 * dat$Npredict[dat$year==(t_start-1)]
for (t in (t_start+1):2016){
	Ft1 <- (1 + dat$R2[dat$year==(t)]^a)^(1/a)
	Ft2 <- (1 + dat$R2[dat$year==(t-1)]^a)^(1/a)
	dat$Npredict[dat$year==t] <- qq * Ft1 * (phiPB3 * (1-psiPB_B2) * dat$Npredict[dat$year==(t-4)] + psiPB_B2 * dat$Npredict[dat$year==(t-3)]) + pp * Ft1/Ft2 * dat$Npredict[dat$year==(t-1)]
}
dat
Npredict_point <- dat$Npredict

# compare predictions to observed values
plot(1973:2016,dat$N,xlab='year',ylab='Number of pups',ylim=c(50,800))
points(1973:2016,dat$Npredict,col='red')
legend(1990,800,c('Counts','Predicted values'), pch = c(1,1), col=c('black','red'))


#------------------------------------------------------------------#
#-------------------- PUPS COUNTS PREDICTION ----------------------#
#---------------------- WITH CONF INTERVALS ----------------------#
#------------------------------------------------------------------#

# get confidence intervals
SE_link <- sqrt(diag(solve(tmpmin$hessian)))
p_CI <- plogis(c(b[6]-1.96*SE_link[6],b[6]+1.96*SE_link[6]))
(a_CI <- -exp(c(b[15]-1.96*SE_link[15],b[15]+1.96*SE_link[15]))) # -0.8259838 -1.1803642
phiPB0_CI <- plogis(c(b[1]-1.96*SE_link[1],b[1]+1.96*SE_link[1]))
phiPB1_CI <- plogis(c(b[2]-1.96*SE_link[2],b[2]+1.96*SE_link[2]))
phiPB2_CI <- plogis(c(b[3]-1.96*SE_link[3],b[3]+1.96*SE_link[3]))
phiPB3_CI <- plogis(c(b[4]-1.96*SE_link[4],b[4]+1.96*SE_link[4]))
psiPB_B2_CI <- plogis(c(b[7]-1.96*SE_link[7],b[7]+1.96*SE_link[7]))
(sigma_CI <- exp(c(b[16]-1.96*SE_link[16],b[16]+1.96*SE_link[16]))) # 0.6114746 4.3539058

# generate boostrap values in conf intervals (use of beta distributions would be more elegant)
nbboot <- 250
pseudo_p <- runif(nbboot,p_CI[1],p_CI[2])
pseudo_a <- runif(nbboot,a_CI[2],a_CI[1])
pseudo_phiPB0 <- runif(nbboot, phiPB0_CI[1], phiPB0_CI[2])
pseudo_phiPB1 <- runif(nbboot, phiPB1_CI[1], phiPB1_CI[2])
pseudo_phiPB2 <- runif(nbboot, phiPB2_CI[1], phiPB2_CI[2])
pseudo_phiPB3 <- runif(nbboot, phiPB3_CI[1], phiPB3_CI[2])
pseudo_psiPB_B2 <- runif(nbboot, psiPB_B2_CI[1], psiPB_B2_CI[2])

res <- NULL
for (i in 1:nbboot){
# predict pups abundance using estimated parameters
t_start <- 1978
qq <- 2 * mean(as.matrix(bsr)) * pseudo_phiPB0[i] * (pseudo_phiPB1[i] * pseudo_phiPB2[i])
dat$Npredict <- dat$N
Ft1 <- (1 + dat$R2[dat$year==(t_start)]^pseudo_a[i])^(1/pseudo_a[i])
Ft2 <- (1 + dat$R2[dat$year==(t_start-1)]^pseudo_a[i])^(1/pseudo_a[i])
dat$Npredict[dat$year==t_start] <- qq * Ft1 * (pseudo_phiPB3[i] * (1-pseudo_psiPB_B2[i]) * dat$Npredict[dat$year==(t_start-4)] + pseudo_psiPB_B2[i] * dat$Npredict[dat$year==(t_start-3)]) + pseudo_p[i] * Ft1/Ft2 * dat$Npredict[dat$year==(t_start-1)]
for (t in (t_start+1):2016){
	Ft1 <- (1 + dat$R2[dat$year==(t)]^pseudo_a[i])^(1/pseudo_a[i])
	Ft2 <- (1 + dat$R2[dat$year==(t-1)]^pseudo_a[i])^(1/pseudo_a[i])
dat$Npredict[dat$year==t] <- qq * Ft1 * (pseudo_phiPB3[i] * (1-pseudo_psiPB_B2[i]) * dat$Npredict[dat$year==(t-4)] + pseudo_psiPB_B2[i] * dat$Npredict[dat$year==(t-3)]) + pseudo_p[i] * Ft1/Ft2 * dat$Npredict[dat$year==(t-1)]
}
res <- cbind(res,dat$Npredict)
}


dim(res)
L <- apply(res,1,quantile,probs=2.5/100) # lower bound boostrap conf interval on Npredict
U <- apply(res,1,quantile,probs=97.5/100) # upper bound boostrap conf interval on Npredict

# compare predictions to observed values
png("fig.png",res=300,width=8,height=8,unit='in') 
plotrix::plotCI(1973:2016, Npredict_point, ui=U, li=L,col='red',ylim=c(50,500),xlab='year',ylab='Number of pups',cex.lab=1.3, cex.axis=1.3)
points(1973:2016,dat$N,xlab='year',ylab='Number of pups')
legend(2000,500,c('Counts','Predicted values'), pch = c(1,1), col=c('black','red'))
dev.off()
