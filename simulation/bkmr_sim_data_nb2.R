
############################################
# NB BKMR: script for simulating data 
# 10/15/2022                                 
##########################################

library(lars); library(lasso2); library(mvtnorm); library(SuppDists); library(MCMCpack); library(splines2);
library(grplasso); library(magic); library(kernlab); library(MASS); library(fields); library(stats); library(msm); library(bkmr);library(Hmisc)
library(BayesLogit);library(tidyverse);  




# setwd("~/")

setwd("C:\\Users\\mutis\\Documents\\Dissertation\\paper 2\\R code paper 2\\nb bkmr - simulation")

# set.seed(21010)

# set.seed(010122)

# set.seed(021622)

# set.seed(1107)

# set.seed(3142023)

set.seed(021622)

########################
#  Spatial Data        #
########################


# A <- matrix(scan("~/ordered_adjmat_sc.txt"), 46, 46)      # SC adjacency matrix

A <- matrix(scan("ordered_adjmat_sc.txt"), 46, 46)

m<-apply(A,1,sum)	                                                      # No. of neighbors
ncounty<-n <- ncol(A)                                                   # No. of counties or spatial units
#Q<-as.spam(diag(m))-as.spam(A)+diag(0.0001,ncounty)                    # Ridge by 0.0001 to invert and generate data (negligible impact on results)
kappa<-.999999			  	                                                # Spatial dependence parameter ~ 1 for intrinsic CAR
true.s2phi <- s2phi <- 0.2
Q<-as.spam(diag(m))-kappa*as.spam(A)
covphi<-solve(Q)*s2phi			                                            # ncounty x ncounty covariance matrix of phis
phi<-c(rmvnorm(1, sigma=covphi))		                                    # ncounty vector of spatial effects

true.phi <- phi<-phi-mean(phi)                                          # ncounty x 1 phi1 vector -- Centered


##################
# Adjacency data #
##################

adjid <- rep(1:ncounty, m)                                              # Tells which elements of A belong to which block
adj <- rep(0, sum(m))                                                   # initialize a vector- size of the total # of neighbors
for (i in 1:ncounty) adj[adjid==i] <- which(A[i, ]==1)                  # returns the row indices for the neighbor units


#################
# Temporal Data	#
#################         
lday<-365                  # Number of time points (lday = "last day")
nis<-rep(lday,ncounty)     # Number of obs for each county
id<-rep(1:ncounty,nis)     # County id indicator for each day
N <-length(id)             # Total sample size (N in paper)
days<-rep(1:lday,ncounty)  # Repeat 1:lday for each ncounty 
pop<-rep(sample(10000:1000000,ncounty,replace=T),nis)
lpop<-log(pop)



######################################
# SC svi data                        #
######################################

svi_sc <- read_csv("SVI2018_US_COUNTY.csv")%>%filter(ST_ABBR%in%"SC")

svi_sc_sub <- svi_sc%>%select(EP_POV, EP_UNEMP, EP_PCI, EP_NOHSDP, EP_AGE65, EP_AGE17, EP_DISABL, EP_SNGPNT, EP_MINRTY, EP_LIMENG, EP_MUNIT, EP_MOBILE, EP_CROWD, EP_NOVEH, EP_GROUPQ)
M = 15                       # number of svi components
zi <- svi_sc_scaled <- apply(as.matrix(svi_sc_sub), 2, scale)      # n x M
Z <- svi_sc_scaled[rep(seq_along(nis), nis), ]                     # N x M


# M <- 5

# Sigma <-  matrix(data = c(1,0.2,0.3,0.45,0.15, 0.2,1,0.4,0.45,0.3, 0.3,0.4,1,0.3,0.15, 0.45,0.45,0.3,1,0.2, 0.15,0.3,0.15,0.2,1), nrow=M)

# Sigma <- diag(c(0.25, 0.3, 0.1, 0.4, 0.7))

# Z.data <- matrix(mvrnorm(n = ncounty, mu = c(rep(0, M)), Sigma = Sigma), byrow = FALSE, nrow = ncounty )

# Z.data <-  matrix(mvrnorm(n = ncounty, mu=c(rep(0,M)), Sigma = matrix(data = c(1,0.2,0.3,0.45,0.15, 0.2,1,0.4,0.45,0.3, 0.3,0.4,1,0.3,0.15, 0.45,0.45,0.3,1,0.2, 0.15,0.3,0.15,0.2,1), nrow=M)), byrow=FALSE, nrow=n)

# zi <- apply(Z.data, 2, scale)

# Z <- zi[rep(seq_along(nis), nis), ]





#####################################################
# generate the true exposure response function: h   #
#####################################################

# true.h1 = 0.25*(Z[,1]+Z[, 2] +Z[, 4]+Z[,9]+ 0.25*(Z[,5]^2-Z[,3]^2) - Z[,3]+Z[, 5]+Z[, 14]+0.25*Z[, 1]*Z[, 14])

# true.h1 = 0.25*(Z[,1]+Z[, 2] +Z[, 4]+Z[,9]+ 0.25*(Z[,5]^2-Z[,3]^2) - Z[,3]+Z[, 5]+Z[, 13]+0.25*Z[, 1]*Z[, 4])

# true.h1 = 0.25*(Z[,1] - Z[,3]+Z[, 5]+ Z[, 9]+Z[, 13]+ Z[,12]+Z[, 14]+0.25*(Z[,5]^2-Z[,3]^2) +0.25*Z[, 1]*Z[, 14])

# true.h1 <- 0.5*(Z[,1]^2 + Z[,2]^2 - Z[,4]^2 - Z[,5]^2 + 0.5*Z[,4]*Z[,5] + Z[,4] + Z[,5])           # Shelly's h function

h1 <- rep(0, ncounty)
true.h1 <- rep(h1, nis)         # null model



#############################
# Fixed-effect Spline Coefs #
# for Overall Time Trend    #
#############################

knots<-seq(14, 351, by=14)                                #  interior knots
Wtmp<-bSpline(1:lday,knots=knots,intercept = F)       
k<-ncol(Wtmp)                                           
W<-apply(Wtmp,2,rep,ncounty)


# spline coefficients

gamma <- 1*c(1.05, -0.26,  1.47,  0.06,  0.19, -0.43, -0.31, -1.36, -0.68, -0.90, -1.81, -1.48, -1.81, -2.27, -1.92, -1.68,
           -0.36,  0.78,  0.71,  1.25,  0.70,  0.59, -0.43, -0.54, -1.13, -0.33, -1.92, -0.16)               


beta <- c(-0.5, gamma)               # bind the fixed effect and spline coefficients into one vector

x1 <- rep(rnorm(ncounty, mean = 0, sd = 1), nis)
X <- cbind(x1, W)                    # design matrix for the fixed effect and splines

p <- ncol(X)

eta <- c(X%*%beta+true.h1+rep(phi, nis)+lpop-log(1000000))
psi <- exp(eta)/(1+exp(eta))
prob <- 1-psi
r <- 1

y <- Y <- rnbinom(N, size = r, prob = prob)



################################
# MLEs from Fixed Effect Model #
################################

nb<-glm.nb(Y~X-1+offset(lpop-log(1000000)))


#######################
# Histogram of counts #
#######################

range(Y)
tmp <- table(Y)/N*100
barplot(tmp, ylab = 'Percent', xlab = 'Count', col = "slateblue4")





