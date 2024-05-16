###################################
# Initialize          
###################################


sig.shape <- 5
sig.scale <- 10
sig.sq <- rinvgamma(1, shape = sig.shape, scale = sig.scale)


nu <- 5
delta <- 2
tau <- 1
q = 1                # random effect col dimension: (=1 here) since we add random intercept only


###########
# h
###########

rho <- 0.07
gauss <- rbfdot(sigma = rho)                 # Initialize the Gaussian radial basis function kernel


Sigma.h <- kernelMatrix(gauss, zi)*tau

if (doginv) {
  Sigma.h.inv		  = ginv(Sigma.h)
} else if (docofactor) {
  Sigma.h.inv		= solve(Sigma.h + diag(n) * cofactor)   # q = 1
} else Sigma.h.inv 	= solve(Sigma.h)
mean.h <- rep(0, n)

h <- mvrnorm(1, mu=mean.h, Sigma=Sigma.h)

h <- rep(0, ncounty)                  # null simulation

###################
# nb parameter r  #
###################

s <- 0.0018                           # tuning parameter for proposal variance of the nb parameter r  
r <- 1.2
accr <- 0                             # keep track of acceptance for the nb parameter r
accrho <- 0                           # keep track of acceptance for sigma (smoothing parameter for rbf kernel)

l <- rep(0, N)                        # for Gibbs update of r


##########################
# spatial effects        #
##########################


phi_init<-c(rmvnorm(1, sigma=diag(.01, ncounty)))	          # Random effects -- can just generate this from a univariate normal as well i.e. rnorm(ncounty, 0, 0.1)
phi<-phi_init-mean(phi_init)
s2phi <- var(phi)
sphi <- sqrt(s2phi)
tauphi<-1/s2phi

phibar <- rep(0, ncounty)                                   # mean of adjacent phi's for updating phi_i
for (j in 1:ncounty) phibar[j]<-mean(phi[adj[which(adjid==j)]])


#############
# priors
#############

beta0 <- beta <- rep(0, p)
beta.num <- length(beta)
T0 <- diag(0.01, p)

c<-d<-1		                                                # Gamma hyper priors for tauphi


#########################################
# FOR STORING POSTERIOR ESTIMATES       #
#########################################

sigma.kernel <- r.nb <- Sphi.sq <-  rep(0, num.reps)
h.post <- matrix(rep(0, num.reps*n), ncol = ncounty)
tau.post <- matrix(rep(0, num.reps), ncol = 1)
beta.post <- matrix(rep(0, num.reps*beta.num), nrow = num.reps)
h <- matrix(h, nrow = n)

burn <- num.reps/2
thinning = 5
lastit <- ((num.reps-burn)/(thinning))


z.mat <- omega.mat <- matrix(rep(0, lastit*length(id)), nrow = lastit)
phi.mat <- matrix(rep(0, lastit*ncounty), nrow = lastit)
Lwaic <- matrix(0, lastit, N)




