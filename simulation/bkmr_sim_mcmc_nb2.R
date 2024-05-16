
#######################################
# NB BKMR: MCMCM simulation code
#######################################


for (iter in 1:num.reps)  {
 
  # updating the nb parameter using  Gibbs step
  
  # Update latent counts, l, using Chinese restaurant table distribution (c.f. Zhou and Carin, Dadaneh)
  for(j in 1:N) l[j]<-sum(rbinom(y[j],1,round(r/(r+1:y[j]-1),6))) # Could try to avoid loop; in rounding avoids numerical stability

  # Update r from conjugate gamma distribution given l and psi
  eta <- X%*%beta+rep(h, nis)+rep(phi, nis)+lpop-log(1000000)
  psi<-exp(eta)/(1+exp(eta))
  r<-rgamma(1,.01+sum(l),0.01-sum(log(1-psi)))

  
  
  
  ####################################
  # update z - the latent response   #
  ####################################
 
  eta <- X%*%beta+rep(h, nis)+rep(phi, nis)+lpop-log(1000000)
  omega <- rpg(N, Y+r, eta)
  z <- (Y-r)/(2*omega)-lpop+log(1000000)
  
  
  ################
  # update beta  #
  ################
  
  vb <- solve(crossprod(X*sqrt(omega)))
  mb <- vb%*%(t(sqrt(omega)*X)%*%(sqrt(omega)*(z-rep(h, nis)-rep(phi, nis))))
  beta <-  mvrnorm(1, mu = mb, Sigma = vb)
  beta.post[iter,] <- beta 
  beta <-  matrix(beta, nrow=p)
  

  
  ###############
  # update phi  #
  ###############
  

  priorprec<-1/(sphi^2)*Q                                         # Prior precision of phi1
  priormean<- 0                                                   # Prior mean of phi1 is 0 
  prec<-priorprec+as.spam(diag(tapply(omega,id,sum),ncounty, ncounty))
  mb<-c(tapply(omega*(z-X%*%beta-rep(h, nis)),id,sum))                      # note priormean is 0 and only data contributes to mb
  phi<-rmvnorm.canonical(1, mb, prec)[1,]

  # center phi and update tauphi

  phi <- phi-mean(phi)

  tauphi<-rgamma(1, c+(ncounty-1)/2, d+t(phi)%*%Q%*%phi/2)                  # n-1 since rank of Q is n-1
  sphi <- sqrt(1/tauphi)

  
  Sphi.sq[iter] <- sphi^2                                                   # spatial random effect variance


 
  ###########################################
  # update tau - the covariance parameter   #
  ###########################################
  
  if (doginv) {
      term.gam		= t(h) %*% ginv(kernelMatrix(gauss, zi)) %*% h
    } else if (docofactor) {
      term.gam= t(h) %*% solve(kernelMatrix(gauss, zi) + diag(n)*cofactor) %*% h
    } else term.gam = t(h) %*% solve(kernelMatrix(gauss, zi)) %*% h
    
  tau <- rinvgamma(n = 1, shape = n/2+5, scale = (term.gam+2*10)/2)              # assuming IG (5, 10) prior for tau square - since we have one component--> no lambda.sq parameter

  tau.post[iter,] <-  tau
  
 
  ###############################
  # update rho for rbf kernel   #
  ###############################
  
  # Updating this parameter implies updating Sigma.h! 
  
  # proposal value
  
  rho.var <- 0.0003
  rhos <- rtnorm(1, mean = rho, sd = sqrt(rho.var), lower = 0)
  gauss.p <- rbfdot(sigma = rhos)
  Sigma.hs <- kernelMatrix(gauss.p, zi)*tau


  if (doginv) {
    Sigma.hs.inv		= ginv(Sigma.hs)
  } else if (docofactor) {
    Sigma.hs.inv		= solve(Sigma.hs + diag(n) * cofactor)
  } else Sigma.hs.inv 	= solve(Sigma.hs)

  cov.hs <-  solve(diag(tapply(omega, id, sum))+Sigma.h.inv)
  
  ratio<-sum(dmvnorm(c(h), mean = rep(0, n), sigma = Sigma.hs, log = T))-sum(dmvnorm(c(h), mean = rep(0, n), sigma = Sigma.h, log = T))+  
    log(dinvgamma(rhos, shape = 50, scale = 4))-log(dinvgamma(rho, shape = 50, scale = 4))+                              # assuming an IG prior for rho
    dtnorm(rho, rhos, sd = sqrt(rho.var),  lower = 0, log = T)-dtnorm(rhos, rho, sd = sqrt(rho.var),  lower = 0, log = T)


  if (log(runif(1))<ratio) {
    rho<-rhos
    accrho<-accrho+1
    gauss <- rbfdot(sigma = rho)
  }

  sigma.kernel[iter] <- rho               # store
  
  
  
  
  ##################################
  # update h                      #
  ##################################
  
  Sigma.h <- kernelMatrix(gauss, zi)*tau

  if (doginv) {
    Sigma.h.inv		= ginv(Sigma.h)
  } else if (docofactor) {
    Sigma.h.inv		= solve(Sigma.h + diag(n) * cofactor)
  } else Sigma.h.inv 	= solve(Sigma.h)                                                  # need to invert the Sigma.h in the posterior covariance of h


  cov.h <- solve(diag(tapply(omega, id, sum))+Sigma.h.inv)                              # posterior covariance of h
  mean.h <- cov.h%*%(tapply(omega*(z-X%*%beta-rep(phi, nis)), id, sum))                 # posterior mean of h
  h	<- matrix(mvrnorm(1, mu=mean.h, Sigma=cov.h), nrow = n)
  h.post[iter,] <-  c(h)
 
  
  
 
  
  ###################################################################
  # update r - the negative binomial parameter using a MH step      #
  ###################################################################
   
  # eta <- X%*%beta+rep(h, nis)+rep(phi, nis)+lpop-log(1000000)
  # psi <- exp(eta)/(1+exp(eta))
  # prob <- 1-psi
  # 
  # rnew<-rtnorm(1, r, sqrt(s),lower=0)       # Treat r as continuous
  # ratio<-sum(dnbinom(Y, rnew, prob,log=T))-sum(dnbinom(Y, r, prob,log=T))+
  #   dtnorm(rnew, r, sqrt(s),0,log=T)-dtnorm(r, rnew, sqrt(s), 0, log=T)  # Diffuse (e.g., uniform) prior for r
  # if (log(runif(1))<ratio) {
  #   r<-rnew
  #   accr<-accr+1
  # }

  r.nb[iter] <- r
  
  
  #####################################
  # store                             #
  #####################################

  if(iter>burn & iter%%thinning==0) {
    i<-(iter-burn)/thinning
    z.mat[i, ] <- z
    omega.mat[i, ] <- omega
    phi.mat[i, ] <- phi
    Lwaic[i, ] <- dnbinom(Y, size = r, prob = prob)
  }
  
  
  #########################
  # print betas           #
  #########################
  
  if(iter%%10==0){
    print(iter)
    print(c(beta))
    print(c(phi[1:15]))
    print(1/tauphi)
  }
  
}



