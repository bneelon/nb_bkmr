#####################################
# Prediction code for new exposures #
# Postmean estimation functions     #
#####################################



newh.postmean <- function(znew = znew, sel = sel) {
  
  if(is.null(dim(znew))) znew = matrix(znew, nrow=1)
  if(class(znew)  == "data.frame") znew <- data.matrix(znew)
  znew = as.matrix(znew)
  pnew = nrow(znew)
  
  
  thin <- sel%%10==0
  
  Tau.vec <- MCMC$Tau[sel][thin]
  rho.vec <- MCMC$Sigma.Kernel[sel][thin]
  Beta.mat <- MCMC$Beta[sel, ][thin, ]
  X <- MCMC$X
  z.mat <- MCMC$z.mat
  omega.mat <- MCMC$omega.mat
  phi.mat <- MCMC$phi.mat
  zi <- zi               # all the 15 svi variables -- scaled version
  id <- MCMC$id
  nis <- MCMC$nis
  
  cofactor <- 1e-7
  
  
  dim <- dim(z.mat)[1]
  postmean.pred <- postvar.pred <- vector('list', dim)               # store the posterior mean and variance for the predicted (znew)
  
  for (j in seq_len(dim)){
    
    z <- c(z.mat[j, ])
    phi <- c(phi.mat[j, ])
    beta <- c(Beta.mat[j, ])
    omega <- c(omega.mat[j, ])
    rho <- rho.vec[j]
    Tau <- Tau.vec[j]
   
    gauss.kernel <- rbfdot(sigma = rho)
    Sigma.h <- kernelMatrix(kernel = gauss.kernel, zi)*Tau                           # The K matrix

    Sigma.h.new = "[<-"(matrix(0, (q*n+q*pnew), (q*n+q*pnew)), 1:nrow(Sigma.h), 1:ncol(Sigma.h), value = Sigma.h) # adding rows/cols of zeroes to square matrix
    

    G <- kernelMatrix(kernel = gauss.kernel, rbind(zi, znew))*Tau
    
    
    G.1.11 = G[ 1:ncounty, 1:ncounty ]
    G.1.22 = G[ (ncounty+1):(ncounty+pnew), (ncounty+1):(ncounty+pnew) ]
    G.1.21 = G[ (ncounty+1):(ncounty+pnew), 1:ncounty ]
    G.1.12 = G[ 1:ncounty, (ncounty+1):(ncounty+pnew) ]
    

    Sigma11 <- G.1.11                           
    Sigma12 <- G.1.12
    Sigma21 <- G.1.21
    Sigma22 <- G.1.22

    Sigma.h.inv = solve(Sigma.h + cofactor*diag(nrow(Sigma.h))) 
    s.h <- solve(diag(tapply(omega, id, sum))+Sigma.h.inv)                             # posterior variance of the observed
    m.h <-  s.h%*%(tapply(omega*(z-X%*%beta-rep(phi, nis)), id, sum))                  # posterior mean of the observed
    
    
    Sigma11.inv = Sigma.h.inv 
    
    Mu.hnew.updated = Sigma21 %*% Sigma11.inv %*% m.h
    Sig11.mh <- Sigma11.inv %*% m.h
    Sig21.mh <- Sigma21 %*% m.h
    Sigma.hnew.updated = Sigma22 - Sigma21 %*% Sigma11.inv %*% Sigma12 + Sigma21 %*% Sigma11.inv %*% s.h %*% Sigma11.inv %*% Sigma12
    
    
    postmean.pred[[j]] <- drop(Mu.hnew.updated)                                      # posterior mean of the predicted
   
  }
  pred.postmean <- apply(matrix(unlist(postmean.pred), ncol = pnew, byrow = T), 2, mean)
  # pred.postvar <-  matrix(apply(matrix(unlist(postvar.pred), ncol = 4, byrow = T), 2, mean), 2, byrow = T)
  return(pred.postmean)
}


