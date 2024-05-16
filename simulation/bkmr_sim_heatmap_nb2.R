 ########################################################
# Script for creating the heatmaps for simulation study #
#########################################################

source("bkmr_sim_pred_nb2.R")

ngrid <- 50
qs <- c(0.25, 0.5, 0.75, 0.05, 0.95)
j = 2 # median levels 
min.plot.dist <- 2.5 ## color plot points within this distance from an observed data point

ylim = c(-5, 5)

######################################################
# Estimated and true exposure-response relationships #
######################################################

z1.1 <- seq(quantile(zi[,1], qs[4]), quantile(zi[,1], qs[5]), length=ngrid)
z1.2 <- seq(quantile(zi[,2], qs[4]), quantile(zi[,2], qs[5]), length=ngrid)

svi_sc_sub <- scale(as.matrix(svi_sc_sub))


cross.sec.1 = cbind(expand.grid(z1=z1.1, z2 = z1.2), z3=median(svi_sc_sub[, 3]), z4=median(svi_sc_sub[, 4]), z5=median(svi_sc_sub[, 5]), z6=median(svi_sc_sub[, 6]),
                    z7=median(svi_sc_sub[, 7]), z8=median(svi_sc_sub[, 8]), z9=median(svi_sc_sub[, 9]), z10=median(svi_sc_sub[, 10]), z11=median(svi_sc_sub[, 11]),
                    z12=median(svi_sc_sub[, 12]), z13=median(svi_sc_sub[, 13]), z14=median(svi_sc_sub[, 14]), z15=median(svi_sc_sub[, 15]))

grid = expand.grid(z1=z1.1, z2 = z1.2)
mindists = rep(NA, nrow(grid))
for(i in seq_along(mindists)) {
  pt <- as.numeric(grid[i,])
  dists <- as.matrix(dist(rbind(pt, zi[,1:2])))["pt",-1]
  mindists[i] <- min(dists)
}
rm(grid, pt, dists)

# Specific cross-sectional graphs: 

hgrid.T.1 <- newh.postmean(znew=cross.sec.1, sel=sel)  # this function returns a vector of posterior mean estimates
post.means <- hgrid.T.1
range(hgrid.T.1)

hgrid.T.1[mindists > min.plot.dist] <- NA
hgrid.T.1 = matrix(hgrid.T.1, nrow=ngrid)

# Truth:

# true.h1 = 0.25*(Z[,1]+Z[, 2] +Z[, 4]+Z[,9]+ 0.25*(Z[,5]^2-Z[,3]^2) - Z[,3]+Z[, 5]+Z[, 13]+0.25*Z[, 1]*Z[, 4])      

true.1 = matrix(0.25*(cross.sec.1[,1] + cross.sec.1[,2] + cross.sec.1[, 4]+cross.sec.1[, 9]+ 0.25*(cross.sec.1[,5]^2 - cross.sec.1[,3]^2) - 
                        cross.sec.1[,3]+cross.sec.1[, 5]+cross.sec.1[, 13]+0.25*cross.sec.1[, 1]*cross.sec.1[, 4] ), nrow = ngrid)

range(true.1)


######################################################################################
# Plotting a panel of the estimated NB BKMR heatmap and the true heatmap 
######################################################################################

set.panel()

zlim = c(-.25, .2)

zlim = c(-.61, 1.11)

par(oma=c( 0,0,0,3)) # margin of 3 spaces width at right hand side
set.panel( 1, 2)     # 1X2 matrix of plots

# draw plots using image command (or ggplot)

image(z1.1, z1.2, true.1, xlab="z1", ylab="z2", col=tim.colors(), main="True h1(z1, z2)", zlim=zlim)

image(z1.1, z1.2, hgrid.T.1, xlab="z1", ylab="z2", col=tim.colors(), main="Estimated h1(z1, z2)", zlim=zlim)

par(oma=c( 0,0,0,1))    # reset margin

image.plot( legend.only = TRUE, zlim=zlim) 

set.panel() # reset plotting device
