###############################
# NB BKMR Simulation Run Code 
###############################

setwd("~/")

rm(list = ls())

doMCMC 			= TRUE
doginv			= FALSE
docofactor	= TRUE
cofactor		= 1e-7	# randomly generated kernel matrix is not always invertible.

if (doMCMC)
{
  # note: the order in which scripts are loaded matters
  source("bkmr_sim_data_nb2.R")
  num.reps 	= 10000
  sel 		  = (num.reps/2+1):num.reps
  source("bkmr_sim_initialize_nb2.R")
  source("bkmr_sim_mcmc_nb2.R")
  MCMC = list(H = h.post, Tau = tau.post, Beta = beta.post, Z = Z, y = y, X = X, W = diag(ncounty), sel = sel, cofactor = cofactor, q = q, M = M, n = n, R = r.nb,
              Sigma.Kernel = sigma.kernel, S2PHI = Sphi.sq, Accrho = accrho, Accr = accr, id = id, nis = nis, zi = zi, z.mat = z.mat, omega.mat = omega.mat, phi.mat = phi.mat, Lwaic = Lwaic)  

}



##########################################################################
# fit lm() for estimated h on true.h1 to check how well h was estimated  #
##########################################################################

true.h1 <- 0.25*(zi[,1] + zi[,14] + zi[, 4]+zi[, 9]+0.25*(zi[,5]^2-zi[,3]^2) - zi[,3]+zi[, 5]+zi[,13]+0.25*zi[, 1]*zi[, 4]) # note:  zi is n x M

true.h1 = 0.25*(Z[,1]+Z[, 2] +Z[, 4]+Z[,9]+ 0.25*(Z[,5]^2-Z[,3]^2) - Z[,3]+Z[, 5]+Z[, 14]+0.25*Z[, 1]*Z[, 4])

true.h1 <- 0.25*(zi[,1] + zi[,14] + zi[, 4]+zi[, 9]+0.25*(zi[,5]^2-zi[,3]^2) - zi[,3]+ zi[, 5]+zi[, 14]+0.25*zi[, 1]*zi[, 4]) # note:  zi is n x M

true.h1 = 0.25*(zi[,1]+zi[, 2] +zi[, 4]+zi[,9]+ 0.25*(zi[,5]^2-zi[,3]^2)-zi[,3]+zi[, 5]+zi[, 13]+0.25*zi[, 1]*zi[, 4])

true.h1 = 0.25*(zi[,1]+zi[, 2] +zi[, 4]+zi[,9]+ 0.25*(zi[,5]^2-zi[,12]^2)-zi[,12]+zi[, 5]+zi[, 13]+0.25*zi[, 1]*zi[, 4])

true.h1 = 0.25*(zi[,1] - zi[,3]+zi[, 5]+ zi[, 9]+zi[, 13]+zi[,12]+ zi[,14]+0.25*(zi[,5]^2-zi[,3]^2) +0.25*zi[, 1]*zi[, 14])


true.h1 <- 0.5*(zi[,1]^2 + zi[,2]^2 - zi[,4]^2 - zi[,5]^2 + 0.5*zi[,4]*zi[,5] + zi[,4] + zi[,5])           # Shelly's h function


# true.h1 <- true.h1-mean(true.h1)
# sel <- 15001:30000
hhat.1 <- apply(MCMC$H[sel, ], 2, mean)
model.1 <- lm(hhat.1~true.h1)
summary(model.1)


###################
# thinned samples
##################

thin <- sel%%5==0

Beta <- MCMC$Beta[sel, ][thin, ]
Sigsq <- MCMC$Sigsq[sel][thin]
Sigma.Kernel <- MCMC$Sigma.Kernel[sel][thin]
Tau <- MCMC$Tau[sel][thin]
R <- MCMC$R[sel][thin]
H <- MCMC$H[sel, ][thin, ]
S2PHI <- MCMC$S2PHI[sel][thin]


Tau.vec <- MCMC$Tau
rho.vec <- MCMC$Sigma.Kernel
Beta.mat <- MCMC$Beta
X <- MCMC$X

Sigsq <- MCMC$Sigsq
Sigma.Kernel <- MCMC$Sigma.Kernel
R <- MCMC$R
H <- MCMC$H
S2PHI <- MCMC$S2PHI

z.mat <- MCMC$z.mat
omega.mat <- MCMC$omega.mat
phi.mat <- MCMC$phi.mat


save.image(file = "null_simulation_rev.RData")

#########################
# posterior summaries   #
#########################


colMeans(Beta)

apply(Beta, 2, quantile, probs = c(0.025, 0.975))

mean(Sigma.Kernel)
quantile(Sigma.Kernel, probs = c(0.025, 0.975))

mean(Tau)
quantile(Tau, probs = c(0.025, 0.975))

mean(S2PHI)
quantile(S2PHI, probs = c(0.025, 0.975))

mean(R)
quantile(R, probs = c(0.025, 0.975))


##############################
# Trace plots
##############################

par(mfrow = c(2, 2))

plot(Beta[, 1], type = 'l', col = "lightgreen", xlab = 'iteration', ylab = expression(beta[11]))
abline(h = mean(Beta[, 1]), col = 'darkblue')
plot(Beta[, 2], type = 'l', col = "lightgreen", xlab = 'iteration', ylab = expression(beta[12]))
abline(h = mean(Beta[, 2]), col = 'darkblue')
plot(Sigma.Kernel, type = 'l', col = "lightgreen", xlab = 'iteration', ylab = expression(rho))
abline(h = mean(Sigma.Kernel), col = 'darkblue')

plot(Tau, type = 'l', col = "lightgreen", xlab = 'iteration', ylab = expression(tau**2))
abline(h = mean(Tau), col = 'darkblue')

plot(Lambda, type = 'l', col = "lightgreen", xlab = 'iteration', ylab = expression(lambda))
abline(h = mean(Lambda), col = 'darkblue')

plot(R, type = 'l', col = "lightgreen", xlab = 'iteration', ylab = 'r')
abline(h = mean(R), col = 'darkblue')

plot(S2PHI, type = 'l', col = "lightgreen", xlab = 'iteration', ylab = expression((sigma[phi])**2))
abline(h = mean(S2PHI), col = 'darkblue')


dev.off()


############
# WAIC     #
############

lhat<-apply(Lwaic, 2, mean)
lpd<-sum(log(lhat))
pwaic<-sum(apply(log(Lwaic), 2, var))
waic<- -2*(lpd-pwaic)         
waic



# Trace for all betas

# par(mfrow = c(5, 6))
# 
# for ( i in 1: 29){
#   plot(Beta[, i], type = 'l', col = "lightgreen", xlab = 'iteration', ylab = expression(beta[i]))
#   abline(h = mean(Beta[, i]), col = 'darkblue')
# }
# 
# dev.off()



###########################################
# convergence diagnostics: Geweke p-vals  #
###########################################

zlist <- apply(Beta, 2, geweke.diag)
len <- length(zlist)
p.val <- vector(mode = 'double',length = len)
for (i in 1:len) {
  z <- zlist[[i]]$z
  p.val[i] <- 2*(ifelse(z>0, 1-pnorm(z), pnorm(z)))
}

sum(p.val>=0.05)


####################################
# plotting the splines             #
####################################


true.spline <- Wtmp%*%gamma

cls <- t(Wtmp%*%t(Beta[, 2:29]))

lcl <- apply(cls, 2, quantile, prob = 0.025)
ucl <- apply(cls, 2, quantile, prob = 0.975)

est.spline <- Wtmp%*%c(colMeans(Beta[, 2:29]))


plot(1:lday, true.spline, type = 'l', col = 'slateblue4', xlab = 'day', ylab = '', ylim = c(-2.5, 1.5))
lines(1:lday, est.spline, type = 'l', lty = 2, col = 'magenta')
lines(1:lday, lcl, type = 'l', lty = 2, col = 'blue')
lines(1:lday, ucl, type = 'l', lty = 2, col = 'blue')



splines_data <- tibble(day = c(1:lday), ft.true = c(true.spline), ft.est = c(est.spline), lcl = lcl, ucl = ucl)

splines <- ggplot(splines_data, aes(x = day, y = ft.true ))+
          geom_line(aes(color = "True Curve"), size = 0.5, linetype = 1)+
          geom_line(aes(day, ft.est, color="Estimated Curve"),size=0.5,linetype=2)+
          geom_ribbon(aes(ymin = lcl, ymax = ucl, color="95% Credible Interval"), alpha=0.1, show.legend = F)+
          scale_color_manual(breaks = c(
            "True Curve", 
            "Estimated Curve",
            "95% Credible Interval"
          ), 
          values = c(
            "True Curve"="blue4",
            "Estimated Curve" = "darkorange",
            "95% Credible Interval"="gray40"
          ))+
          ylab(expression(paste(f(t))))+
          ylim(-2.5, 1.4)+
          coord_fixed(ratio = 150/(2.5+1.4))+
          scale_x_continuous(name="Day",breaks=c(1, 32, 60, 91, 121, 152, 182, 213, 244, 274, 305, 335, 365)
          )+
          guides( color = guide_legend(title="Temporal Trends",
                                       override.aes = list(
                                         linetype = c(1, 2, 1),
                                         shape = c(NA, NA, NA)),
                                       reverse = F))+
          # ggtitle("(b)")+
          theme(legend.key = element_rect(fill = "white"),
                legend.position = c(0.35, 0.8),
                legend.text = element_text(size = 8),
                legend.title = element_text(size = 8, face = 'bold'),
                axis.text = element_text(size = 8, face = "bold"),
                axis.title = element_text(size = 8),
                plot.title = element_text(size = 8, hjust = 0.5, face = "bold"),
                panel.background = element_rect(colour = "grey50"),
                plot.margin = unit(c(0,0,0,0), "cm"))



ggsave(
  "~/splines_sim.pdf",
  device = pdf(),
  scale = 1,
  dpi = "print"(300),
  width = 7,
  height = 3.5,
  units = "in"
)

dev.off()



##################################
# random effect maps             #
##################################


##############################################################################################
# Uses shape files from TIGER/Line shapefiles of the US Census available from tigris package #
##############################################################################################

library(tigris)
library(rgdal)
library(sf)
library(tmap)
library(RColorBrewer)

sc_counties_map <- counties(state = c('45'))

GEOID <- sc_counties_map$GEOID



###################################
#  h and phi maps  
###################################



H <- MCMC$H
hhat <- apply(H, 2, mean)
hphi.data <- tibble(NAME = covid_data_sc$County[!duplicated(id)], h.true = true.h1, h.post = hhat, phi.true = true.phi, phi.post = c(apply(MCMC$phi.mat, 2, mean)))


maps.data <- left_join(sc_counties_map, hphi.data, by = 'NAME')



#######################################
# True vs Predicted: Count part       #
#######################################

pal <- brewer.pal(5,"BuGn")
pal



###############################################
# maps for posterior estimates for h and phi  #
###############################################



tm_pmh <- tm_shape(maps.data)+
  tm_fill(c("h.post"), midpoint = c(NA), title = c('h'), palette = pal, style = "quantile")+
  #tm_style("col_blind")+
  tm_facets(free.scales = TRUE, nrow = 1)+
  tm_layout(title = "Quintile Map",
            title.snap.to.legend = TRUE,
            title.size = 0.8,
            title.position = c("right", "bottom"),
            legend.outside = FALSE,
            legend.position = c("right", "top"),
            # legend.position = c(0.8, 0.85),
            main.title = "(a)",
            main.title.fontface = "bold",
            main.title.position = "center")+
  tm_borders(alpha = 0.3, lwd = 1)


tmap_save(tm_pmh, filename = "~/tm_pmh_app2.jpg", dpi = 300 )



tm_pmphi <- tm_shape(maps.data)+
  tm_fill(c("phi.post"), midpoint = c(NA), title = c(expression(paste(phi))), palette = pal, style = "quantile")+
  #tm_style("col_blind")+
  tm_facets(free.scales = TRUE, nrow = 1)+
  # tm_facets('hcut', nrow = 1)+
  tm_layout(title = "Quintile Map",
            title.snap.to.legend = TRUE,
            title.size = 0.8,
            title.position = c("right", "bottom"),
            legend.outside = FALSE,
            legend.position = c("right", "top"),
            # legend.position = c(0.8, 0.85),
            main.title = "(b)",
            main.title.fontface = "bold",
            main.title.position = "center")+
  tm_borders(alpha = 0.3, lwd = 1)


tmap_save(tm_pmphi, filename = "~/tm_pmphi_app2.jpg", dpi = 300 )












