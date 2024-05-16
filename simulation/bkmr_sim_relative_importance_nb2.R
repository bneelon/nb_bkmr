############################################################
# NB BKMR: Relative importance simulation code             #
############################################################


source("bkmr_sim_pred_nb2_rr.R")

mat <-  matrix(0, nrow=q, ncol=M)                           # stores the posterior mean
true.mat <- mat.ci <-  matrix(0, nrow = 2, ncol = M)        # stores the  95% Credible intervals for the posterior mean rr. Note: can just store both mean and CrIs in one matrix.

# Zs <- svi_sc_scaled
Zs <- zi

qs = c(0.25, 0.75, 0.50)

for (j in 1:M) {
  cross.sec 		= rbind(apply(Zs, 2, median), apply(Zs, 2, median))
  cross.sec[,j] 	= c(quantile(Zs[,j], qs[2]), quantile(Zs[,j], qs[1]))
  hgrid.cross.sec = newh.postmean_rr(znew = cross.sec, sel = sel)    
  rr_vector <- hgrid.cross.sec[, 1]/hgrid.cross.sec[, 2]             # a vector of risk ratios across all iterations
  mat[1,j] 	= mean(rr_vector)                                        # posterior mean rr estimate
  mat.ci[, j] <- c(unname(quantile(rr_vector, probs = c(0.025, 0.975))))
}

# plot

ylim <- c(0.5, 2)
plot(1:M, mat[1,], xaxt="n",
     ylim=ylim,
     pch=15, xlab="SVI", ylab="Main effect of each SVI variable",
)

arrows(1:M, mat.ci[1, ], 1:M, mat.ci[2, ], length=0.05, angle=90, code=3)
axis(1, at=1:M, labels=c("pov", "unemp", "pci", "nhsd", "age65", "age17", "disabl", "sngpnt", "mnrty", "limeng", "munit", "mobile", "crowd", "noveh", "grpq"))
abline(h=1)


#############################
# relative importance plots
#############################

shape <-  c("diamond = included in h", "diamond = included in h", "diamond = included in h", "diamond = included in h", "diamond = included in h", 
            "circle = not included in h", "circle = not included in h", "circle = not included in h", "diamond = included in h", "circle = not included in h",
            "circle = not included in h", "circle = not included in h", "diamond = included in h", "circle = not included in h", "circle = not included in h")

shape <- factor(shape, levels = unique(shape))

ri_data <- tibble(component = c("pov", "unemp", "pci", "nhsd", "age65", "age17", "disabl", "sngpnt", "mnrty", "limeng", "munit", "mobile", "crowd", "noveh", "grpq"),
                          phi_pm = c(mat[,]), phi_lcl = c(mat.ci[1,]), phi_ucl = c(mat.ci[2, ]),
                   shape = shape)%>%
                  mutate(component_fct = factor(component, levels = unique(component)),
                         comp = paste0('Z', 1:15),
                         zcomp = factor(comp, levels = unique(comp)))



ri_data <- tibble(component = c("pov", "unemp", "pci", "nhsd", "age65", "age17", "disabl", "sngpnt", "mnrty", "limeng", "munit", "mobile", "crowd", "noveh", "grpq"),
                  phi_pm = c(mat[,]), phi_lcl = c(mat.ci[1,]), phi_ucl = c(mat.ci[2, ]))%>%
            mutate(component_fct = factor(component, levels = unique(component)),
                   comp = paste0('Z', 1:15),
                   zcomp = factor(comp, levels = unique(comp)))

plot1 <- ggplot(ri_data, aes(zcomp,  phi_pm, ymin = phi_lcl, ymax = phi_ucl)) + 
          scale_shape_manual(values = c(5, 16))+
          geom_pointrange(position = position_dodge(width = 0.1))+
          geom_hline(yintercept = 1)+
          # geom_point(mapping = aes(zcomp, true.ri, color = 'red', show.legend = FALSE))+
          ylab("Relative importance")+
          ylim(0.5, 1.75)+
          xlab("SVI variable")+
          theme(plot.title = element_text(size = 8, face = 'bold', hjust = 0.5),
                axis.text = element_text(size = 8, face = 'bold'),
                axis.title = element_text(size = 8),
               legend.position = c(0.4, 0.85),
               legend.title = element_blank())




ggsave(
  "relative_importance_null sim.pdf",
  device = pdf(),
  scale = 1,
  dpi = "print"(300),
  width = 8,
  height = 3.5,
  units = "in"
)

dev.off()


