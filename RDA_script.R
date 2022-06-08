
#install.packages('psych')
library(psych)    
#install.packages('vegan')
library(vegan)  
#install.packages('adegenet')
library(adegenet)
#install.packages('envirem')
library(envirem)

# R version 3.6.1 (2019-07-05) -- "Action of the Toes"

# This directory will need to be changed based on where you put the scripts and data files
setwd("~/my_docs/Carex/hybrid/Carex_hybridization")

# This example is based on the alpha outliers using hybrids + C. nelsonii
# You will need to change n.ind and n.loc as needed when using other 
# combinations of loci and species

input <- read.structure('nel_hyb_alpha_high.stru', 
                        n.ind=101, n.loc = 1027, 
                        onerowperind=F,col.lab=1,col.pop=2,
                        row.marknames=1,NA.char='-9',ask=F)

popz <- genind2genpop(input)
freq <- makefreq(popz,missing='NA')
ncol(freq)
nrow(freq)

# below, the to = value should match the ncol(freq) value from above
s <- seq(from = 1, to = 1547,by=2)
freq_l1 <- freq[,s]

# this following input file will need to be changed if different 
# combinations of species are used
new_env_carex_pop <- read.table("nel_hyb_pca_notprojected_ordered.txt"
                                ,header=T)
colnames(new_env_carex_pop) <- c("Population", "PC1", "PC2", "PC3")

pred_carex_pop <- subset(new_env_carex_pop, select=c(PC1,PC2,PC3))
pairs.panels(pred_carex_pop, scale=T)

nel_nov.rda_pop <- rda(freq_l1~ ., data=pred_carex_pop, scale=T)
nel_nov.rda_pop

RsquareAdj(nel_nov.rda_pop)
# these are the values you should get when using 'nel_hyb_alpha_high.stru'
# and 'nel_hyb_pca_notprojected_ordered.txt'

#$r.squared
#[1] 0.4313525

#$adj.r.squared
#[1] 0.2418033

signif.full_carex <- anova.cca(nel_nov.rda_pop, parallel=getOption("mc.cores")) # default is permutation=999
signif.full_carex     

signif.axis_carex <- anova.cca(nel_nov.rda_pop, by="axis", parallel=getOption("mc.cores"))
signif.axis_carex

bg_carex_pop <- c("black","purple3","white",
                  "gray65","mediumpurple1",
                  "red","red","red","red",
                  "red","red","red","red")

eco_carex_pop <- new_env_carex_pop$Population
bg_carex_pop[eco_carex_pop]

pdf("high_RDA12_nel_hyb_PC1-PC2-PC3.pdf", height=7,width=10)
plot(nel_nov.rda_pop, type="n", scaling=3)
points(nel_nov.rda_pop, display="species", pch=20, cex=0.7, col="gray32", scaling=3)           # the SNPs
points(nel_nov.rda_pop, display="sites", pch=21, cex=1.5, col="gray32", scaling=3, bg=bg_carex_pop[eco_carex_pop]) # the pops
text(nel_nov.rda_pop, scaling=3, display="bp", col="#0868ac", cex=1)                           # the predictors
legend("bottomright", legend=levels(eco_carex_pop), bty="n", col="gray32", pch=21, cex=0.7, pt.bg=bg_carex_pop)
title("High alphas")
dev.off()

#####

load.rda_pop <- scores(nel_nov.rda_pop, choices=c(1:3), display="species")  # Species scores for the first three constrained axes

# the following file will be used in the other script, 'Loadings_script.R'
write.table(load.rda_pop, file="loadings_high_loci_rda_nel+hyb.txt", row.names = T,quote=F)



