
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

# Load in the RDA results
all_load <- read.table("loadings_all_loci_rda_nel+hyb.txt", header=T)
low_load <- read.table("loadings_low_loci_rda_nel+hyb.txt", header=T)
high_load <- read.table("loadings_high_loci_rda_nel+hyb.txt", header=T)
in_load <- read.table("loadings_inlier_loci_rda_nel+hyb.txt", header=T)

# Load in the candidate adaptive loci -- outliers with respect to the 
# RDA space of inlier loci
in_cand <- read.table("cand_pop_inlier_nel+hyb.txt", header=T)    
head(in_cand)
head(in_load)
in_merge <- merge(in_cand,in_load, all=T)

low_cand <- read.table("cand_pop_low_nel+hyb.txt", header=T)    
head(low_cand)
head(low_load)
low_merge <- merge(low_cand,low_load, all=T)

high_cand <- read.table("cand_pop_high_nel+hyb.txt", header=T)    
head(high_cand)
head(high_load)
high_merge <- merge(high_cand,high_load, by="snp")

all_cand <- read.table("cand_pop_all_nel+hyb.txt", header=T)


load_merge <- data.frame(rbind(in_merge,high_merge,low_merge))

# Keep only loci that are at least 2 standard deviations from the 
# mean of the inlier RDA space

high_sig <- subset(high_merge, high_merge$loading>=2*sd(in_load$RDA1) |   
                     high_merge$loading<=(-2*sd(in_load$RDA1)))
low_sig <- subset(low_merge, low_merge$loading>=2*sd(in_load$RDA1) |      
                    low_merge$loading<=(-2*sd(in_load$RDA1)))

high_sig1 <- subset(high_sig, high_sig$predictor=="PC1") 
high_sig2 <- subset(high_sig, high_sig$predictor=="PC2") 
high_sig3 <- subset(high_sig, high_sig$predictor=="PC3")

low_sig1 <- subset(low_sig, low_sig$predictor=="PC1") 
low_sig2 <- subset(low_sig, low_sig$predictor=="PC2") 
low_sig3 <- subset(low_sig, low_sig$predictor=="PC3")

plot.new()
par(oma=c(1,1,1,1));
par(mar=c(6,8,3,3));
plot.window(xlim=c(-0.4,0.7), ylim=c(-0.35,0.35))

points(in_load$RDA2~in_load$RDA1, pch=21,col="black", bg="white", cex=0.8)

points(high_sig2$RDA2~high_sig2$RDA1, pch=21, col="black",bg="red") 
points(low_sig2$RDA2~low_sig2$RDA1, pch=21, col="black",bg="darkblue")
points(high_sig1$RDA2~high_sig1$RDA1, pch=21, col="black",bg="maroon") 
points(low_sig1$RDA2~low_sig1$RDA1, pch=21, col="black",bg="cornflowerblue") 

# This variable is defined in the other script 'RDA_script.R'
text(nel_nov.rda_pop, scaling=3, display="bp", col="blue", cex=1)      # the predictors

abline(v= 3*sd(in_load$RDA1), col='black', lty=2)
abline(v= -3*sd(in_load$RDA1), col='black', lty=2)
abline(v= 2*sd(in_load$RDA1), col='black', lty=2)
abline(v= -2*sd(in_load$RDA1), col='black', lty=2)
abline(v= 2.5*sd(in_load$RDA1), col='black', lty=2)
abline(v= -2.5*sd(in_load$RDA1), col='black', lty=2)
axis(1, at=seq(-0.4, 0.4, by=0.1), cex.axis=1.15);
axis(2, at=seq(-0.35, 0.35, by=0.1), cex.axis=1.15);
mtext(side=1, text='RDA1',line=3, cex=1.25, adj = 0.4)
mtext(side=2, text='RDA2', line=4, cex=1.25)


quartz.save(file="SNPs_locus_corr_plot.pdf", type = "pdf", 
            height=7, width=10, dpi=300)
