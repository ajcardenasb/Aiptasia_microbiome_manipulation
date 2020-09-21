library(ggplot2)
library(dplyr)
library(phyloseq)
library(gridExtra)
library(plyr)

setwd("~/Documents/Bioinformatics_scripts/R_scripts/Coral_microbiome_transplants/50_samples_only/Input/")

############################################################
##################### phyloseq #############################
############################################################

otu = read.table("OTU_noContaminants.raw.txt", header = TRUE,  row.names=1)
tax = read.table("Tax_noContaminants.txt", header = TRUE,  row.names=1)
met = read.table("meta_50.txt", header = TRUE,  row.names=1)
met$final_names=factor(met$final_names, levels = c("APO", "APO_AB", "APO+APOinoc", "APO+ACRinoc", "APO+PORinoc", "SYM", "SYM_AB", "SYM+SYMinoc", "SYM+ACRinoc", "SYM+PORinoc"))

##create Phyloseq object

otu.t= otu_table(otu, taxa_are_rows=TRUE)
sam.t= sample_data(met)
tax.t= tax_table(as.matrix(tax))

phy.all= phyloseq(otu.t, tax.t,  sam.t)

##transform data and subset phyloseq objects
phy.t=microbiome::transform(phy.all, transform = "log10", target = "OTU", shift = 0, scale = 1)
apo.phy=subset_samples(phy.t, Symbiosis == "APO")
sym.phy=subset_samples(phy.t, Symbiosis == "SYM")


##Plot ordination

COL=c("#e69f00",  "#00a075", "#f0e442", "#0072b2", "#cc79a7") 

apo.PCOA_br = ordinate(apo.phy, method = "PCoA", distance = "bray")
a = plot_ordination(apo.phy, apo.PCOA_br , color = "final_names", )  + geom_point(size = 3, alpha = 1) + theme_bw()  + theme(plot.title = element_text(hjust = 0.5)) + scale_colour_manual(values=COL)+ theme(legend.title = element_blank()) + ggtitle("APO") 

sym.PCOA_br = ordinate(sym.phy, method = "PCoA", distance = "bray")
b = plot_ordination(sym.phy, sym.PCOA_br , color = "final_names", )  + geom_point(size = 3, alpha = 1) + theme_bw()  + theme(plot.title = element_text(hjust = 0.5)) + scale_colour_manual(values=COL)+ theme(legend.title = element_blank()) + ggtitle("SYM") 

pdf("../outputs/transplants50_PCoAs.pdf", onefile = TRUE, width=10,height=5, pointsize = 12)
grid.arrange(a, b, ncol = 2)
dev.off()

