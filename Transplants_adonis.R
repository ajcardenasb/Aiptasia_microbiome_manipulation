library(pairwiseAdonis)
library(vegan)

setwd("~/Documents/Bioinformatics_scripts/R_scripts/Coral_microbiome_transplants/50_samples_only/Input/")

##data in
otu = as.data.frame(read.table("OTU_noContaminants.raw.txt", header = TRUE,  row.names=1))
otu= as.data.frame(t(sweep(otu,2,colSums(otu),"/")))
met = read.table("meta_50.txt", header = TRUE,  row.names=1)

#### full model all
otu$symbiosis=met$Symbiosis[match(rownames(otu), rownames(met))]
otu$treatment=met$new[match(rownames(otu), rownames(met))]
adonis(otu[,1:448]~  otu$treatment +otu$symbiosis , p.adjust.m ='fdr', sim.method = 'bray')
adonis(otu[,1:448]~ otu$symbiosis+ otu$treatment  , p.adjust.m ='fdr', sim.method = 'bray')

#### pairwise adonis all
otu$adonis=paste(met$Symbiosis,met$barplots, sep = "_")[match(rownames(otu), rownames(met))]
ado.all=pairwise.adonis(otu[,1:448],otu$adonis, p.adjust.m ='fdr', sim.method = 'bray')
write.table(ado.all, "../outputs/adonis_50.txt", quote = FALSE, row.names = FALSE)
pairwise.adonis(otu[,1:448],otu$adonis, p.adjust.m ='fdr', sim.method = 'euclidean')
