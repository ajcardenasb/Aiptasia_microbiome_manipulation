###################################################################
#####Identifying and removing contaminant OTUs normalized data#####
###################################################################

#read in OTU table generated in mothur
otu = t(read.table("transplants.final.OTU_table", header = FALSE)) 
colnames(otu) = otu[2, ]
otu=otu[-c(1,2,3), ] 

tax = read.csv("transplants.final.taxonomy", header = TRUE, sep = "\t")
tax$Taxonomy = gsub("[(0-9)]", "", tax$Taxonomy) 
tax$Taxonomy = gsub("[a-z]__", "", tax$Taxonomy) 
tax[c('Kingdom', 'Phylum', 'Class', 'Order', 'Family', 'Genus', 'Species')] <- colsplit(tax$Taxonomy,';',c('Kingdom', 'Phylum', 'Class', 'Order', 'Family', 'Genus', 'Species'))
tax=tax[, -c(3)]

#Identify and removing contaminant OTUs raw data
otu.n=apply(otu[,2:ncol(otu)], 2, as.numeric)
row.names(otu.n)=otu[,1]
otu.r=as.data.frame(sweep(otu.n,2,colSums(otu.n),`/`))
#otu.r$negSum = cumsum(otu.r$Neg_PCR)/sum((otu.r$Neg_ext))
otu.r=transform(otu.r,  Sum = rowSums(otu.r[,2:ncol(otu.r)]))
names(otu.r)
otu.r=transform(otu.r,  SumNegs = rowSums(otu.r[,c(54,55)])) # define negative controls here
names(otu.r)
otu.r=transform(otu.r,  contaFactor=(otu.r$SumNegs/otu.r$Sum)*100)
Conta=subset(otu.r, otu.r$contaFactor > 10)
Conta$Family=tax$Family[match(rownames(Conta), tax$OTU)]
cat("number of potentially contaminant OTUs: ", length(rownames(Conta)), "\n", rownames(Conta),"\n", Conta$Family)

# Export normalized and raw OTU tables 
OTU.noConta.n=subset(as.data.frame(otu.n), !rownames(otu.n) %in% rownames(Conta)) # not normalized
OTU.noConta.r=subset(as.data.frame(otu.r), !rownames(otu.r) %in% rownames(Conta)) # normalized
write.table(OTU.noConta.n[,c(1:53,56:107)], "OTU_noContaminants.raw.txt",  quote = FALSE, row.names=TRUE, sep = "\t") #define sample range
write.table(OTU.noConta.r[,c(1:53,56:107)], "OTU_noContaminants.normalized.txt",  quote = FALSE, row.names=TRUE, sep = "\t") #define sample range

# Exclude contaminants from taxa file
tax.noConta=subset(as.data.frame(tax), !tax$OTU %in% rownames(Conta)) 
write.table(tax.noConta, "Tax_noContaminants.txt",  quote = FALSE, row.names=FALSE, sep = "\t") 


###########################################
#####Identifying and removing outliers#####
###########################################

library("phyloseq")
library(gridExtra)
library(plyr)

otu = read.table("OTU_noContaminants.raw.txt", header = TRUE,  row.names=1)
tax = read.table("Tax_noContaminants.txt", header = TRUE,  row.names=1)
met = read.table("transplants_metadata.txt", header = TRUE,  row.names=1)

tax.1=tax[,-1]
otu.t= otu_table(otu, taxa_are_rows=TRUE)
sam.t= sample_data(met)
tax.t= tax_table(as.matrix(tax.1))

phy.all= phyloseq(otu.t, tax.t,  sam.t)

subset=subset_samples(phy.all, names(otu) %like% "APO_C_")#6
subset=subset_samples(phy.all, names(otu) %like% "SYM_C_")#1 or 3
subset=subset_samples(phy.all, names(otu) %like% "APO_Aip_7D")# 4 or 5
subset=subset_samples(phy.all, names(otu) %like% "SYM_Aip_7D_")#1
subset=subset_samples(phy.all, names(otu) %like% "APO_CB_7D_")#6,2 or 4
subset=subset_samples(phy.all, names(otu) %like% "SYM_CB_7D_")#4
subset=subset_samples(phy.all, names(otu) %like% "APO_CA_7D_")#2
subset=subset_samples(phy.all, names(otu) %like% "SYM_CA_7D_")#5
subset=subset_samples(phy.all, names(otu) %like% "APO_AB_")#2
subset=subset_samples(phy.all, names(otu) %like% "SYM_AB_")#1


PCOA = ordinate(subset, method = "PCoA", distance = "bray") 
plot_ordination(subset, PCOA ,  label="Replicate", )  + geom_text(mapping = aes(label = Replicate), size = 5) 

outliers=c("APO_C_6", "SYM_C_1","APO_Aip_7D_5", "SYM_Aip_7D_1", "APO_CB_7D_6", "SYM_CB_7D_4", "APO_CA_7D_2", "SYM_CA_7D_5", "APO_AB_2", "SYM_AB_1")

# Exclude outliers from otu and met file
otu.noOut=otu[, !names(otu) %in% outliers]
write.table(otu.noOut, "otu_noConta_noOutliers.txt",  quote = FALSE, row.names=TRUE, sep = "\t") 
met=read.table("transplants_metadata.txt", header = TRUE)
met.noOut=subset(met, !met$Sample %in% outliers)
write.table(met.noOut, "transplants_metadata_noOutliers.txt",  quote = FALSE, row.names=FALSE, sep = "\t") 


otu2 = read.table("OTU_noContaminants.normalized.txt", header = TRUE,  row.names=1)
otu2.noOut=otu2[, !names(otu2) %in% outliers]
write.table(otu2.noOut, "otu_normalized_noConta_noOutliers.txt",  quote = FALSE, row.names=TRUE, sep = "\t") 
