library(reshape2)

####################################################
#####Identifying and removing contaminant OTUs #####
####################################################

setwd("~/Documents/Bioinformatics_scripts/R_scripts/Coral_microbiome_transplants/50_samples_only/Input/")

#read in OTU table generated in mothur
otu = t(read.table("transplants.final.OTU_table", header = FALSE)) 
colnames(otu) = otu[2, ]
otu=otu[-c(1,2,3), ] 

tax = read.table("transplants.final.taxonomy", header = TRUE)
tax$Taxonomy = gsub("[(0-9)]", "", tax$Taxonomy) 
tax$Taxonomy = gsub("[a-z]__", "", tax$Taxonomy) 
tax[c('Kingdom', 'Phylum', 'Class', 'Order', 'Family', 'Genus', 'Species')] <- colsplit(tax$Taxonomy,';',c('Kingdom', 'Phylum', 'Class', 'Order', 'Family', 'Genus', 'Species'))
tax=tax[, -c(3)]

#Identify and removing contaminant OTUs raw data
otu.n=apply(otu[,2:ncol(otu)], 2, as.numeric)
row.names(otu.n)=otu[,1]
otu.r=sweep(otu.n,2,colSums(otu.n),`/`)
#otu.r$negSum = cumsum(otu.r$Neg_PCR)/sum((otu.r$Neg_ext))
otu.r=transform(otu.r,  Sum = rowSums(otu.r[,2:ncol(otu.r)]))
names(otu.r)
otu.r=transform(otu.r,  SumNegs = rowSums(otu.r[,c(26,27)])) # define negative controls here
names(otu.r)
otu.r=transform(otu.r,  contaFactor=(otu.r$SumNegs/otu.r$Sum)*100)
Conta=subset(otu.r, otu.r$contaFactor > 10)
Conta$Family=tax$Family[match(rownames(Conta), tax$OTU)]
cat("number of potentially contaminant OTUs: ", length(rownames(Conta)), "\n", rownames(Conta),"\n", Conta$Family)
#write.table(rownames(Conta), "contaminants.txt",  quote = FALSE, row.names=FALSE, sep = "\t") #define sample range


# Export normalized and raw OTU tables 
OTU.noConta.n=subset(as.data.frame(otu.n), !rownames(otu.n) %in% rownames(Conta)) # not normalized
OTU.noConta.r=subset(as.data.frame(otu.r), !rownames(otu.r) %in% rownames(Conta)) # normalized
names(OTU.noConta.n)
#write.table(OTU.noConta.n[,c(1:25,28:52)], "OTU_noContaminants.raw.txt",  quote = FALSE, row.names=TRUE, sep = "\t") #define sample range
names(OTU.noConta.r)
#write.table(OTU.noConta.r[,c(1:25,28:52)], "OTU_noContaminants.normalized.txt",  quote = FALSE, row.names=TRUE, sep = "\t") #define sample range

# Export CLR transformed data

require(Seurat)
clr.otu=NormalizeData(OTU.noConta.r[,c(1:25,28:52)], normalization.method = "CLR")

# Exclude contaminants from taxa file
tax.noConta=subset(as.data.frame(tax), !tax$OTU %in% rownames(Conta)) 
write.table(tax.noConta, "Tax_noContaminants.txt",  quote = FALSE, row.names=FALSE, sep = "\t") 

## Exclude outliers and negatives from metafile
# meta = read.table("transplants_metadata.txt", header = TRUE)
# met.noConta=subset(meta, meta$Sample %in% names(OTU.noConta.r)) 
# write.table(met.noConta, "meta_50.txt",  quote = FALSE, row.names=FALSE, sep = "\t") 

### make final OTU table
library(phylotools)
fin=OTU.noConta.n[,c(1:25,28:52)]
colnames(fin)=gsub("APO_Aip_7D", "APO+APOinoc",colnames(fin))
colnames(fin)=gsub("APO_CA_7D", "APO+ACRinoc",colnames(fin))
colnames(fin)=gsub("APO_CB_7D", "APO+PORinoc",colnames(fin))
colnames(fin)=gsub("APO_C", "APO",colnames(fin))
colnames(fin)=gsub("SYM_Aip_7D", "SYM+SYMinoc",colnames(fin))
colnames(fin)=gsub("SYM_CA_7D", "SYM+ACRinoc",colnames(fin))
colnames(fin)=gsub("SYM_CB_7D", "SYM+PORinoc",colnames(fin))
colnames(fin)=gsub("SYM_C", "SYM",colnames(fin))
colnames(fin)

fasta=phylotools::read.fasta("transplants.final.fasta")
tax = read.table("transplants.final.taxonomy", header = TRUE)
tax_2=read.table("transplants.final.nr_v138.wang.taxonomy")

fin$Sum=rowSums(fin[,1:50])
fin$GreenGenes=tax$Taxonomy[match(rownames(fin), tax$OTU)]
fin$SILVA=tax_2$V2[match(rownames(fin), tax_2$V1)]
fin$fasta=fasta$seq.text[match(rownames(fin), fasta$seq.name)]
write.table(fin, "Final_OTU_table.txt", quote = F, sep = "\t")
