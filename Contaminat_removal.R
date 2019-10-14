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
