
#setwd("~/Documents/Projects/clam_microbiome/") #change to specify the location of your directory

#1. import data
otu = read.table("OTU_noContaminants.raw.txt", header = TRUE)  # Replace with your file name, using Mothur's output directly (e.g tab-delimited file)
otu$otu=rownames(otu)
#2. create lists
#otu[,2:ncol(otu)] = as.numeric(as.character(otu[,2:ncol(otu)]))
#otu.df[2:ncol(otu.df)] <- lapply(otu.df[2:ncol(otu.df)], as.numeric) # transform factors into numeric
otu.l = melt(otu, id.vars = c("otu"))
pres=subset(otu.l, value > 0)
met = read.table("transplants_metadata.txt", header = TRUE)
pres$condition1=met$new[match(pres$variable, met$Sample)]# check what is the column you wan to add in your metadata and replace V1 and V2
pres$condition2=met$Symbiosis[match(pres$variable, met$Sample)]#
grp1=subset(pres, condition1 == "3D-CB" & condition2 == "APO") ## add here your group of interest
grp2=subset(pres, condition1 == "3D-CB" & condition2 == "SYM")
grp3=subset(pres, condition1 == "Inoc" & condition2 == "B") 
grp4=subset(pres, condition1 == "3D-CA" & condition2 == "APO") 
grp5=subset(pres, condition1 == "3D-CA" & condition2 == "SYM") 
grp6=subset(pres, condition1 == "Inoc" & condition2 == "A") 
grp7=subset(pres, condition1 == "3D-Aip" & condition2 == "APO") 
grp8=subset(pres, condition1 == "3D-Aip" & condition2 == "SYM") 
grp9=subset(pres, condition1 == "Inoc" & condition2 == "APO") 
grp10=subset(pres, condition1 == "Inoc" & condition2 == "SYM") 

#3. plot Venn diagrams
library(gplots)

pdf("transplants_venn.pdf", onefile = TRUE, width=15,height=15)
venn(list(APO=unique(grp1$otu), SYM=unique(grp2$otu), Inoc=unique(grp3$otu)))
venn(list(APO=unique(grp4$otu), SYM=unique(grp5$otu), Inoc=unique(grp6$otu)))
venn(list(APO=unique(grp7$otu), SYM=unique(grp8$otu), InocAPO=unique(grp9$otu), InocSYM=unique(grp10$otu)))

dev.off()


#4. (optional) exporting lists

write.table(grp1, "grp1_list.txt", quote = FALSE, row.names = FALSE)
write.table(grp2, "grp2_list.txt", quote = FALSE, row.names = FALSE) 
write.table(grp3, "grp3_list.txt", quote = FALSE, row.names = FALSE)
write.table(grp4, "grp4_list.txt", quote = FALSE, row.names = FALSE)