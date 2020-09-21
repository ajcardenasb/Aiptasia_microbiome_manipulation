library(data.table)
library(VennDiagram)
library(gridExtra)

setwd("~/Documents/Bioinformatics_scripts/R_scripts/Coral_microbiome_transplants/50_samples_only/Input/")

#1. import data
otu = read.table("OTU_noContaminants.raw.txt", header = TRUE)  # Replace with your file name, using Mothur's output directly (e.g tab-delimited file)
otu$otu=rownames(otu)
otu.l = reshape2::melt(otu, id.vars = c("otu"))
pres=subset(otu.l, value >= 3)
met = read.table("meta_50.txt", header = TRUE)

#2. consider OTUs that are present in all samples per group

apoT.m=otu[grepl("APO_C_", names(otu))]
symT.m=otu[grepl("SYM_C_", names(otu))]
inocApo.m=otu[grepl("APO_Aip_7D", names(otu))]
inocSym.m=otu[grepl("SYM_Aip_7D_", names(otu))]
inocPorApo.m=otu[grepl("APO_CB_7D_", names(otu))]
inocPorSym.m=otu[grepl("SYM_CB_7D_", names(otu))]
inocAcrApo.m=otu[grepl("APO_CA_7D_", names(otu))]
inocAcrSym.m=otu[grepl("SYM_CA_7D_", names(otu))]

apoT.f=apoT.m[apply(apoT.m, MARGIN = 1, function(x) all(x > 0)), ]
symT.f=symT.m[apply(symT.m, MARGIN = 1, function(x) all(x > 0)), ]
inocApo.f=inocApo.m[apply(inocApo.m, MARGIN = 1, function(x) all(x > 0)), ]
inocSym.f=inocSym.m[apply(inocSym.m, MARGIN = 1, function(x) all(x > 0)), ]
inocPorApo.f=inocPorApo.m[apply(inocPorApo.m, MARGIN = 1, function(x) all(x > 0)), ]
inocPorSym.f=inocPorSym.m[apply(inocPorSym.m, MARGIN = 1, function(x) all(x > 0)), ]
inocAcrApo.f=inocAcrApo.m[apply(inocAcrApo.m, MARGIN = 1, function(x) all(x > 0)), ]
inocAcrSym.f=inocAcrSym.m[apply(inocAcrSym.m, MARGIN = 1, function(x) all(x > 0)), ]


#3.exporting lists

merge(apoT.f, symT.f, inocApo.f, inocSym.f, inocPorApo.f, inocPorSym.f, inocAcrApo.f, inocAcrSym.f)
write.table(unique(rownames(apoT.f)), "Apo_tank_list.txt", quote = FALSE, row.names = FALSE)
write.table(unique(rownames(symT.f)), "Sym_tank_list.txt", quote = FALSE, row.names = FALSE)
write.table(unique(rownames(inocApo.f)), "Apo_inoc_Apo_list.txt", quote = FALSE, row.names = FALSE)
write.table(unique(rownames(inocSym.f)), "Sym_inoc_Sym_list.txt", quote = FALSE, row.names = FALSE)
write.table(unique(rownames(inocPorApo.f)), "Apo_inoc_Por_list.txt", quote = FALSE, row.names = FALSE)
write.table(unique(rownames(inocPorSym.f)), "Sym_inoc_Por_list.txt", quote = FALSE, row.names = FALSE)
write.table(unique(rownames(inocAcrApo.f)), "Apo_inoc_Acr_list.txt", quote = FALSE, row.names = FALSE)
write.table(unique(rownames(inocAcrSym.f)), "Sym_inoc_Acr_list.txt", quote = FALSE, row.names = FALSE)
                                
                                
## 4. Creating table

tax = read.table("Tax_noContaminants.txt", header = TRUE,  row.names=1)
tab=data.frame(OTUs = unique(c(rownames(apoT.f),rownames(symT.f), rownames(inocApo.f), rownames(inocSym.f), rownames(inocPorApo.f), rownames(inocPorSym.f), rownames(inocAcrApo.f), rownames(inocAcrSym.f))))
tab$APO=ifelse(tab$OTUs %in% rownames(apoT.f), "present", "not present")
tab$APO_APOinoc=ifelse(tab$OTUs %in% rownames(inocApo.f) , "present", "not present")
tab$APO_ACRinoc=ifelse(tab$OTUs %in%  rownames(inocAcrApo.f) , "present", "not present")
tab$APO_PORinoc=ifelse(tab$OTUs %in% rownames(inocPorApo.f), "present", "not present")
tab$SYM=ifelse(tab$OTUs %in% rownames(symT.f) , "present", "not present")
tab$SYM_SYMinoc=ifelse(tab$OTUs %in% rownames(inocSym.f) , "present", "not present")
tab$SYM_ACRinoc=ifelse(tab$OTUs %in% rownames(inocAcrSym.f) , "present", "not present")
tab$SYM_PORinoc=ifelse(tab$OTUs %in% rownames(inocPorSym.f) , "present", "not present")
tab$tax=paste(tax$Genus, tax$Family, sep = "-")[match(tab$OTUs, rownames(tax))]
write.table(tab, "VennDiagram_table.txt", quote = F, row.names=F, sep = "\t" )


###5.  Ploting Venn diagrams

vennD1=venn.diagram(list(rownames(apoT.f), rownames(inocApo.f)), NULL, category.names = c("APO", "APO+APOinoc"), rotation.degree = 90, fill = c("#e69f00", "#f0e442"), lty="blank", cat.pos = c(180,0),  alpha = rep(0.7, 2)) #cat.fontfamily = "Arial",  fontfamily="Arial" ,
vennD2=venn.diagram(list(rownames(symT.f), rownames(inocSym.f)), NULL, category.names = c("SYM", "SYM+SYMinoc"), rotation.degree = 90, fill = c("#e69f00", "#f0e442"), lty="blank", cat.pos = c(180,0) , alpha = rep(0.7, 2))
vennD3=venn.diagram(list(rownames(apoT.f),rownames(inocAcrSym.f),rownames(symT.f), rownames(inocAcrApo.f) ),NULL, category.names = c("APO", "SYM+ACRinoc","SYM", "APO+ACRinoc" ),  fill = c("#e69f00", "#0072b2", "#e69f00", "#0072b2"  ), lty="blank")
vennD4=venn.diagram(list(rownames(apoT.f),  rownames(inocPorSym.f),rownames(symT.f),rownames(inocPorApo.f)),NULL, category.names = c("APO", "SYM+PORinoc","SYM" , "APO+PORinoc"), fill = c("#e69f00",   "#0072b2", "#e69f00", "#0072b2" ), lty="blank")

pdf("../outputs/transplants50_venn2.pdf", width=25,height=5, pointsize = 16)
grid.arrange(grobTree(vennD1), grobTree(vennD2),grobTree(vennD3), grobTree(vennD4), ncol = 5, nrow = 1)
dev.off()



