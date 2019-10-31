library(data.table)


#1. import data
otu = read.table("otu_noConta_noOutliers.txt", header = TRUE)  # Replace with your file name, using Mothur's output directly (e.g tab-delimited file)
otu$otu=rownames(otu)
otu.l = melt(otu, id.vars = c("otu"))
pres=subset(otu.l, value >= 3)
met = read.table("transplants_metadata_noOutliers.txt", header = TRUE)

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


#3. checking plot Venn diagrams, final ones were created on http://bioinformatics.psb.ugent.be/webtools/Venn/

library(VennDiagram)
library(Vennerable)
library(gridExtra)

#color set 1
vennD1=venn.diagram(list(rownames(apoT.f), rownames(inocApo.f)), NULL, category.names = c("APO", "APO+APOinoc"), rotation.degree = 90, fill = c("#A3C38E", "#577545"), lty="blank", cat.pos = c(180,0), cat.fontfamily = "Arial",  fontfamily="Arial" , alpha = rep(0.7, 2))
vennD2=venn.diagram(list(rownames(symT.f), rownames(inocSym.f)), NULL, category.names = c("SYM", "SYM+SYMinoc"), rotation.degree = 90, fill = c("#A3C38E", "#577545"), lty="blank", cat.pos = c(180,0), cat.fontfamily = "Arial",  fontfamily="Arial" , alpha = rep(0.7, 2))
vennD3=venn.diagram(list(rownames(apoT.f), rownames(inocAcrApo.f)), NULL, category.names = c("APO", "APO+ACRinoc"), rotation.degree = 90, fill = c("#A3C38E", "#39597C"), lty="blank", cat.pos = c(180,0), cat.fontfamily = "Arial",  fontfamily="Arial" , alpha = rep(0.7, 2))
vennD4=venn.diagram(list(rownames(inocAcrSym.f), rownames(symT.f)), NULL, category.names = c("SYM", "SYM+ACRinoc"), rotation.degree = 90, fill = c("#A3C38E", "#39597C"), lty="blank", cat.pos = c(180,0), cat.fontfamily = "Arial",  fontfamily="Arial" , alpha = rep(0.7, 2))
vennD5=venn.diagram(list(rownames(apoT.f), rownames(inocPorApo.f)), NULL, category.names = c("APO", "APO+PORinoc"), rotation.degree = 90, fill = c("#A3C38E", "#7E9EC1"), lty="blank", cat.pos = c(180,0), cat.fontfamily = "Arial",  fontfamily="Arial" , alpha = rep(0.7, 2))
vennD6=venn.diagram(list(rownames(symT.f), rownames(inocPorSym.f)), NULL, category.names = c("SYM", "SYM+PORinoc"), rotation.degree = 90, fill = c("#A3C38E", "#7E9EC1"), lty="blank", cat.pos = c(180,0), cat.fontfamily = "Arial",  fontfamily="Arial" , alpha = rep(0.7, 2))

#color set 2
vennD1=venn.diagram(list(rownames(apoT.f), rownames(inocApo.f)), NULL, category.names = c("APO", "APO+APOinoc"), rotation.degree = 90, fill = c("#FF6549", "#B82B12"), lty="blank", cat.pos = c(180,0), cat.fontfamily = "Arial",  fontfamily="Arial" , alpha = rep(0.7, 2))
vennD2=venn.diagram(list(rownames(symT.f), rownames(inocSym.f)), NULL, category.names = c("SYM", "SYM+SYMinoc"), rotation.degree = 90, fill = c("#FF6549", "#B82B12"), lty="blank", cat.pos = c(180,0), cat.fontfamily = "Arial",  fontfamily="Arial" , alpha = rep(0.7, 2))
vennD3=venn.diagram(list(rownames(apoT.f), rownames(inocAcrApo.f)), NULL, category.names = c("APO", "APO+ACRinoc"), rotation.degree = 90, fill = c("#FF6549", "#39597C"), lty="blank", cat.pos = c(180,0), cat.fontfamily = "Arial",  fontfamily="Arial" , alpha = rep(0.7, 2))
vennD4=venn.diagram(list(rownames(inocAcrSym.f), rownames(symT.f)), NULL, category.names = c("SYM", "SYM+ACRinoc"), rotation.degree = 90, fill = c("#FF6549", "#39597C"), lty="blank", cat.pos = c(180,0), cat.fontfamily = "Arial",  fontfamily="Arial" , alpha = rep(0.7, 2))
vennD5=venn.diagram(list(rownames(apoT.f), rownames(inocPorApo.f)), NULL, category.names = c("APO", "APO+PORinoc"), rotation.degree = 90, fill = c("#FF6549E", "#7E9EC1"), lty="blank", cat.pos = c(180,0), cat.fontfamily = "Arial",  fontfamily="Arial" , alpha = rep(0.7, 2))
vennD6=venn.diagram(list(rownames(symT.f), rownames(inocPorSym.f)), NULL, category.names = c("SYM", "SYM+PORinoc"), rotation.degree = 90, fill = c("#FF6549", "#7E9EC1"), lty="blank", cat.pos = c(180,0), cat.fontfamily = "Arial",  fontfamily="Arial" , alpha = rep(0.7, 2))


png("transplants_venn.png", width=1264,height=250, res = 300, pointsize = 4)
svg("transplants_venn.svg", width = 10, height = 3, pointsize = 12)
grid.arrange(grobTree(vennD1), grobTree(vennD2),grobTree(vennD3), grobTree(vennD4),grobTree(vennD5), grobTree(vennD6),  ncol = 6, nrow = 1)
dev.off()

#4.exporting lists
# 
# write.table(unique(rownames(apoT.f)), "Apo_tank_list.txt", quote = FALSE, row.names = FALSE)
# write.table(unique(rownames(symT.f)), "Sym_tank_list.txt", quote = FALSE, row.names = FALSE)
# write.table(unique(rownames(inocApo.f)), "Apo_inoc_Apo_list.txt", quote = FALSE, row.names = FALSE)
# write.table(unique(rownames(inocSym.f)), "Sym_inoc_Sym_list.txt", quote = FALSE, row.names = FALSE)
# write.table(unique(rownames(inocPorApo.f)), "Apo_inoc_Por_list.txt", quote = FALSE, row.names = FALSE)
# write.table(unique(rownames(inocPorSym.f)), "Sym_inoc_Por_list.txt", quote = FALSE, row.names = FALSE)
# write.table(unique(rownames(inocAcrApo.f)), "Apo_inoc_Acr_list.txt", quote = FALSE, row.names = FALSE)
# write.table(unique(rownames(inocAcrSym.f)), "Sym_inoc_Acr_list.txt", quote = FALSE, row.names = FALSE)
