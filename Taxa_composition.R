library(reshape2)
library(ggplot2)
library(scales)
library(gridExtra)
library(stringr)

#############################################################
#####Taxonomic profiles of the 20 most abundant families#####
#############################################################

otu = read.table("otu_noConta_noOutliers.txt", header = TRUE,  row.names=1)
tax = read.table("Tax_noContaminants.txt", header = TRUE,  row.names=1)

#########################################################
#####Taxonomic profiles of the 20 most abundant OTUs#####
#########################################################
otu.n=apply(otu[,1:ncol(otu)], 2, as.numeric)
names=rownames(otu)
rownames(otu.n)=names
otu.n=transform(otu.n, OTU = rownames(otu))
names(otu.n)
otu.n2=otu.n[,c(!grepl("Inoc|water|Slurry",names(otu.n))), drop=F]# exclude slurries and water and others
names(otu.n2)

topOTUs=otu.n2[order(rowSums(otu.n2[, c(1:86)]),decreasing = TRUE),][1:20,87] 
otu.top=subset(otu.n2, rownames(otu.n) %in% topOTUs) 
otu.bot=subset(otu.n2, !rownames(otu.n) %in% topOTUs) 
otu.bot$OTU=gsub(".*","Others", otu.bot$OTU)
names(otu.bot)
others=aggregate(otu.bot[, c(1:86)], by = list(OTU=otu.bot[, 87]), FUN =  sum)
names(others)
others.2=others[, c(2:87,1)]
all.2 =rbind(otu.top, others)
all.l=melt(all.2, id.vars=c("OTU"), variable.name = "OTU", value.name = "Abundance")
colnames(all.l)=c("OTU","Sample","Abundance")

## Add sample information
tax = read.table("Tax_noContaminants.txt", header = TRUE,  row.names=1)
meta = read.table("transplants_metadata_noOutliers.txt", header = TRUE)
all.l$Symbiosis=meta$Symbiosis[match(all.l$Sample, meta$Sample)]
all.l$plot=meta$barplots[match(all.l$Sample, meta$Sample)]
all.l$Replicate=meta$Replicate2[match(all.l$Sample, meta$Sample)]
all.l$genus=tax$Genus[match(all.l$OTU, rownames(tax))]
all.l$taxa=paste(all.l$OTU,all.l$genus, sep = " - ")
all.l.por=subset(all.l, !plot == "ignore")
all.l.por$plot=factor(all.l.por$plot, levels = c("Tank", "AB", "Aiptasia", "Acropora", "Porites"))


P21=c( "#C0C0C0", "#771155", "#AA4488", "#CC99BB", "#114477", "#4477AA", "#77AADD", "#117777", "#44AAAA", "#77CCCC", "#117744", "#44AA77", "#88CCAA", "#777711", "#AAAA44", "#DDDD77", "#774411", "#AA7744", "#DDAA77", "#771122", "#AA4455")
bar=ggplot() +geom_bar(aes(y = Abundance, x = Replicate, fill = taxa), data = all.l.por, stat="identity", position = "fill") +  scale_y_continuous(labels = percent_format(), expand = c(0, 0)) + theme(axis.text.x=element_text(angle=90,hjust=1)) + labs( y= "Percentage of 16S rRNA sequences", x="", title= "") + scale_fill_manual(values=P21)  + theme( legend.key.size = unit(0.2, "cm"),legend.key.width = unit(0.2,"cm")) + facet_grid(Symbiosis~plot, space = "free")+ theme_bw() +  theme(strip.background =element_rect(fill="white"), text = element_text(size=18))

pdf("transplants_otu.pdf", onefile = TRUE, width=15,height=7)
plot(bar)
dev.off()



# otu$Family=tax$Family[match(rownames(otu),rownames(tax))]
# names(otu)
# otu.ag=aggregate(otu[, 1:97], by = list(otu[, 98]), FUN =  sum) #define sample range and group factor
# topFamilies=otu.ag[order(rowSums(otu.ag[, 2:ncol(otu.ag)]),decreasing = TRUE),][1:20,1]
# fam.top=subset(otu.ag, otu.ag$Group.1 %in% topFamilies) 
# fam.bot=subset(otu.ag, !otu.ag$Group.1 %in% topFamilies) 
# fam.bot$Group.1=gsub(".*","Others", fam.bot$Group.1)
# others=aggregate(fam.bot[, 2:ncol(fam.bot)], by = list(fam.bot[, 1]), FUN =  sum)
# all.2 =rbind(fam.top, others)
# all.l=melt(all.2, id.vars=c("Group.1"), variable.name = "Family", value.name = "Abundance")
# colnames(all.l)=c("Family","Sample","Abundance")
# 
# ## Add sample information
# meta = read.table("transplants_metadata_noOutliers.txt", header = TRUE)
# all.l$Symbiosis=meta$Symbiosis[match(all.l$Sample, meta$Sample)]
# all.l$Condition=meta$Condition[match(all.l$Sample, meta$Sample)]
# all.l$Replicate=meta$Replicate[match(all.l$Sample, meta$Sample)]
# all.l$Time=meta$Time[match(all.l$Sample, meta$Sample)]
# #all.l.apo.por=subset(all.l, Symbiosis == "APO" & Condition == "CB" | Condition == "AB" | Condition == "C")
# #all.l.sym.por=subset(all.l, Symbiosis == "SYM" & Condition == "CB" | Condition == "AB" | Condition == "C")
# all.l.por=subset(all.l, Condition == "CB" | Condition == "AB" | Condition == "C" )
# all.l.acr=subset(all.l, Condition == "CA" | Condition == "AB" | Condition == "C" )
# all.l.aip=subset(all.l, Condition == "Aip" | Condition == "AB" | Condition == "C" )
# 
# 
# ## Plot
# P21=c("#771155", "#AA4488", "#CC99BB", "#114477", "#4477AA", "#77AADD", "#117777", "#44AAAA", "#77CCCC", "#117744", "#44AA77", "#88CCAA", "#777711", "#AAAA44", "#DDDD77", "#774411", "#AA7744", "#DDAA77", "#771122", "#AA4455", "#C0C0C0")
# 
# #png(filename = "Transplants_barplots_POR.png", res = 300, width = 2500, height = 1200, units = "px")
# #svg(filename = "Transplants_barplots_POR.svg",  width = 10, height = 5, pointsize = 12) 
# 
# por=ggplot() +geom_bar(aes(y = Abundance, x = Replicate, fill = Family), data = all.l.por, stat="identity", position = "fill") +  scale_y_continuous(labels = percent_format(), expand = c(0, 0)) + theme(axis.text.x=element_text(angle=90,hjust=1)) + labs( y= "Percentage of 16S rRNA sequences", x="", title= "Porites") + scale_fill_manual(values=P21)  + theme( legend.key.size = unit(0.2, "cm"),legend.key.width = unit(0.2,"cm")) + facet_grid(Symbiosis~Time, space = "free")
# 
# #png(filename = "Transplants_barplots_ACR.png", res = 300, width = 2500, height = 1200, units = "px")
# #svg(filename = "Transplants_barplots_ACR.svg",  width = 10, height = 5, pointsize = 12) 
# 
# acro=ggplot() +geom_bar(aes(y = Abundance, x = Replicate, fill = Family), data = all.l.acr, stat="identity", position = "fill") +  scale_y_continuous(labels = percent_format(), expand = c(0, 0)) + theme(axis.text.x=element_text(angle=90,hjust=1)) + labs( y= "Percentage of 16S rRNA sequences", x="", title= "Acropora") + scale_fill_manual(values=P21)  + theme( legend.key.size = unit(0.2, "cm"),legend.key.width = unit(0.2,"cm")) + facet_grid(Symbiosis~Time, space = "free")
# 
# aip=ggplot() +geom_bar(aes(y = Abundance, x = Replicate, fill = Family), data = all.l.aip, stat="identity", position = "fill") +  scale_y_continuous(labels = percent_format(), expand = c(0, 0)) + theme(axis.text.x=element_text(angle=90,hjust=1)) + labs( y= "Percentage of 16S rRNA sequences", x="", title= "Aiptasio") + scale_fill_manual(values=P21)  + theme( legend.key.size = unit(0.2, "cm"),legend.key.width = unit(0.2,"cm")) + facet_grid(Symbiosis~Time, space = "free")
# 
# 
# pdf("transplants_family.pdf", onefile = TRUE, width=15,height=15)
# grid.arrange(por, acro, aip)
# dev.off()

