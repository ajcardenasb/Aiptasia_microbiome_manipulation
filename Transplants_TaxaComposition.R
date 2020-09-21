library(reshape2)
library(ggplot2)
library(scales)
library(gridExtra)
library(stringr)
library(ggpubr)

setwd("~/Documents/GitHub/Coral_microbiome_transplants/50_samples_only/Input/")

#############################################################
#####Taxonomic profiles of the 20 most abundant families#####
#############################################################

otu = read.table("OTU_noContaminants.raw.txt", header = TRUE,  row.names=1)
tax = read.table("Tax_noContaminants.txt", header = TRUE,  row.names=1)

#########################################################
#####Taxonomic profiles of the 20 most abundant OTUs#####
#########################################################
otu.n=apply(otu[,1:ncol(otu)], 2, as.numeric)
names=rownames(otu)
rownames(otu.n)=names
otu.n=transform(otu.n, OTU = rownames(otu))
names(otu.n)
topOTUs=otu.n[order(rowSums(otu.n[, c(1:50)]),decreasing = TRUE),][1:20,51] 
otu.top=subset(otu.n, rownames(otu.n) %in% topOTUs) 
otu.bot=subset(otu.n, !rownames(otu.n) %in% topOTUs) 
otu.bot$OTU=gsub(".*","Others", otu.bot$OTU)
names(otu.bot)
others=aggregate(otu.bot[, c(1:50)], by = list(OTU=otu.bot[, 51]), FUN =  sum)
names(others)
others.2=others[, c(2:51,1)]
all.2 =rbind(otu.top, others)
all.l=melt(all.2, id.vars=c("OTU"), variable.name = "OTU", value.name = "Abundance")
colnames(all.l)=c("OTU","Sample","Abundance")

## Add sample information
tax = read.table("Tax_noContaminants.txt", header = TRUE,  row.names=1)
meta = read.table("meta_50.txt", header = TRUE)
all.l$Symbiosis=meta$Symbiosis[match(all.l$Sample, meta$Sample)]
all.l$plot1=meta$final_names[match(all.l$Sample, meta$Sample)]
all.l$plot2=meta$barplots[match(all.l$Sample, meta$Sample)]
all.l$Replicate=meta$Replicate2[match(all.l$Sample, meta$Sample)]
all.l$genus=tax$Genus[match(all.l$OTU, rownames(tax))]
all.l$taxa=paste(all.l$OTU,all.l$genus, sep = " - ")
all.l$plot1=factor(all.l$plot1, levels = c("APO", "APO_AB", "APO+APOinoc", "APO+ACRinoc", "APO+PORinoc", "SYM", "SYM_AB", "SYM+SYMinoc", "SYM+ACRinoc", "SYM+PORinoc"))

apo.plot=subset(all.l, all.l$Symbiosis == "APO")
sym.plot=subset(all.l,  all.l$Symbiosis == "SYM")


P21=c( "#C0C0C0","#771155", "#AA4488", "#CC99BB", "#114477", "#4477AA", "#77AADD", "#117777", "#44AAAA", "#77CCCC", "#117744", "#44AA77", "#88CCAA", "#777711", "#AAAA44", "#DDDD77", "#774411", "#AA7744", "#DDAA77", "#771122", "#AA4455")
#bar.apo=ggplot() +geom_bar(aes(y = Abundance, x = Replicate, fill = taxa), data = apo.plot, stat="identity", position = "fill") +  scale_y_continuous(labels = percent_format(), expand = c(0, 0)) + theme(axis.text.x=element_text(angle=90,hjust=1)) + labs( y= "Percentage of 16S rRNA sequences", x="", title= "") + scale_fill_manual(values=P21)  + theme( legend.key.size = unit(0.2, "cm"),legend.key.width = unit(0.2,"cm"), strip.background=element_rect(fill="white"), legend.position="none", text = element_text(size=12)) + theme_bw() + facet_grid(~plot1 )
#bar.sym=ggplot() +geom_bar(aes(y = Abundance, x = Replicate, fill = taxa), data = sym.plot, stat="identity", position = "fill") +  scale_y_continuous(labels = percent_format(), expand = c(0, 0)) + theme(axis.text.x=element_text(angle=90,hjust=1)) + labs( y= "Percentage of 16S rRNA sequences", x="", title= "") + scale_fill_manual(values=P21) + theme_bw() +  theme(strip.background =element_rect(fill="white"), text = element_text(size=12), legend.position="bottom", legend.key.size = unit(0.2, "cm"),legend.key.width = unit(0.2,"cm")) + facet_grid(~plot1 )

bar.apo=ggplot() +geom_bar(aes(y = Abundance, x = Replicate, fill = taxa), data = apo.plot, stat="identity", position = "fill") +  scale_y_continuous(labels = percent_format(), expand = c(0, 0)) + theme(axis.text.x=element_text(angle=90,hjust=1)) + labs( y= "Percentage of 16S rRNA sequences", x="", title= "") + scale_fill_manual(values=P21)  + theme( legend.key.size = unit(0.2, "cm"),legend.key.width = unit(0.2,"cm")) + facet_grid(~plot1, space = "free")+ theme_bw() +  theme(strip.background =element_rect(fill="white"), text = element_text(size=18),  legend.position="none") 
bar.sym=ggplot() +geom_bar(aes(y = Abundance, x = Replicate, fill = taxa), data = sym.plot, stat="identity", position = "fill") +  scale_y_continuous(labels = percent_format(), expand = c(0, 0)) + theme(axis.text.x=element_text(angle=90,hjust=1)) + labs( y= "Percentage of 16S rRNA sequences", x="", title= "") + scale_fill_manual(values=P21)  + theme( legend.key.size = unit(0.2, "cm"),legend.key.width = unit(0.2,"cm")) + facet_grid(~plot1, space = "free")+ theme_bw() +  theme(strip.background =element_rect(fill="white"), text = element_text(size=18),  legend.position="none") 
names=ggplot() +geom_bar(aes(y = Abundance, x = Replicate, fill = taxa), data = sym.plot, stat="identity", position = "fill") +  scale_y_continuous(labels = percent_format(), expand = c(0, 0)) + theme(axis.text.x=element_text(angle=90,hjust=1)) + labs( y= "Percentage of 16S rRNA sequences", x="", title= "") + scale_fill_manual(values=P21)  + theme( legend.key.size = unit(0.2, "cm"),legend.key.width = unit(0.2,"cm")) + facet_grid(~plot1, space = "free")+ theme_bw() +  theme(strip.background =element_rect(fill="white"), text = element_text(size=18),  legend.position="bottom", legend.title = element_blank())
leg= get_legend(names)
as_ggplot(leg)

pdf("../outputs/barplots_50_otu.pdf", onefile = TRUE, width=20,height=15)
grid.arrange(bar.apo, bar.sym,leg,  ncol = 1)
dev.off()

