library("phyloseq")
library(gridExtra)
library(plyr)

############################################################
##################### phyloseq #############################
############################################################

otu = read.table("OTU_noContaminants.raw.txt", header = TRUE,  row.names=1)
tax = read.table("Tax_noContaminants.txt", header = TRUE,  row.names=1)
met = read.table("transplants_metadata.txt", header = TRUE,  row.names=1)

tax.1=tax[,-1]
otu.t= otu_table(otu, taxa_are_rows=TRUE)
sam.t= sample_data(met)
tax.t= tax_table(as.matrix(tax.1))

phy.all= phyloseq(otu.t, tax.t,  sam.t)

SYM=subset_samples(phy.all, !Symbiosis == "APO") 
APO=subset_samples(phy.all, !Symbiosis == "SYM") 

#COL7= c( "#E7869E", "#DA4167", "#8B2A42","#aed35d", "#3C121D","#5e8408","#0254df") #AB, Aiptasia, Control, CB, Inoc(aipt), InocB, water)
#water #0254df green
# Aiptaisia #E7869E, #8B2A42 #DA4167 #3C121D lighter to darker
#coral #aed35d, #5e8408 lighter to darker
COL3.p=c("#FF3C38", "#D8A038", "#3E7CB1") #Aiptasia, Control, Porites
COL3.a=c("#5B772C", "#FF3C38", "#D8A038") #Acropora, Aiptasia, Control

por.list=c("C", "Aip", "CB", "Inoc", "InocB", "AB") #water
acr.list=c("C", "Aip", "CA", "Inoc", "InocA", "AB") #water
SYM.por=subset_samples(SYM, Condition %in% por.list)
APO.por=subset_samples(APO, Condition %in% por.list)
SYM.acr=subset_samples(SYM, Condition %in% acr.list)
APO.acr=subset_samples(APO, Condition %in% acr.list)


SYM.por.PCOA_br = ordinate(SYM.por, method = "PCoA", distance = "bray")
a=plot_ordination(SYM.por, SYM.por.PCOA_br, color = "Inoculum", shape = "category")  + geom_point(size = 3, alpha = 1) + theme_bw()  + ggtitle("Porites SYM") + theme(plot.title = element_text(hjust = 0.5)) + scale_colour_manual(values=COL3.p)

APO.por.PCOA_br = ordinate(APO.por, method = "PCoA", distance = "bray")
b=plot_ordination(APO.por, APO.por.PCOA_br, color = "Inoculum", shape = "category")  + geom_point(size = 3, alpha = 1) + theme_bw()  + ggtitle("Porites APO") + theme(plot.title = element_text(hjust = 0.5)) + scale_colour_manual(values=COL3.p)

SYM.acr.PCOA_br = ordinate(SYM.acr, method = "PCoA", distance = "bray")
c=plot_ordination(SYM.acr, SYM.acr.PCOA_br, color = "Inoculum", shape = "category")  + geom_point(size = 3, alpha = 1) + theme_bw()  + ggtitle("Acropora SYM") + theme(plot.title = element_text(hjust = 0.5)) + scale_colour_manual(values=COL3.a)

APO.acr.PCOA_br = ordinate(APO.acr, method = "PCoA", distance = "bray")
d=plot_ordination(APO.acr, APO.acr.PCOA_br, color = "Inoculum", shape = "category")  + geom_point(size = 3, alpha = 1) + theme_bw()  + ggtitle("Acropora APO") + theme(plot.title = element_text(hjust = 0.5)) + scale_colour_manual(values=COL3.a)

pdf("transplants_PCoAs.pdf", onefile = TRUE, width=15,height=15)
grid.arrange(a, b, c, d)
dev.off()

### testing id APO and SYMs are similar using the same inoculum
por.list.all=c( "CB",  "C", "Aip") #water
acro.list.all=c( "CA", "C", "Aip") #water
all.por=subset_samples(phy.all, Condition %in% por.list.all)
all.acr=subset_samples(phy.all, Condition %in% acro.list.all)
all.por.PCOA_br = ordinate(all.por, method = "PCoA", distance = "bray")
all.acr.PCOA_br = ordinate(all.acr, method = "PCoA", distance = "bray")

COL2=c("#66717E", "#77C627")

e=plot_ordination(all.por, all.por.PCOA_br, color = "new", shape = "Symbiosis")  + geom_point(size = 3, alpha = 1) + theme_bw()  + ggtitle("Porites") + theme(plot.title = element_text(hjust = 0.5)) #+ scale_colour_manual(values=COL2)
f=plot_ordination(all.acr, all.acr.PCOA_br,  color = "new", shape = "Symbiosis")  + geom_point(size = 3, alpha = 1) + theme_bw()  + ggtitle("Acropora") + theme(plot.title = element_text(hjust = 0.5)) #+ scale_colour_manual(values=COL2)

pdf("transplants_PCoAs_2.pdf", onefile = TRUE, width=15,height=15)
grid.arrange(e,f)
dev.off()

############################################################
##################### Alpha-diversity ######################
############################################################

alpha.SYM.por=estimate_richness(SYM.por, split = TRUE, measures = c("Observed", "Chao1"))
alpha.SYM.por$Inoculum=met$Inoculum[match(rownames(alpha.SYM.por), rownames(met))]
alpha.SYM.por$category=met$category[match(rownames(alpha.SYM.por), rownames(met))]
alpha.sym.por.sum=ddply(alpha.SYM.por, c("Inoculum", "category"), summarise, N = length(Chao1), mean = mean(Chao1), sd= sd(Chao1), se= sd / sqrt(N))
alpha.sym.por.sum$lower=alpha.sym.por.sum$mean-alpha.sym.por.sum$sd
alpha.sym.por.sum$upper=alpha.sym.por.sum$mean+alpha.sym.por.sum$sd
alpha.sym.por.sum[is.na(alpha.sym.por.sum)] <- 0
alpha.sym.por.sum$category=factor(alpha.sym.por.sum$category, levels = c("Tank","AB", "Inoc", "3D", "7D"))
e=ggplot(alpha.sym.por.sum, aes(x=category, weight=mean, ymin=lower, ymax=upper, fill=Inoculum)) + geom_bar(position=position_dodge(), aes(y=mean), stat="identity", alpha=0.7) + geom_errorbar(position=position_dodge(width=0.9), colour="black") + scale_fill_manual(values=COL3.p) + labs( y= "Chao1", x="", title= "Porites SYM")

alpha.APO.por=estimate_richness(APO.por, split = TRUE, measures = c("Observed", "Chao1"))
alpha.APO.por$Inoculum=met$Inoculum[match(rownames(alpha.APO.por), rownames(met))]
alpha.APO.por$category=met$category[match(rownames(alpha.APO.por), rownames(met))]
alpha.APO.por.sum=ddply(alpha.APO.por, c("Inoculum", "category"), summarise, N = length(Chao1), mean = mean(Chao1), sd= sd(Chao1), se= sd / sqrt(N))
alpha.APO.por.sum$lower=alpha.APO.por.sum$mean-alpha.APO.por.sum$sd
alpha.APO.por.sum$upper=alpha.APO.por.sum$mean+alpha.APO.por.sum$sd
alpha.APO.por.sum[is.na(alpha.APO.por.sum)] <- 0
alpha.APO.por.sum$category=factor(alpha.APO.por.sum$category, levels = c("Tank","AB", "Inoc", "3D", "7D"))
f=ggplot(alpha.APO.por.sum, aes(x=category, weight=mean, ymin=lower, ymax=upper, fill=Inoculum)) + geom_bar(position=position_dodge(), aes(y=mean), stat="identity", alpha=0.7) + geom_errorbar(position=position_dodge(width=0.9), colour="black") + scale_fill_manual(values=COL3.p) + labs( y= "Chao1", x="", title= "Porites APO")

alpha.SYM.acr=estimate_richness(SYM.acr, split = TRUE, measures = c("Observed", "Chao1"))
alpha.SYM.acr$Inoculum=met$Inoculum[match(rownames(alpha.SYM.acr), rownames(met))]
alpha.SYM.acr$category=met$category[match(rownames(alpha.SYM.acr), rownames(met))]
alpha.sym.acr.sum=ddply(alpha.SYM.acr, c("Inoculum", "category"), summarise, N = length(Chao1), mean = mean(Chao1), sd= sd(Chao1), se= sd / sqrt(N))
alpha.sym.acr.sum$lower=alpha.sym.acr.sum$mean-alpha.sym.acr.sum$sd
alpha.sym.acr.sum$upper=alpha.sym.acr.sum$mean+alpha.sym.acr.sum$sd
alpha.sym.acr.sum[is.na(alpha.sym.acr.sum)] <- 0
alpha.sym.acr.sum$category=factor(alpha.sym.acr.sum$category, levels = c("Tank","AB", "Inoc", "3D", "7D"))
g=ggplot(alpha.sym.acr.sum, aes(x=category, weight=mean, ymin=lower, ymax=upper, fill=Inoculum)) + geom_bar(position=position_dodge(), aes(y=mean), stat="identity", alpha=0.7) + geom_errorbar(position=position_dodge(width=0.9), colour="black") + scale_fill_manual(values=COL3.a) + labs( y= "Chao1", x="", title= "Acropora SYM")

alpha.APO.acr=estimate_richness(APO.acr, split = TRUE, measures = c("Observed", "Chao1"))
alpha.APO.acr$Inoculum=met$Inoculum[match(rownames(alpha.APO.acr), rownames(met))]
alpha.APO.acr$category=met$category[match(rownames(alpha.APO.acr), rownames(met))]
alpha.APO.acr.sum=ddply(alpha.APO.acr, c("Inoculum", "category"), summarise, N = length(Chao1), mean = mean(Chao1), sd= sd(Chao1), se= sd / sqrt(N))
alpha.APO.acr.sum$lower=alpha.APO.acr.sum$mean-alpha.APO.acr.sum$sd
alpha.APO.acr.sum$upper=alpha.APO.acr.sum$mean+alpha.APO.acr.sum$sd
alpha.APO.acr.sum[is.na(alpha.APO.acr.sum)] <- 0
alpha.APO.acr.sum$category=factor(alpha.APO.acr.sum$category, levels = c("Tank","AB", "Inoc", "3D", "7D"))
h=ggplot(alpha.APO.acr.sum, aes(x=category, weight=mean, ymin=lower, ymax=upper, fill=Inoculum)) + geom_bar(position=position_dodge(), aes(y=mean), stat="identity", alpha=0.7) + geom_errorbar(position=position_dodge(width=0.9), colour="black") + scale_fill_manual(values=COL3.a) + labs( y= "Chao1", x="", title= "Acropora APO")

pdf("transplants_Chao1.pdf", onefile = TRUE, width=15,height=15)
grid.arrange(e, f, g, h)
dev.off()

############################################################
##################### Venn Diagrams ########################
############################################################
library(gplots)
grp1=as.list(unique(grp1$Group))
venn(list(grp3, grp1))
write.table(grp1, "grp1_list.txt", quote = FALSE, row.names = FALSE)