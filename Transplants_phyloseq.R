library("phyloseq")
library(gridExtra)
library(plyr)

############################################################
##################### phyloseq #############################
############################################################

otu = read.table("otu_noConta_noOutliers.txt", header = TRUE,  row.names=1)
tax = read.table("Tax_noContaminants.txt", header = TRUE,  row.names=1)
met = read.table("transplants_metadata_noOutliers.txt", header = TRUE,  row.names=1)
met$barplots=factor(met$barplots, levels = c("Tank", "AB", "Aiptasia", "Acropora", "Porites"))

tax.1=tax[,-1]
otu.t= otu_table(otu, taxa_are_rows=TRUE)
sam.t= sample_data(met)
tax.t= tax_table(as.matrix(tax.1))

phy.all= phyloseq(otu.t, tax.t,  sam.t)


COL=c("#A3C38E", "#C1A222",  "#577545",  "#39597C", "#7E9EC1") #cont, contInoc, Acro, Por, AB
#AB  "#C1A222"
#ACROinoc  "#39597C"
#PORinoc "#7E9EC1"
#AIPinoc "#577545"
#Tank controls "#A3C38E"

#list=c("C","Aip", "CB",  "AB", "CA") #water

APO=subset_samples(phy.all, !Symbiosis == "SYM") 
APO.2=subset_samples(APO,  !barplots == "ignore")
APO.PCOA_br = ordinate(APO.2, method = "PCoA", distance = "bray")
a=plot_ordination(APO.2, APO.PCOA_br, color = "barplots", )  + geom_point(size = 3, alpha = 1) + theme_bw()  + ggtitle("APO") + theme(plot.title = element_text(hjust = 0.5)) + scale_colour_manual(values=COL)+ theme(legend.title = element_blank())


SYM=subset_samples(phy.all, !Symbiosis == "APO") 
SYM.2=subset_samples(SYM, !barplots == "ignore")
SYM.PCOA_br = ordinate(SYM.2, method = "PCoA", distance = "bray")
b=plot_ordination(SYM, SYM.PCOA_br, color = "barplots", )  + geom_point(size = 3, alpha = 1) + theme_bw()  + ggtitle("SYM") + theme(plot.title = element_text(hjust = 0.5)) + scale_colour_manual(values=COL)+ theme(legend.title = element_blank())


pdf("transplants_PCoAs.pdf", onefile = TRUE, width=10,height=5)
pcoas=grid.arrange(a, b, ncol = 2)
dev.off()


############################################################
##################### Alpha-diversity ######################
############################################################

alpha.SYM=estimate_richness(SYM.2, split = TRUE, measures = c("Observed", "Chao1"))
alpha.SYM$barplots=met$barplots[match(rownames(alpha.SYM), rownames(met))]
alpha.sym.sum=ddply(alpha.SYM, c("barplots"), summarise, N = length(Chao1), mean = mean(Chao1), sd= sd(Chao1), se= sd / sqrt(N))
alpha.sym.sum$lower=alpha.sym.sum$mean-alpha.sym.sum$sd
alpha.sym.sum$upper=alpha.sym.sum$mean+alpha.sym.sum$sd
alpha.sym.sum[is.na(alpha.sym.sum)] <- 0
alpha.sym.sum$barplots=factor(alpha.sym.sum$barplots, levels =  c("Tank", "AB", "Aiptasia", "Acropora", "Porites"))
e=ggplot(alpha.sym.sum, aes(x=barplots, weight=mean, ymin=lower, ymax=upper, fill=barplots)) + geom_bar(position=position_dodge(), aes(y=mean), stat="identity", alpha=0.7) + geom_errorbar(position=position_dodge(width=0.5), colour="black") + scale_fill_manual(values=COL) + labs( y= "Chao1", x="", title= "SYM") + theme(legend.title = element_blank())

alpha.apo=estimate_richness(APO.2, split = TRUE, measures = c("Observed", "Chao1"))
alpha.apo$barplots=met$barplots[match(rownames(alpha.apo), rownames(met))]
alpha.apo.sum=ddply(alpha.apo, c("barplots"), summarise, N = length(Chao1), mean = mean(Chao1), sd= sd(Chao1), se= sd / sqrt(N))
alpha.apo.sum$lower=alpha.apo.sum$mean-alpha.apo.sum$sd
alpha.apo.sum$upper=alpha.apo.sum$mean+alpha.apo.sum$sd
alpha.apo.sum[is.na(alpha.apo.sum)] <- 0
alpha.apo.sum$barplots=factor(alpha.apo.sum$barplots, levels =  c("Tank", "AB", "Aiptasia", "Acropora", "Porites"))
f=ggplot(alpha.apo.sum, aes(x=barplots, weight=mean, ymin=lower, ymax=upper, fill=barplots)) + geom_bar(position=position_dodge(), aes(y=mean), stat="identity", alpha=0.7) + geom_errorbar(position=position_dodge(width=0.5), colour="black") + scale_fill_manual(values=COL) + labs( y= "Chao1", x="", title= "APO")+ theme(legend.title = element_blank())

pdf("transplants_Chao1.pdf", onefile = TRUE, width=12,height=5)
grid.arrange(e, f, ncol = 2)
dev.off()


# alpha.SYM.acr=estimate_richness(SYM.acr, split = TRUE, measures = c("Observed", "Chao1"))
# alpha.SYM.acr$Inoculum=met$Inoculum[match(rownames(alpha.SYM.acr), rownames(met))]
# alpha.SYM.acr$category=met$category[match(rownames(alpha.SYM.acr), rownames(met))]
# alpha.sym.acr.sum=ddply(alpha.SYM.acr, c("Inoculum", "category"), summarise, N = length(Chao1), mean = mean(Chao1), sd= sd(Chao1), se= sd / sqrt(N))
# alpha.sym.acr.sum$lower=alpha.sym.acr.sum$mean-alpha.sym.acr.sum$sd
# alpha.sym.acr.sum$upper=alpha.sym.acr.sum$mean+alpha.sym.acr.sum$sd
# alpha.sym.acr.sum[is.na(alpha.sym.acr.sum)] <- 0
# alpha.sym.acr.sum$category=factor(alpha.sym.acr.sum$category, levels = c("Tank","AB", "Inoc", "3D", "7D"))
# g=ggplot(alpha.sym.acr.sum, aes(x=category, weight=mean, ymin=lower, ymax=upper, fill=Inoculum)) + geom_bar(position=position_dodge(), aes(y=mean), stat="identity", alpha=0.7) + geom_errorbar(position=position_dodge(width=0.9), colour="black") + scale_fill_manual(values=COL3.a) + labs( y= "Chao1", x="", title= "Acropora SYM")
# 
# alpha.APO.acr=estimate_richness(APO.acr, split = TRUE, measures = c("Observed", "Chao1"))
# alpha.APO.acr$Inoculum=met$Inoculum[match(rownames(alpha.APO.acr), rownames(met))]
# alpha.APO.acr$category=met$category[match(rownames(alpha.APO.acr), rownames(met))]
# alpha.APO.acr.sum=ddply(alpha.APO.acr, c("Inoculum", "category"), summarise, N = length(Chao1), mean = mean(Chao1), sd= sd(Chao1), se= sd / sqrt(N))
# alpha.APO.acr.sum$lower=alpha.APO.acr.sum$mean-alpha.APO.acr.sum$sd
# alpha.APO.acr.sum$upper=alpha.APO.acr.sum$mean+alpha.APO.acr.sum$sd
# alpha.APO.acr.sum[is.na(alpha.APO.acr.sum)] <- 0
# alpha.APO.acr.sum$category=factor(alpha.APO.acr.sum$category, levels = c("Tank","AB", "Inoc", "3D", "7D"))
# h=ggplot(alpha.APO.acr.sum, aes(x=category, weight=mean, ymin=lower, ymax=upper, fill=Inoculum)) + geom_bar(position=position_dodge(), aes(y=mean), stat="identity", alpha=0.7) + geom_errorbar(position=position_dodge(width=0.9), colour="black") + scale_fill_manual(values=COL3.a) + labs( y= "Chao1", x="", title= "Acropora APO")


