#Overview of this script:
    #Computes the changes in physicochemical properties for all single mutants of deep mutational scan
    #Relates experimentally determined DMS functional scores to changes in physicochemical properties
    #Generates violin plots to identify global trends in host-specific enrichment of single mutants

#Read in the file that contains table of physicochemical properties of each amino acid 
pchem    <- read.csv("../reference/AA_physicochemical_properties.csv", header=T, stringsAsFactors=F)
#Read in the file that contains the table of experimentally determined DMS functional scores
fscore   <- read.csv("../reference/DMS_Fscores.csv", header=T, stringsAsFactors=F)

#Compute the change in each physicochemical property for all single mutants 
#----------------------------------------------------------------------------------------------
for(i in 1:length(fscore[,1]))
{  
  wt  <- subset(pchem, pchem[,1]==substr(fscore[i,1],1,1))
  mut <- subset(pchem, pchem[,1]==substr(fscore[i,3],1,1))
  
  for(property in 1:9)
  {
    fscore[i,property+8] <- mut[1,property+1] - wt[1,property+1]
  }
}
for(property in 1:9)
{
  colnames(fscore)[property+8] <- paste0("delta_",colnames(pchem)[property+1])
}
write.csv(fscore, "../analysis/delta_physicochemical_properties.csv", row.names=F, quote=F)
#----------------------------------------------------------------------------------------------

#Based on Fscores, assign a phenotype ("depleted", "Tolerated", or "Enriched") to each single mutant 
#----------------------------------------------------------------------------------------------
pheno <- as.data.frame(matrix(0,1,11))
colnames(pheno)[1:9] <- colnames(fscore)[9:17]
colnames(pheno)[10:11] <- c("strain","group")

c <- 0
for(strain in 1:5)
{
  for(i in 1:length(fscore[,1]))
  {
    c <- c + 1
    pheno[c,1:9] <- fscore[i,9:17]
    pheno[c,10] <- colnames(fscore)[3+strain]
    if(is.na(fscore[i,3+strain]))
    {
      pheno[c,11] <- NA
    } 
    else if(fscore[i,3+strain] < 0.1)
    {
      pheno[c,11] <- "Depleted"
    } 
    else if(fscore[i,3+strain] < 10)
    {
      pheno[c,11] <- "Tolerated"
    } 
    else
    {
      pheno[c,11] <- "Enriched"
    }
  }
}
#----------------------------------------------------------------------------------------------

#Generate violin plots to identify global trends in host-specific enrichment of single mutants
#----------------------------------------------------------------------------------------------
library(ggplot2)
library(ggpubr)
my_comparisons <- list(c("Depleted", "Tolerated"), c("Tolerated", "Enriched"), c("Depleted", "Enriched"))
pheno <- subset(pheno, !is.na(pheno[,11]))
pdf(file="../analysis/pdfs/violin_delta_hydrophobicity.pdf", width=12, height=4.5, useDingbats=F)
p <- ggviolin(pheno, x="group", y="delta_hydrophobicity", facet.by="strain", short.panel.labs=T, 
              add.params = list(fill = "white"), fill="group", 
              palette=c("#c2a5cf","gray90","#a6dba0"), add="boxplot", repel=T,
              ylab=expression(paste(Delta,Hydrophobicity)), xlab="", ylim=c(-5,7.5),
              order=c("Depleted","Tolerated","Enriched"))
p + stat_compare_means(method="kruskal.test", label="p.format", label.x=1.6, label.y=7.3) + stat_compare_means(comparisons=my_comparisons) + theme_classic() + facet_grid(. ~ strain) + theme(legend.position="none", text=element_text(size=20), axis.text.x=element_text(angle=45, colour="black", vjust=0.5), axis.text.y=element_text(colour="black"))
dev.off()

pdf(file="../analysis/pdfs/violin_delta_hydrophilicity.pdf", width=12, height=4.5, useDingbats=F)
p <- ggviolin(pheno, x="group", y="delta_hydrophilicity", facet.by="strain", short.panel.labs=T, 
              add.params = list(fill = "white"), fill="group", 
              palette=c("#c2a5cf","gray90","#a6dba0"), add="boxplot", repel=T,
              ylab=expression(paste(Delta,Hydrophilicity)), xlab="", ylim=c(-10,13),
              order=c("Depleted","Tolerated","Enriched"))
p + stat_compare_means(method="kruskal.test", label="p.format", label.x=1.6, label.y=12.5) + stat_compare_means(comparisons=my_comparisons) + theme_classic() + facet_grid(. ~ strain) + theme(legend.position="none", text=element_text(size=20), axis.text.x=element_text(angle=45, colour="black", vjust=0.5), axis.text.y=element_text(colour="black"))
dev.off()

pdf(file="../analysis/pdfs/violin_delta_polarity.pdf", width=12, height=4.5, useDingbats=F)
p <- ggviolin(pheno, x="group", y="delta_polarity", facet.by="strain", short.panel.labs=T, 
              add.params = list(fill = "white"), fill="group", 
              palette=c("#c2a5cf","gray90","#a6dba0"), add="boxplot", repel=T,
              ylab=expression(paste(Delta,Polarity)), xlab="", ylim=c(-12,16),
              order=c("Depleted","Tolerated","Enriched"))
p + stat_compare_means(method="kruskal.test", label="p.format", label.x=1.6, label.y=15.5) + stat_compare_means(comparisons=my_comparisons) + theme_classic() + facet_grid(. ~ strain) + theme(legend.position="none", text=element_text(size=20), axis.text.x=element_text(angle=45, colour="black", vjust=0.5), axis.text.y=element_text(colour="black"))
dev.off()

pdf(file="../analysis/pdfs/violin_delta_polarizability.pdf", width=12, height=4.5, useDingbats=F)
p <- ggviolin(pheno, x="group", y="delta_polarizability", facet.by="strain", short.panel.labs=T, 
              add.params = list(fill = "white"), fill="group", 
              palette=c("#c2a5cf","gray90","#a6dba0"), add="boxplot", repel=T,
              ylab=expression(paste(Delta,Polarizability)), xlab="", ylim=c(-0.5,0.8),
              order=c("Depleted","Tolerated","Enriched"))
p + stat_compare_means(method="kruskal.test", label="p.format", label.x=1.6, label.y=0.78) + stat_compare_means(comparisons=my_comparisons) + theme_classic() + facet_grid(. ~ strain) + theme(legend.position="none", text=element_text(size=20), axis.text.x=element_text(angle=45, colour="black", vjust=0.5), axis.text.y=element_text(colour="black"))
dev.off()

pdf(file="../analysis/pdfs/violin_delta_mass.pdf", width=12, height=4.5, useDingbats=F)
p <- ggviolin(pheno, x="group", y="delta_mass", facet.by="strain", short.panel.labs=T, 
              add.params = list(fill = "white"), fill="group", 
              palette=c("#c2a5cf","gray90","#a6dba0"), add="boxplot", repel=T,
              ylab=expression(paste(Delta,Mass (g/mol))), xlab="", ylim=c(-180,250),
              order=c("Depleted","Tolerated","Enriched"))
p + stat_compare_means(method="kruskal.test", label="p.format", label.x=1.6, label.y=242) + stat_compare_means(comparisons=my_comparisons) + theme_classic() + facet_grid(. ~ strain) + theme(legend.position="none", text=element_text(size=20), axis.text.x=element_text(angle=45, colour="black", vjust=0.5), axis.text.y=element_text(colour="black"))
dev.off()
#----------------------------------------------------------------------------------------------

