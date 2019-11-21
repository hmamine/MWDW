####################################################################
#R scripts written by M. Amine Hassani - hassani.medamine@gmail.com#
####################################################################
#required packages
pkg=c("plyr", "ggplot2", "phyloseq", "data.table", "reshape2", "picante",  "PMCMR")
lapply(pkg, library, character.only = TRUE)
rm(list = ls()[grep("*.*", ls())])
alphaInd = c("Shannon", "Observed")
color_palette<-c("#cc2a36","#166590","#11902d","#000000")
p<-list()
theme_new <- theme (
	panel.grid.major = element_blank(),
	panel.grid.minor = element_blank(),
	legend.position="none",
	axis.text.x = element_text(angle=90, vjust=1)
	)
sd=read.table("sample_data_merge.txt", sep="\t", row.names=1, header=1)
SD=sample_data(sd)
#####################################
#uploading and preparing phyloseq ###
#object for the bacterial microbiota#
#####################################
	mat.B=read.table( "otu_table_merge_bacteria.txt", sep="\t", row.names=1, header=T)
	mat.B=as.matrix(mat.B)
	OTU.B=otu_table(mat.B, taxa_are_rows=T) 
	tax.B=read.table("taxonomy_bacteria.txt", sep="\t", row.names=1, header=1)
	tax.B=as.matrix(tax.B)
	ColPhylum<-c("Acidobacteria","Actinobacteria","Bacteroidetes","Firmicutes","Proteobacteria")
	tax.B[,"Phylum"] <- ifelse(!tax.B[,"Phylum"] %in% ColPhylum, "Other", ifelse(tax.B[,"Phylum"] == "Bacteroidetes", "Bacteroidetes", ifelse(tax.B[,"Phylum"] == "Firmicutes", "Firmicutes", ifelse (tax.B[,"Phylum"] == "Actinobacteria", "Actinobacteria", ifelse (tax.B[,"Phylum"] == "Acidobacteria", "Acidobacteria", "Proteobacteria")))))
	TAXA.B=tax_table(tax.B)
	TREE=read.tree("tree.nwk")
	physeq.B=phyloseq(OTU.B, TAXA.B, SD, TREE) 
#rarefication of reads
	rf=min(sample_sums(physeq.B))
	physeq.rrf=rarefy_even_depth(physeq.B, rf, replace=TRUE, rngseed = 131)
SUBSET<-list("Acidobacteria","Actinobacteria","Bacteroidetes","Firmicutes","Proteobacteria", "Other")
for ( i in SUBSET) 
{
	print("########################################")
	print(i)
	physeq.rrf.sub= subset_taxa(physeq.rrf, Phylum %in% i )
	#Ploting species richness	
	p_nat=plot_richness(physeq.rrf.sub,"Concat2","Geno" ,measures=alphaInd)
	p_nat$layers <- p_nat$layers[-1]
	p1=p_nat+geom_boxplot(data=p_nat$data, aes(x=Concat2, y=value, color=Geno))+
	geom_point(size=4)+theme_bw()+
	theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
	scale_colour_manual(values=color_palette)	
	subp1 =	ggplot(data=p_nat$data[p_nat$data$variable=="Observed",], aes(x=Concat2, y=value, color=Geno))+
	geom_boxplot(width=1)+
	geom_point(size=1, color="black")+theme_bw()+ theme_new +
	scale_colour_manual(values=color_palette) + 
	facet_wrap(~variable) +
	ylab("Richness")
	#testing for significance
	print("########################################")
	print("leaf samples")
	DT.R.L<-subp1$data[subp1$data$Comp2 == "Leaf",]
	print(kruskal.test(data=DT.R.L, value ~ Concat2))
	print(posthoc.kruskal.conover.test(data=DT.R.L, value ~ Concat2, method="BH"))
	print("########################################")
	print("root samples")
	DT.R.R<-subp1$data[subp1$data$Comp2 == "Root",]
	print(kruskal.test(data=DT.R.R, value ~ Concat2))
	print(posthoc.kruskal.conover.test(data=DT.R.R, value ~ Concat2, method="BH"))
	print("########################################")
	print("rhizo. samples")
	DT.R.Rz<-subp1$data[subp1$data$Comp == "Rhizo",]
	print(kruskal.test(data=DT.R.Rz, value ~ Concat2))
	print(posthoc.kruskal.conover.test(data=DT.R.Rz, value ~ Concat2, method="BH"))
	p[[i]]<-subp1
}
	#pdf("Supplementary_Fig_1.pdf",paper="A4" ,useDingbats=FALSE)
	gridExtra::grid.arrange(p$Acidobacteria, p$Actinobacteria, p$Bacteroidetes, p$Firmicutes, p$Proteobacteria, p$Other, ncol=2)
	#dev.off()
###################################################
### figures were edited to increase readability ###
###################################################
           
            
          
 









