####################################################################
#R scripts written by M. Amine Hassani - hassani.medamine@gmail.com#
####################################################################

#required packages for the analysis the analysis
pkg=c("plyr", "ggplot2", "phyloseq", "data.table", "reshape2", "picante",  "PMCMR")
lapply(pkg, library, character.only = TRUE)

rm(list = ls()[grep("*.*", ls())])

alphaInd = c("Shannon", "Observed")
color_palette<-c("#cc2a36","#166590","#11902d","#000000")

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
	TAXA.B=tax_table(tax.B)
	TREE=read.tree("tree.nwk")

	physeq.B=phyloseq(OTU.B, TAXA.B, SD, TREE) 

#rarefication of reads
	rf=min(sample_sums(physeq.B))
	physeq.rrf=rarefy_even_depth(physeq.B, rf, replace=TRUE, rngseed = 131)

#computing species richness and Shanon index	
	p_nat=plot_richness(physeq.rrf,"Concat2","Geno" ,measures=alphaInd)
	p_nat$layers <- p_nat$layers[-1]
#	p_nat$data$Concat <- ordered(p_nat$data$Concat, levels= neworder)
	p1=p_nat+geom_boxplot(data=p_nat$data, aes(x=Concat2, y=value, color=Geno))+
	geom_point(size=4)+theme_bw()+
	theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
	scale_colour_manual(values=color_palette)

#Preparing plots

### Species richness - Fig. 1a	
subp1 =	ggplot(data=p_nat$data[p_nat$data$variable=="Observed",],aes(x=Concat2, y=value, color=Geno))+
	geom_boxplot(width=1)+
	geom_point(size=1, color="black")+theme_bw()+ theme_new +
	scale_colour_manual(values=color_palette) + 
	facet_wrap(~variable) +
	ylab("Richness")
#testing for significance - richness
DT.Ob.L<-subp1$data[subp1$data$Comp2 == "Leaf",]
#kruskal.test(data=DT.Ob.L, value ~ Concat2)
#posthoc.kruskal.conover.test(data=DT.Ob.L, value ~ Concat2, method="BH")
DT.Ob.R<-subp1$data[subp1$data$Comp2 == "Root",]
#kruskal.test(data=DT.Ob.R, value ~ Concat2)
DT.Ob.Rz<-subp1$data[subp1$data$Comp == "Rhizo",]
#kruskal.test(data=DT.Ob.Rz, value ~ Concat2)
#posthoc.kruskal.conover.test(data=DT.Ob.Rz, value ~ Concat2, method="BH")
panel_a<-subp1

### Shannon index - Fig. 1b	
subp2 =	ggplot(data=p_nat$data[p_nat$data$variable=="Shannon",],aes(x=Concat2, y=value, color=Geno))+
	geom_boxplot(width=1)+
	geom_point(size=1, color="black")+theme_bw()+ theme_new +
	scale_colour_manual(values=color_palette) + 
	facet_wrap(~variable) +
	ylab("Shannon Index")
#testing significance - Evenness
DT.Sh.L<-subp2$data[subp2$data$Comp2 == "Leaf",]
#kruskal.test(data=DT.Sh.L, value ~ Concat2)
DT.Sh.R<-subp2$data[subp2$data$Comp2 == "Root",]
#kruskal.test(data=DT.Sh.R, value ~ Concat2)
DT.Sh.Rz<-subp2$data[subp2$data$Comp == "Rhizo",]
#kruskal.test(data=DT.Sh.Rz, value ~ Concat2)
#posthoc.kruskal.conover.test(data=DT.Sh.Rz, value ~ Concat2, method="BH")
panel_b<-subp2

### Faith's phylogenetic index - Fig. 1c 
	otu.mat.rrf=as( otu_table( physeq.rrf), "matrix" )
	Faith.PD<-data.table(pd(t(otu.mat.rrf), TREE), keep.rownames=T, key="rn")
	setnames(Faith.PD, "rn", "samples")
	DT_PD<-merge(subp2$data[,c(-7:-8)], Faith.PD)
	DT_PD$variable <- "FaithPD"
#Preparing PD index plots
subp3 = ggplot(data=DT_PD, aes(x=Concat2, y=PD, color=Geno))+
	geom_boxplot(width=1)+
	geom_point(size=1, color="black")+theme_bw()+ theme_new +
	scale_colour_manual(values=color_palette) + 
	facet_wrap(~variable) +
	ylab("Phylogenetic diversity")
#testing significance - Faith PD
DT.PD.L<-subp3$data[subp3$data$Comp2 == "Leaf",]
#kruskal.test(data=DT.PD.L, PD ~ Concat2)
DT.PD.R<-subp3$data[subp3$data$Comp2 == "Root",]
#kruskal.test(data=DT.PD.R, PD ~ Concat2)
DT.PD.Rz<-subp3$data[subp3$data$Comp == "Rhizo",]
#kruskal.test(data=DT.PD.Rz, PD ~ Concat2)
panel_c<-subp3

##################################
#uploading and preparing phyloseq#
#object for the fungal microbiota#
##################################
	mat.F=read.table( "otu_table_merge_fungi.txt", sep="\t", row.names=1, header=T)
	mat.F=as.matrix(mat.F)
	OTU.F=otu_table(mat.F, taxa_are_rows=T) 
	tax.F=read.table("taxonomy_fungi.txt", sep="\t", row.names=1, header=1)
	tax.F=as.matrix(tax.F)
	TAXA.F=tax_table(tax.F)
	physeq.F= phyloseq(OTU.F, TAXA.F, SD) 

#Ploting species richness	
	p_nat=plot_richness(physeq.F,"Concat2","Geno" ,measures=alphaInd)
	p_nat$layers <- p_nat$layers[-1]
	p1=p_nat+geom_boxplot(data=p_nat$data, aes(x=Concat2, y=value, color=Geno))+
	geom_point(size=4)+theme_bw()+
	theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
	scale_colour_manual(values=color_palette)

### Fungal species richness - Fig. 1d	
subp4 =	ggplot(data=p_nat$data[p_nat$data$variable=="Observed",],aes(x=Concat2, y=value, color=Geno))+
	geom_boxplot(width=1)+
	geom_point(size=1, color="black")+theme_bw()+ theme_new +
	scale_colour_manual(values=color_palette) + 
	facet_wrap(~variable) +
	ylab("Richness")
#testing significance - Richness
DT.Ob.L<-subp4$data[subp4$data$Comp2 == "Leaf",]
#kruskal.test(data=DT.Ob.L, value ~ Concat2)
DT.Ob.R<-subp4$data[subp4$data$Comp2 == "Root",]
#kruskal.test(data=DT.Ob.R, value ~ Concat2)
#posthoc.kruskal.conover.test(data=DT.Ob.R, value ~ Concat2, method="BH")
DT.Ob.Rz<-subp4$data[subp4$data$Comp == "Rhizo",]
#kruskal.test(data=DT.Ob.Rz, value ~ Concat2)
panel_d<-subp4

### Fungal Shannon index - Fig. 1e	
subp5 =	ggplot(data=p_nat$data[p_nat$data$variable=="Shannon",],aes(x=Concat2, y=value, color=Geno))+
	geom_boxplot(width=1)+
	geom_point(size=1, color="black")+theme_bw()+ theme_new +
	scale_colour_manual(values=color_palette) + 
	facet_wrap(~variable) +
	ylab("Shannon Index") 
#testing significance - Evenness
DT.Sh.L<-subp5$data[subp5$data$Comp2 == "Leaf",]
#kruskal.test(data=DT.Sh.L, value ~ Concat2)
DT.Sh.R<-subp5$data[subp5$data$Comp2 == "Root",]
#kruskal.test(data=DT.Sh.R, value ~ Concat2)
#posthoc.kruskal.conover.test(data=DT.Sh.R, value ~ Concat2, method="BH")
DT.Sh.Rz<-subp5$data[subp5$data$Comp == "Rhizo",]
#kruskal.test(data=DT.Sh.Rz, value ~ Concat2)
panel_e<-subp5

	#pdf("Fig1.pdf", paper="A4", useDingbats=FALSE)
	gridExtra::grid.arrange(panel_a, panel_b, panel_c, panel_d, panel_e, nrow=3, ncol=2)
	#dev.off()

###################################################
### figures were edited to increase readability ###
###################################################
           
            
          
 









