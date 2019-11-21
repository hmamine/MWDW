####################################################################
#R scripts written by M. Amine Hassani - hassani.medamine@gmail.com#
####################################################################
#required packages
pkg=c("plyr", "ggplot2", "phyloseq", "data.table", "reshape2", "GUniFrac", "ape", "phytools", "metagenomeSeq","PMCMR")
lapply(pkg, library, character.only = TRUE)
rm(list = ls()[grep("*.*", ls())])
color_palette.1<-c("#cc2a36","#166590","#11902d","#000000")
color_palette.2<-c("#cc2a36","#166590","#11902d")
shape.1=c(17,19,15,13)
shape.2=c(17)
shape.3=c(15)
shape.4=c(19,13)
theme_new <- theme (
	panel.border = element_rect(),
	panel.grid.major = element_blank(),
	panel.grid.minor = element_blank(),
	legend.position="none"
	)
meth1=("PCoA")
meth2=("CAP")
dist1=("bray")

# uploading data
	mat=read.table( "otu_table_merge_fungi.txt", sep="\t", row.names=1, header=T)
	mat=as.matrix(mat)
	OTU=otu_table(mat, taxa_are_rows=T) 
	tax=read.table("taxonomy_fungi.txt", sep="\t", row.names=1, header=1)
	tax=as.matrix(tax)
	TAXA=tax_table(tax)
	sd=read.table("sample_data_merge.txt", sep="\t", row.names=1, header=1)
	SD=sample_data(sd)
	physeq= phyloseq(OTU, TAXA, SD)
#normalization of count reads using CSS 
	otumat=as( otu_table( physeq ), "matrix" )
	mp=newMRexperiment( ( otumat ) )
	physeq_norm=phyloseq( otu_table( MRcounts( cumNorm( mp, p=cumNormStat( mp ) ),
	norm=TRUE, log=TRUE ),taxa_are_rows=TRUE ), TAXA, SD )
#computing Bray-Curtis distances
	dist_BC=distance( physeq_norm, dist1 )

#computing unconstrained PCoA
	pcoa_BC=ordinate( physeq_norm, meth1 , dist_BC )
	pBC<-plot_ordination( physeq_norm, pcoa_BC, shape="Comp", color="Geno")
	pBC$layers<-pBC$layers[-1]	
	p1=pBC+geom_point(size=3)+theme_bw()+theme_new+
	scale_shape_manual(values=shape.1)+ 
	scale_color_manual(values=color_palette.1)
print("####################################################################")
print("Permanova test for the term habitat (i.e. leaf, root, rhizo., Unpl.)")
DT.all=data.table( as( sample_data( physeq_norm ), "data.frame" ),keep.rownames=T, key="rn" )
print(with( DT.all, adonis ( dist_BC ~ Comp ) ) )
print("####################################################################")

#computing constrained PCoA
## for leaf samples	
	physeq_norm.L=subset_samples( physeq_norm, Comp == "Leaf" )
	dist_BC.L=distance( physeq_norm.L, dist1 )
	cpcoa_BC.L=ordinate( physeq_norm.L, meth2 , dist_BC.L, ~Geno )
	cpBC.L<-plot_ordination( physeq_norm.L, cpcoa_BC.L, color="Geno", shape="Comp")
	cpBC.L$layers<-cpBC.L$layers[-1]
	p2=cpBC.L+geom_point(size=4)+theme_bw()+theme_new+
	scale_shape_manual(values=shape.2)+ 
	scale_color_manual(values=color_palette.2)
print("####################################################################")
print("Permanova test for the term genotype (i.e. Tae, Tbo, Tur) for leaf samples")
print(anova.cca(cpcoa_BC.L, by="term", perm.max=1000))
DT.L=data.table( as( sample_data( physeq_norm.L ), "data.frame" ),keep.rownames=T, key="rn" )
with( DT.L, adonis ( dist_BC.L ~ Geno ) )
print("####################################################################")

##for root samples	
	physeq_norm.R=subset_samples( physeq_norm, Comp == "Root" )
	dist_BC.R=distance( physeq_norm.R, dist1 )
	cpcoa_BC.R=ordinate( physeq_norm.R, meth2 , dist_BC.R, ~Geno )
	cpBC.R<-plot_ordination( physeq_norm.R, cpcoa_BC.R, color="Geno", shape="Comp")
	cpBC.R$layers<-cpBC.R$layers[-1]
	p3=cpBC.R+geom_point(size=4)+theme_bw()+theme_new+
	scale_shape_manual(values=shape.3)+ 
	scale_color_manual(values=color_palette.2)
print("####################################################################")
print("Permanova test for the term genotype (i.e. Tae, Tbo, Tur) for root samples")
print(anova.cca(cpcoa_BC.R, by="term", perm.max=1000))
DT.R=data.table( as( sample_data( physeq_norm.R ), "data.frame" ),keep.rownames=T, key="rn" )
with( DT.R, adonis ( dist_BC.R ~ Geno ) )
print("####################################################################")

##for rhizo. and Unpl. soil samples	
	physeq_norm.S=subset_samples( physeq_norm, !Comp %in% c("Root", "Leaf") )
	dist_BC.S=distance( physeq_norm.S, dist1 )
	cpcoa_BC.S=ordinate( physeq_norm.S, meth2 , dist_BC.S, ~Geno )
	cpBC.S<-plot_ordination( physeq_norm.S, cpcoa_BC.S, color="Geno", shape="Comp")
	cpBC.S$layers<-cpBC.S$layers[-1]
	p4=cpBC.S+geom_point(size=4)+theme_bw()+theme_new+
	scale_shape_manual(values=shape.4)+ 
	scale_color_manual(values=color_palette.1)
print("####################################################################")
print("Permanova test for the term genotype (i.e. Tae, Tbo, Tur) for root samples")
print(anova.cca(cpcoa_BC.S, by="term", perm.max=1000))
DT.S=data.table( as( sample_data( physeq_norm.S ), "data.frame" ),keep.rownames=T, key="rn" )
with( DT.S, adonis ( dist_BC.S ~ Geno ) )
print("####################################################################")

#pdf("Supplementary_Fig_4.pdf",paper="A4" ,useDingbats=FALSE)
gridExtra::grid.arrange(p1,p2,p3,p4, nrow=2, ncol=2)
#dev.off()

