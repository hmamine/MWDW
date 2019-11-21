####################################################################
#R scripts written by M. Amine Hassani - hassani.medamine@gmail.com#
####################################################################
#required packages
pkg=c("plyr", "ggplot2", "phyloseq", "data.table", "reshape2", "GUniFrac", "ape", "phytools", "metagenomeSeq","PMCMR", "seqtime")
lapply(pkg, library, character.only = TRUE)

rm(list = ls()[grep("*.*", ls())])

color_palette.1<-c("#166590","#11902d")
color_palette.2<-c("#166590","#11902d", "#000000")
theme_new <- theme (
	panel.border = element_rect(),
	panel.grid.major = element_blank(),
	panel.grid.minor = element_blank(),
	legend.position="none",
	axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank()
	)
dist1=(	"bray"	)
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
DT<-data.table(rn=c("S1","S2","S3","S4","S5"),  key="rn")
Mat<-list()
SD<-list()
Habitat<-list ( "Leaf","Root", "RS" )
for( i in Habitat ) 
{
	print( i )
	subset.i=c( i )
	physeq_i=subset_samples( physeq, Comp2 %in% subset.i )
	otu.i=otu_table(physeq_i)
	taxa.i=tax_table(physeq_i)
	sd.i=sample_data(physeq_i)
	flt=filterTaxonMatrix(otu.i,minocc=3, keepSum = TRUE, return.filtered.indices = TRUE)	
	otu.i.filtered=flt$mat	
	taxa.i.filtered=taxa.i[setdiff(1:nrow(taxa.i),flt$filtered.indices),]	
	dummyTaxonomy=c("k__dummy","p__","c__","o__","f__","g__","s__")
	TAXA.filt=rbind(taxa.i.filtered, dummyTaxonomy)
	rownames(TAXA.filt)[nrow(TAXA.filt)]="SUM"
	rownames(otu.i.filtered)[nrow(otu.i.filtered)]="SUM"
	NewOTU=otu_table(otu.i.filtered, taxa_are_rows = TRUE)
	NewTaxa=tax_table(TAXA.filt)
	physeq.filter.i=phyloseq(NewOTU, NewTaxa, sd.i)
	#normalization of count reads using CSS 
	otumat=as( otu_table( physeq.filter.i ), "matrix" )	
	mp=newMRexperiment( ( otumat ) )	
	physeq_norm_i=phyloseq( otu_table( MRcounts( cumNorm( mp, p=cumNormStat( mp ) ),
	norm=TRUE, log=TRUE ),taxa_are_rows=TRUE ), NewTaxa, sd.i )
	#computing Bray-Curtis distances
	dist_BC=distance( physeq_norm_i, dist1 )
	Mat[[i]]<-as.matrix( dist_BC )
	SD[[i]]<-data.frame( sd.i )
}
#Preparing the plots
##leaf samples
	print("##################################################")
	Mat_leaf<-as.matrix(Mat$Leaf)
	LTae<-c("LTae1","LTae2","LTae3","LTae4","LTae5","LTae6")
	LTbo<-c("LTbo1","LTbo2","LTbo3","LTbo4","LTbo5","LTbo6")
	LTur<-c("LTur1","LTur2","LTur3","LTur4")
	Mat_leaf.sub<-Mat_leaf[LTae, c(LTbo, LTur)]
	DT_leaf=data.table( melt(Mat_leaf.sub), keep.rownames=T, key=c("Var1","Var2") )
	DT_leaf=data.table( DT_leaf, keep.rownames=T, key=c("Var1","Var2") )
	SD_leaf=data.table( SD$Leaf, keep.rownames=T, key="rn" )
	DT_leaf=merge( DT_leaf, SD_leaf, by.x="Var2", by.y="rn" )
	DT_leaf=merge( DT_leaf, SD_leaf, by.x="Var1", by.y="rn" )
	p1=ggplot( DT_leaf, aes( x=Geno.y, y=value, color=Geno.x ) )+ 
	geom_boxplot()+ theme_bw()+ 
	theme_new + scale_color_manual(values=color_palette.1)+
	ylab("Bray-Curtis distances") + xlab("T. aestivum")
	print("Testing significance - leaf samples")
	print(kruskal.test(data=p1$data, value ~ Geno.x))
	print("##################################################")
##root samples
	print("##################################################")
	Mat_root<-as.matrix(Mat$Root)
	RTae<-c("RTae1","RTae2","RTae3","RTae4","RTae5","RTae6")
	RTbo<-c("RTbo1","RTbo2","RTbo3","RTbo4","RTbo5","RTbo6")
	RTur<-c("RTur1","RTur2","RTur3","RTur4")
	Mat_root.sub<-Mat_root[RTae, c(RTbo, RTur)]
	DT_root=data.table( melt(Mat_root.sub), keep.rownames=T, key=c("Var1","Var2") )
	DT_root=data.table( DT_root, keep.rownames=T, key=c("Var1","Var2") )
	SD_root=data.table( SD$Root, keep.rownames=T, key="rn" )
	DT_root=merge( DT_root, SD_root, by.x="Var2", by.y="rn" )
	DT_root=merge( DT_root, SD_root, by.x="Var1", by.y="rn" )
	p2=ggplot( DT_root, aes( x=Geno.y, y=value, color=Geno.x ) )+ 
	geom_boxplot()+ theme_bw()+ 
	theme_new + scale_color_manual(values=color_palette.1)+
	ylab("Bray-Curtis distances") + xlab("T. aestivum")
	print("Testing significance - root samples")
	print(kruskal.test(data=p2$data, value ~ Geno.x))
	print("##################################################")
# Soil and Rhizo. samples
	print("##################################################")
	Mat_rs<-as.matrix(Mat$RS)
	RzTae<-c("RzTae1","RzTae2","RzTae3","RzTae4","RzTae5","RzTae6")
	RzTbo<-c("RzTbo1","RzTbo2","RzTbo3","RzTbo4","RzTbo5","RzTbo6")
	RzTur<-c("RzTur1","RzTur2","RzTur3","RzTur4")
	UnSoil<-c("UnSoil1","UnSoil2","UnSoil3","UnSoil4")
	Mat_rs.sub<-Mat_rs[RzTae, c(UnSoil, RzTbo, RzTur)]
	DT_rs=data.table( melt(Mat_rs.sub), keep.rownames=T, key=c("Var1","Var2") )
	DT_rs=data.table( DT_rs, keep.rownames=T, key=c("Var1","Var2") )
	SD_rs=data.table( SD$RS, keep.rownames=T, key="rn" )
	DT_rs=merge( DT_rs, SD_rs, by.x="Var2", by.y="rn" )
	DT_rs=merge( DT_rs, SD_rs, by.x="Var1", by.y="rn" )
	p3=ggplot( DT_rs, aes( x=Geno.y, y=value, color=Geno.x ) )+ 
	geom_boxplot()+ theme_bw()+ 
	theme_new + scale_color_manual(values=color_palette.2)+
	ylab("Bray-Curtis distances") + xlab("T. aestivum")
	print("Testing significance - soil and rhizo. samples")
	print(kruskal.test(data=p3$data, value ~ Geno.x))
	print(posthoc.kruskal.conover.test(data=p3$data, value ~ Geno.x, method="BH" ))
	print("##################################################")
#	gridExtra::grid.arrange(p1,p2,p3, nrow=2, ncol=2)

	DT_comb=rbind(p1$data, p2$data, p3$data)
	panel_c=ggplot( DT_comb, aes( x=Concat2.y, y=value, color=Geno.x ) )+ 
	geom_boxplot(width=1)+ theme_bw()+ 
	theme_new + scale_color_manual(values=color_palette.2)+
	ylab("Bray-Curtis distances") + xlab("T. aestivum")
#	pdf("Supplementary_Fig_5c",paper="A4" ,useDingbats=FALSE)
	gridExtra::grid.arrange(panel_c, nrow=2, ncol=2)
#	dev.off()

