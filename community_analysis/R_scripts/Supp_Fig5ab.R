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
dist1=("bray")
dist2=("wunifrac")
# uploading the data
	mat=read.table( "otu_table_merge_bacteria.txt", sep="\t", row.names=1, header=T)
	mat=as.matrix(mat)
	OTU=otu_table(mat, taxa_are_rows=T) 
	tax=read.table("taxonomy_bacteria.txt", sep="\t", row.names=1, header=1)
	tax=as.matrix(tax)
	TAXA=tax_table(tax)
	sd=read.table("sample_data_merge.txt", sep="\t", row.names=1, header=1)
	SD=sample_data(sd)
	TREE=read.tree("tree.nwk")
	physeq= phyloseq(OTU, TAXA, SD,TREE)  
#subset samples to habitat
Habitat<-list("Leaf", "Root", "RS" )
BC<-list()
WUF<-list()
SD<-list()
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
	ColPhylum<-c("Actinobacteria","Bacteroidetes","Firmicutes","Proteobacteria")
TAXA.filt[,"Phylum"] <- ifelse(! TAXA.filt[,"Phylum"] %in% ColPhylum, "Other", ifelse(TAXA.filt[,"Phylum"] == "Bacteroidetes", "Bacteroidetes", ifelse(TAXA.filt[,"Phylum"] == "Firmicutes", "Firmicutes", ifelse (TAXA.filt[,"Phylum"] == "Actinobacteria", "Actinobacteria", "Proteobacteria"))))
	rownames(TAXA.filt)[nrow(TAXA.filt)]="SUM"	
	rownames(otu.i.filtered)[nrow(otu.i.filtered)]="SUM"
	NewOTU=otu_table(otu.i.filtered, taxa_are_rows = TRUE)	
	NewTaxa=tax_table(TAXA.filt)
	physeq.filter.i=phyloseq(NewOTU, NewTaxa, sd.i)
	#normalization of count reads using CSS 
	otumat=as( otu_table( physeq.filter.i ), "matrix" )	
	mp=newMRexperiment( ( otumat ) )	
	physeq_norm_i=phyloseq( otu_table( MRcounts( cumNorm( mp, p=cumNormStat( mp ) ),
	norm=TRUE, log=TRUE ),taxa_are_rows=TRUE ), NewTaxa, sd.i, TREE )
	#computing Bray-Curtis distances
	dist_BC=distance( physeq_norm_i, dist1 )
	BC[[i]]<-as.matrix( dist_BC )
	#computing weighted UniFrac distances
	dist_WUF=distance( physeq_norm_i, dist2 )
	WUF[[i]]<-as.matrix( dist_WUF )

	SD[[i]]<-data.frame( sd.i )
}
#Preparing the plots
##leaf samples
	print("##################################################")
	LTae<-c("LTae1","LTae2","LTae3","LTae4","LTae5","LTae6")
	LTbo<-c("LTbo1","LTbo2","LTbo3","LTbo4","LTbo5","LTbo6")
	LTur<-c("LTur1","LTur2","LTur3","LTur4")	
	SD_leaf=data.table( SD$Leaf, keep.rownames=T, key="rn" )
###for Bray Curtis distances
	BC_leaf<-as.matrix(BC$Leaf)
	BC_leaf.sub<-BC_leaf[LTae, c(LTbo, LTur)]
	DT_BC_leaf=data.table( melt(BC_leaf.sub), keep.rownames=T, key=c("Var1","Var2") )
	DT_BC_leaf=data.table( DT_BC_leaf, keep.rownames=T, key=c("Var1","Var2") )
	DT_BC_leaf=merge( DT_BC_leaf, SD_leaf, by.x="Var2", by.y="rn" )
	DT_BC_leaf=merge( DT_BC_leaf, SD_leaf, by.x="Var1", by.y="rn" )
	pLBC=ggplot( DT_BC_leaf, aes( x=Geno.y, y=value, color=Geno.x ) )+ 
	geom_boxplot()+ theme_bw()+ 
	theme_new + scale_color_manual(values=color_palette.1)+
	ylab("Bray-Curtis distances to T. aestivum")
	#testing for significance
	print("testing significance - leaf samples - BC distances")
	print(kruskal.test(data=pLBC$data, value ~ Geno.x))
	print("##################################################")
###for Weighted UniFrac distances
	WUF_leaf<-as.matrix(WUF$Leaf)
	WUF_leaf.sub<-WUF_leaf[LTae, c(LTbo, LTur)]
	DT_WUF_leaf=data.table( melt(WUF_leaf.sub), keep.rownames=T, key=c("Var1","Var2") )
	DT_WUF_leaf=data.table( DT_WUF_leaf, keep.rownames=T, key=c("Var1","Var2") )
	DT_WUF_leaf=merge( DT_WUF_leaf, SD_leaf, by.x="Var2", by.y="rn" )
	DT_WUF_leaf=merge( DT_WUF_leaf, SD_leaf, by.x="Var1", by.y="rn" )
	pLWUF=ggplot( DT_WUF_leaf, aes( x=Geno.y, y=value, color=Geno.x ) )+ 
	geom_boxplot()+ theme_bw()+ 
	theme_new + scale_color_manual(values=color_palette.1)+
	ylab("Weighted UniFrac distances to T. aestivum")
	#testing for significance
	print("testing significance - leaf samples - WUF distances")
	print(kruskal.test(data=pLWUF$data, value ~ Geno.x))
	print("##################################################")
##root samples
	print("##################################################")
	RTae<-c("RTae1","RTae2","RTae3","RTae4","RTae5","RTae6")
	RTbo<-c("RTbo1","RTbo2","RTbo3","RTbo4","RTbo5","RTbo6")
	RTur<-c("RTur1","RTur2","RTur3","RTur4")	
	SD_root=data.table( SD$Root, keep.rownames=T, key="rn" )
###for Bray Curtis distances
	BC_root<-as.matrix(BC$Root)
	BC_root.sub<-BC_root[RTae, c(RTbo, RTur)]
	DT_BC_root=data.table( melt(BC_root.sub), keep.rownames=T, key=c("Var1","Var2") )
	DT_BC_root=data.table( DT_BC_root, keep.rownames=T, key=c("Var1","Var2") )
	DT_BC_root=merge( DT_BC_root, SD_root, by.x="Var2", by.y="rn" )
	DT_BC_root=merge( DT_BC_root, SD_root, by.x="Var1", by.y="rn" )
	pRBC=ggplot( DT_BC_root, aes( x=Geno.y, y=value, color=Geno.x ) )+ 
	geom_boxplot()+ theme_bw()+ 
	theme_new + scale_color_manual(values=color_palette.1)+
	ylab("Bray-Curtis distances to T. aestivum")
	#testing for significance
	print("testing significance - root samples - BC distances")
	print(kruskal.test(data=pRBC$data, value ~ Geno.x))
	print("##################################################")
###for Weighted UniFrac distances
	WUF_root<-as.matrix(WUF$Root)
	WUF_root.sub<-WUF_root[RTae, c(RTbo, RTur)]
	DT_WUF_root=data.table( melt(WUF_root.sub), keep.rownames=T, key=c("Var1","Var2") )
	DT_WUF_root=data.table( DT_WUF_root, keep.rownames=T, key=c("Var1","Var2") )
	DT_WUF_root=merge( DT_WUF_root, SD_root, by.x="Var2", by.y="rn" )
	DT_WUF_root=merge( DT_WUF_root, SD_root, by.x="Var1", by.y="rn" )
	pRWUF=ggplot( DT_WUF_root, aes( x=Geno.y, y=value, color=Geno.x ) )+ 
	geom_boxplot()+ theme_bw()+ 
	theme_new + scale_color_manual(values=color_palette.1)+
	ylab("Weighted UniFrac distances to T. aestivum")
	#testing for significance
	print("testing significance - root samples - WUF distances")
	print(kruskal.test(data=pRWUF$data, value ~ Geno.x))
	print("##################################################")
##rhizo. and unpl. samples
	print("##################################################")
	RzTae<-c("RzTae1","RzTae2","RzTae3","RzTae4","RzTae5","RzTae6")
	RzTbo<-c("RzTbo1","RzTbo2","RzTbo3","RzTbo4","RzTbo5","RzTbo6")
	RzTur<-c("RzTur1","RzTur2","RzTur3","RzTur4")
	UnSoil<-c("UnSoil1","UnSoil2","UnSoil3","UnSoil4")
	SD_rs=data.table( SD$RS, keep.rownames=T, key="rn" )
###for Bray Curtis distances
	BC_rs<-as.matrix(BC$RS)
	BC_rs.sub<-BC_rs[RzTae, c(UnSoil, RzTbo, RzTur)]
	DT_BC_rs=data.table( melt(BC_rs.sub), keep.rownames=T, key=c("Var1","Var2") )
	DT_BC_rs=data.table( DT_BC_rs, keep.rownames=T, key=c("Var1","Var2") )
	DT_BC_rs=merge( DT_BC_rs, SD_rs, by.x="Var2", by.y="rn" )
	DT_BC_rs=merge( DT_BC_rs, SD_rs, by.x="Var1", by.y="rn" )
	pRSBC=ggplot( DT_BC_rs, aes( x=Geno.y, y=value, color=Geno.x ) )+ 
	geom_boxplot()+ theme_bw()+ 
	theme_new + scale_color_manual(values=color_palette.2)+
	ylab("Bray-Curtis distances to T. aestivum")
	#testing for significance
	print("testing significance - root samples - BC distances")
	print(kruskal.test(data=pRSBC$data, value ~ Geno.x))
	print("##################################################")
###for Weighted UniFrac distances
	WUF_rs<-as.matrix(WUF$RS)
	WUF_rs.sub<-WUF_rs[RzTae, c(UnSoil, RzTbo, RzTur)]
	DT_WUF_rs=data.table( melt(WUF_rs.sub), keep.rownames=T, key=c("Var1","Var2") )
	DT_WUF_rs=data.table( DT_WUF_rs, keep.rownames=T, key=c("Var1","Var2") )
	DT_WUF_rs=merge( DT_WUF_rs, SD_rs, by.x="Var2", by.y="rn" )
	DT_WUF_rs=merge( DT_WUF_rs, SD_rs, by.x="Var1", by.y="rn" )
	pRSWUF=ggplot( DT_WUF_rs, aes( x=Geno.y, y=value, color=Geno.x ) )+ 
	geom_boxplot()+ theme_bw()+ 
	theme_new + scale_color_manual(values=color_palette.2)+
	ylab("Weighted UniFrac distances to T. aestivum")
	#testing for significance
	print("testing significance - root samples - WUF distances")
	print(kruskal.test(data=pRSWUF$data, value ~ Geno.x))
	print("##################################################")

#	gridExtra::grid.arrange(pLBC, pLWUF, pRBC, pRWUF, pRSBC, pRSWUF, ncol=2)

	DT_BC_comb=rbind(pLBC$data, pRBC$data, pRSBC$data)
	panel_a=ggplot( DT_BC_comb, aes( x=Concat2.y, y=value, color=Geno.x ) )+ 
	geom_boxplot(width=1)+ theme_bw()+ 
	theme_new + scale_color_manual(values=color_palette.2)+
	ylab("Bray-Curtis distances to T. aestivum")

	DT_WUF_comb=rbind(pLWUF$data, pRWUF$data, pRSWUF$data)
	panel_b=ggplot( DT_WUF_comb, aes( x=Concat2.y, y=value, color=Geno.x ) )+ 
	geom_boxplot(width=1)+ theme_bw()+ 
	theme_new + scale_color_manual(values=color_palette.2)+
	ylab("Weighted UniFrac distances to T. aestivum")
	
	#pdf("Supplementary_Fig_5ab.pdf",paper="A4" ,useDingbats=FALSE)
	gridExtra::grid.arrange(panel_a, panel_b, nrow=2)
	#dev.off()

