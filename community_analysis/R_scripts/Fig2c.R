####################################################################
#R scripts written by M. Amine Hassani - hassani.medamine@gmail.com#
####################################################################

#required packages
pkg=c("plyr", "ggplot2", "phyloseq", "data.table", "reshape2", "GUniFrac", "ape", "phytools", "metagenomeSeq","PMCMR", "seqtime")
lapply(pkg, library, character.only = TRUE)

rm(list = ls()[grep("*.*", ls())])

color_palette<-c("#cc2a36","#166590","#11902d","#000000")

theme_new <- theme (
	panel.grid.major = element_blank(),
	panel.grid.minor = element_blank(),
	legend.position="none",
	axis.text.x = element_text(angle=90, vjust=1)
	)

DT<-data.table(rn=c("S1","S2","S3","S4","S5"),  key="rn")
Mat<-list()
LIST<-list ( "Leaf", "Root", "RS" )
dist1=("bray")

#uploading and preparing phyloseq object
	mat=read.table( "otu_table_merge_fungi.txt", sep="\t", row.names=1, header=T)
	mat=as.matrix(mat)
	OTU=otu_table(mat, taxa_are_rows=T) 
	tax=read.table("taxonomy_fungi.txt", sep="\t", row.names=1, header=1)
	tax=as.matrix(tax)
	TAXA=tax_table(tax)
	sd=read.table("sample_data_merge.txt", sep="\t", row.names=1, header=1)
	SD=sample_data(sd)
	physeq= phyloseq(OTU, TAXA, SD)  

for( i in LIST ) 
	{
	print( i )
	SUBSET=c( i )
	physeq_sub=subset_samples( physeq, Comp2 %in% SUBSET )
	otu.I=otu_table(physeq_sub)
	taxa.I=tax_table(physeq_sub)
	sd.I=sample_data(physeq_sub)
	flt=filterTaxonMatrix(otu.I,minocc=1, keepSum = TRUE, return.filtered.indices = TRUE)	
	otu.I.filtered=flt$mat	
	taxa.I.filtered=taxa.I[setdiff(1:nrow(taxa.I),flt$filtered.indices),]	
	dummyTaxonomy=c("k__dummy","p__","c__","o__","f__","g__","s__")
	TAXA.filt=rbind(taxa.I.filtered, dummyTaxonomy)
	rownames(TAXA.filt)[nrow(TAXA.filt)]="SUM"	
	rownames(otu.I.filtered)[nrow(otu.I.filtered)]="SUM"
	NewOTU=otu_table(otu.I.filtered, taxa_are_rows = TRUE)
	NewTaxa=tax_table(TAXA.filt)	
	physeq.filter.I=phyloseq(NewOTU, NewTaxa, sd.I)

	#normalization of count reads using CSS 
	otumat=as( otu_table( physeq.filter.I ), "matrix" )	
	mp=newMRexperiment( ( otumat ) )	
	physeq_norm=phyloseq( otu_table( MRcounts( cumNorm( mp, p=cumNormStat( mp ) ),
	norm=TRUE, log=TRUE ),taxa_are_rows=TRUE ), NewTaxa, sd.I, TREE )

	#computing Bray-Curtis distances
	sd_tab=data.table( as( sample_data( physeq_norm ), "data.frame" ), keep.rownames=TRUE, key="rn" )
	dist_BC=distance( physeq_norm, dist1 )

	#computing BC distances to centroid
	mod <- with( sd_tab, betadisper( dist_BC, Concat2, type="centroid", bias.adjust =F ) )
	print(TukeyHSD(mod))
	p=boxplot(mod)
	dev.off()
	M1<-p$stats
	colnames(M1)<-p$names
	rownames(M1)<-c("S1","S2","S3","S4","S5")
	DT_i<-as.data.table(M1, keep.rownames=T, key="rn")
	DT=DT[DT_i,]
	Mat[[i]]<-as.matrix( dist_BC )
	}

	DT_sd=data.table( as( sample_data( physeq ), "data.frame" ), keep.rownames=TRUE)
	DT_M=melt(DT)
DT_M$grp<-ifelse(DT_M$variable %in% c("LTae", "RTae", "RzTae"), "Tae", ifelse (DT_M$variable %in% c("LTbo", "RTbo", "RzTbo"), "Tbo", ifelse(DT_M$variable %in% c("LTur", "RTur", "RzTur"), "Tur", "Unsoil")))

#preparing the plot
	panel_c=ggplot(DT_M, aes(x=variable, y=value,  color=grp))+
	geom_boxplot(width=1)+
	geom_point(size=1, color="black")+theme_bw()+ theme_new +
	scale_colour_manual(values=color_palette) + 
	ylab("Bray-Curtis distances to centroid")

	#pdf("Fig2c.pdf", useDingbats=FALSE, paper ="A4")
	gridExtra::grid.arrange(panel_c, nrow=2, ncol=2)
	#dev.off()

##########################################################
### figure was edited and merged with fig2 panels a and b#
##########################################################








