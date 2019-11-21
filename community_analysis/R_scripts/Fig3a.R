####################################################################
#R scripts written by M. Amine Hassani - hassani.medamine@gmail.com#
####################################################################

#required packages
pkg=c("plyr", "ggplot2", "phyloseq", "data.table", "reshape2", "GUniFrac", "ape", "phytools", "metagenomeSeq","PMCMR","minpack.lm", "Hmisc", "reltools", "stats4")

lapply(pkg, library, character.only = TRUE)

rm(list = ls()[grep("*.*", ls())])

Palette.phylum <-c("#EF5656","#47B3DA","#F7A415","#88888a","#2BB065")
Palette.shape.1<-c(21, 21, 16)
Palette.shape.2<-c(17,15,19,13)
alpha=0.05
minZ=1.5
maxZ=8.5

Stat<-list()
p<-list()

theme_change <- theme(
	plot.background = element_blank(),
	panel.grid.minor = element_blank(),
	panel.grid.major = element_blank(),
	panel.border = element_rect(),
	axis.title.y=element_blank(),
	axis.title.x=element_blank(),
	panel.background = element_rect(fill="white", colour ="black")
	)
	
# uploading and preparing phyloseq object
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

#subset data to bacteria	
	physeq_B=subset_taxa( physeq, Kingdom == "Bacteria" )
	OTU_B=t(otu_table(physeq_B))
	TAXA_B=tax_table(physeq_B)
#transformation counts to rel. abund.
	physeq.ra = transform_sample_counts(physeq_B, function(x) x/sum(x))
	DT=data.table(otu_table(physeq.ra), keep.rownames=T, key="rn")

	L.Tae<-c('LTae1','LTae2','LTae3','LTae4','LTae5','LTae6')
	R.Tae<-c('RTae1','RTae2','RTae3','RTae4','RTae5','RTae6')
	Rz.Tae<-c('RzTae1','RzTae2','RzTae3','RzTae4','RzTae5','RzTae6')
	L.Tbo<-c('LTbo1','LTbo2','LTbo3','LTbo4','LTbo5','LTbo6')
	R.Tbo<-c('RTbo1','RTbo2','RTbo3','RTbo4','RTbo5','RTbo6')
	Rz.Tbo<-c('RzTbo1','RzTbo2','RzTbo3','RzTbo4','RzTbo5','RzTbo6')
	L.Tur<-c('LTur1','LTur2','LTur3','LTur4')
	R.Tur<-c('RTur1','RTur2','RTur3','RTur4')
	Rz.Tur<-c('RzTur1','RzTur2','RzTur3','RzTur4')
	UnSoil<-c('UnSoil1', 'UnSoil2', 'UnSoil3', 'UnSoil4')

	All<-list(L.Tae, L.Tbo, L.Tur,R.Tae, R.Tbo, R.Tur,Rz.Tae, Rz.Tbo, Rz.Tur, UnSoil)

	DT.taxa <- data.table(tax, keep.rownames=T, key="rn")
	DT.m=merge(DT,DT.taxa)

	ColTaxa<-colnames(DT.taxa)
	
for (i in All) {

	#Fiting Sloan Neutral Model
	print(i)
	DT.i<-DT.m[,c("rn",..i)]
	DT.i$mean<-DT.i[,.(rowMeans(.SD,na.rm = T)),.SDcols = i]
	LIST=DT.i[DT.i$mean!=0]$rn
		
	print(paste0("Computing SNCM for ",i[[1]], " samples ..."))	
	sncm.out=fit_sncm(spp=OTU_B[i,LIST],taxon=TAXA_B[LIST,])
	
	DT.sncm<-data.table(sncm.out$predictions, keep.rownames=T, key="rn")
	
	ColPhylum<-c("Actinobacteria","Bacteroidetes","Firmicutes","Proteobacteria")
	DT.sncm$ColPhylum <- ifelse(!DT.sncm$Phylum %in% ColPhylum, "Other", ifelse(DT.sncm$Phylum == "Bacteroidetes", "Bacteroidetes", ifelse(DT.sncm$Phylum == "Firmicutes", "Firmicutes", ifelse (DT.sncm$Phylum == "Actinobacteria", "Actinobacteria", "Proteobacteria"))))

	DT.sncm$ShpClass <- ifelse(DT.sncm$fit_class =="Below prediction", "Blw", ifelse(DT.sncm$fit_class == "Above prediction", "Abv", "Prd"))

	Abv<-length(DT.sncm[which(DT.sncm$ShpClass =="Abv"),]$rn)*100
	Blw<-length(DT.sncm[which(DT.sncm$ShpClass =="Blw"),]$rn)*100

	Richness<-sncm.out$fitstats$Richness

	fit.stat<-sncm.out$fitstats

	print(paste0("Percentage OTU Below = ", Blw/Richness))	
	print(paste0("Percentage OTU Above = ", Abv/Richness))

	pp1=ggplot(DT.sncm, aes(x=log(p), y=freq))

	Stat[[i[[1]]]]<-fit.stat
	print("######################################################################")
}

	sncm.stat<-rbind(Stat$LTae1, Stat$LTbo1, Stat$LTur1, Stat$RTae1, Stat$RTbo1, Stat$RTur1, Stat$RzTae1, Stat$RzTbo1, Stat$RzTur1, Stat$UnSoil1)
	rownames(sncm.stat)<-c("LTae", "LTbo", "LTur", "RTae", "RTbo", "RTur", "RzTae", "RzTbo", "RzTur", "UnSoil")
	sncm.stat$comp<-c("L","L","L","R","R","R","Rz","Rz","Rz", "UnSoil")

	mid<-mean(sncm.stat$Richness)

	pp<-ggplot(sncm.stat, aes(x=Rsqr, y=m, color=Richness, shape=comp))

	panel_a=pp+geom_point(size=5)+
	scale_color_gradient2(low="darkcyan", mid="darkorange", high="darkviolet", midpoint=mid)+
	theme_bw() + theme_change +
	scale_shape_manual(values=Palette.shape.2)

#	pdf("Fig3a.pdf", paper="A4", useDingbats=FALSE)
	print(panel_a)
#	dev.off()

##########################################################
###### figure was edited and merged with fig 3b    #######
##########################################################
















