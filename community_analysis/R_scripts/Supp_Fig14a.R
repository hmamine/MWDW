####################################################################
#R scripts written by M. Amine Hassani - hassani.medamine@gmail.com#
####################################################################
pkg=c("plyr", "ggplot2", "phyloseq", "data.table", "reshape2", "GUniFrac", "ape", "phytools", "metagenomeSeq","PMCMR","minpack.lm", "Hmisc", "reltools", "stats4")

lapply(pkg, library, character.only = TRUE)

rm(list = ls()[grep("*.*", ls())])

val.Tae<-c('UnSoil.mean')
val.Tbo<-c('UnSoil.mean')
val.Tur<-c('UnSoil.mean')
color_palette<-c("#cc2a36","#166590","#11902d","#000000")
theme_change <- theme(
	legend.position="none",
	plot.background = element_blank(),
	panel.grid.minor = element_blank(),
	panel.grid.major = element_blank(),
	panel.border = element_rect(),
	axis.title.y=element_blank(),
	axis.title.x=element_blank(),
	panel.background = element_rect(fill="white", colour ="black")
	)

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
	
	sd=read.table("sample_data_merge.txt", sep="\t", row.names=1, header=1)
	SD=sample_data(sd)

#computing disparsal from soil
############
##Bacteria##
############
###uploading data 
	mat.B=read.table( "otu_table_merge_bacteria.txt", sep="\t", row.names=1, header=T)
	mat.B=as.matrix(mat.B)
	OTU.B=otu_table(mat.B, taxa_are_rows=T) 
	tax.B=read.table("taxonomy_bacteria.txt", sep="\t", row.names=1, header=1)
	tax.B=as.matrix(tax.B)
	TAXA.B=tax_table(tax.B)
	physeq_B= phyloseq(OTU.B, TAXA.B, SD)
###subset
	physeq_B=subset_taxa( physeq_B, Kingdom == "Bacteria" )
###transformt counts to rel. abund.
	physeq.B.ra = transform_sample_counts(physeq_B, function(x) x/sum(x))
###computing means
	DT.B.ra=data.table(otu_table(physeq.B.ra), keep.rownames=T, key="rn")
	DT.B.taxa <- data.table(tax.B, keep.rownames=T, key="rn")
	DT.B.m=merge(DT.B.ra, DT.B.taxa)
	DT.B.m$LTae.mean<-DT.B.m[,.(rowMeans(.SD,na.rm = T)),.SDcols = L.Tae]
	DT.B.m$RTae.mean<-DT.B.m[,.(rowMeans(.SD,na.rm = T)),.SDcols = R.Tae]
	DT.B.m$RzTae.mean<-DT.B.m[,.(rowMeans(.SD,na.rm = T)),.SDcols = Rz.Tae]
	DT.B.m$LTbo.mean<-DT.B.m[,.(rowMeans(.SD,na.rm = T)),.SDcols = L.Tbo]
	DT.B.m$RTbo.mean<-DT.B.m[,.(rowMeans(.SD,na.rm = T)),.SDcols = R.Tbo]
	DT.B.m$RzTbo.mean<-DT.B.m[,.(rowMeans(.SD,na.rm = T)),.SDcols = Rz.Tbo]	
	DT.B.m$LTur.mean<-DT.B.m[,.(rowMeans(.SD,na.rm = T)),.SDcols = L.Tur]
	DT.B.m$RTur.mean<-DT.B.m[,.(rowMeans(.SD,na.rm = T)),.SDcols = R.Tur]
	DT.B.m$RzTur.mean<-DT.B.m[,.(rowMeans(.SD,na.rm = T)),.SDcols = Rz.Tur]	
	DT.B.m$UnSoil.mean<-DT.B.m[,.(rowMeans(.SD,na.rm = T)),.SDcols = UnSoil]	
	DT.B<-DT.B.m[,-2:-54]
	DT_B<- DT.B[,-2:-7]
	COLS<-colnames(DT_B)[-1]
	for (col in COLS) set (DT_B, j=col, value = +(DT_B[[col]]>0))
###computing percentages
	BLTae=length(which(DT_B$LTae.mean == "1" & DT_B[,..val.Tae] == "1" ))/length(which(DT_B[,..val.Tae] == "1" ))
	BRTae=length(which(DT_B$RTae.mean == "1" & DT_B[,..val.Tae] == "1" ))/length(which(DT_B[,..val.Tae] == "1" ))
	BRzTae=length(which(DT_B$RzTae.mean == "1" & DT_B[,..val.Tae] == "1" ))/length(which(DT_B[,..val.Tae] == "1" ))
	BLTbo=length(which(DT_B$LTbo.mean == "1" & DT_B[,..val.Tbo] == "1" ))/length(which(DT_B[,..val.Tbo] == "1" ))
	BRTbo=length(which(DT_B$RTbo.mean == "1" & DT_B[,..val.Tbo] == "1" ))/length(which(DT_B[,..val.Tbo] == "1" ))
	BRzTbo=length(which(DT_B$RzTbo.mean == "1" & DT_B[,..val.Tbo] == "1" ))/length(which(DT_B[,..val.Tbo] == "1" ))
	BLTur=length(which(DT_B$LTur.mean == "1" & DT_B[,..val.Tur] == "1" ))/length(which(DT_B[,..val.Tur] == "1" ))
	BRTur=length(which(DT_B$RTur.mean == "1" & DT_B[,..val.Tur] == "1" ))/length(which(DT_B[,..val.Tur] == "1" ))	
	BRzTur=length(which(DT_B$RzTur.mean == "1" & DT_B[,..val.Tur] == "1" ))/length(which(DT_B[,..val.Tur] == "1" ))
###preparing plot
	disp.B<-rbind(BLTae, BLTbo,BLTur,BRTae,BRTbo,BRTur,BRzTae,BRzTbo,BRzTur)
	colnames(disp.B)<-"value"
	disp.B<-data.frame(disp.B)
	disp.B$concat<-c("LTae", "LTbo", "LTur", "RTae", "RTbo", "RTur", "RzTae", "RzTbo", "RzTur")
	disp.B$comp<-c("L","L","L","R","R","R","Rz","Rz","Rz")
	disp.B$kingdom<-c(rep.int("B",9))
	disp.B$geno<-c("Tae","Tbo","Tur","Tae","Tbo","Tur","Tae","Tbo","Tur")

#########
##Fungi##
#########
###uploading data
	mat.F=read.table( "otu_table_merge_fungi.txt", sep="\t", row.names=1, header=T)
	mat.F=as.matrix(mat.F)
	OTU.F=otu_table(mat.F, taxa_are_rows=T) 
	tax.F=read.table("taxonomy_fungi.txt", sep="\t", row.names=1, header=1)
	tax.F=as.matrix(tax.F)
	TAXA.F=tax_table(tax.F)
	physeq_F= phyloseq(OTU.F, TAXA.F, SD)
###subset to Fungi	
	physeq_F=subset_taxa( physeq_F, Kingdom == "Fungi" )
###transformt counts to rel. abund.
	physeq.F.ra = transform_sample_counts(physeq_F, function(x) x/sum(x))
###computing means
	DT.F.ra=data.table(otu_table(physeq.F.ra), keep.rownames=T, key="rn")
	DT.F.taxa <- data.table(tax.F, keep.rownames=T, key="rn")
	DT.F.m=merge(DT.F.ra, DT.F.taxa)
	DT.F.m$LTae.mean<-DT.F.m[,.(rowMeans(.SD,na.rm = T)),.SDcols = L.Tae]
	DT.F.m$RTae.mean<-DT.F.m[,.(rowMeans(.SD,na.rm = T)),.SDcols = R.Tae]
	DT.F.m$RzTae.mean<-DT.F.m[,.(rowMeans(.SD,na.rm = T)),.SDcols = Rz.Tae]
	DT.F.m$LTbo.mean<-DT.F.m[,.(rowMeans(.SD,na.rm = T)),.SDcols = L.Tbo]
	DT.F.m$RTbo.mean<-DT.F.m[,.(rowMeans(.SD,na.rm = T)),.SDcols = R.Tbo]
	DT.F.m$RzTbo.mean<-DT.F.m[,.(rowMeans(.SD,na.rm = T)),.SDcols = Rz.Tbo]
	DT.F.m$LTur.mean<-DT.F.m[,.(rowMeans(.SD,na.rm = T)),.SDcols = L.Tur]
	DT.F.m$RTur.mean<-DT.F.m[,.(rowMeans(.SD,na.rm = T)),.SDcols = R.Tur]
	DT.F.m$RzTur.mean<-DT.F.m[,.(rowMeans(.SD,na.rm = T)),.SDcols = Rz.Tur]
	DT.F.m$UnSoil.mean<-DT.F.m[,.(rowMeans(.SD,na.rm = T)),.SDcols = UnSoil]
	DT.F<-DT.F.m[,-2:-54]
	DT_F<- DT.F[,-2:-7]
	COLS<-colnames(DT_F)[-1]
	for (col in COLS) set (DT_F, j=col, value = +(DT_F[[col]]>0))
###computing percentages
	FLTae=length(which(DT_F$LTae.mean == "1" & DT_F[,..val.Tae] == "1" ))/length(which(DT_F[,..val.Tae] == "1" ))
	FRTae=length(which(DT_F$RTae.mean == "1" & DT_F[,..val.Tae] == "1" ))/length(which(DT_F[,..val.Tae] == "1" ))
	FRzTae=length(which(DT_F$RzTae.mean == "1" & DT_F[,..val.Tae] == "1" ))/length(which(DT_F[,..val.Tae] == "1" ))
	FLTbo=length(which(DT_F$LTbo.mean == "1" & DT_F[,..val.Tbo] == "1" ))/length(which(DT_F[,..val.Tbo] == "1" ))
	FRTbo=length(which(DT_F$RTbo.mean == "1" & DT_F[,..val.Tbo] == "1" ))/length(which(DT_F[,..val.Tbo] == "1" ))
	FRzTbo=length(which(DT_F$RzTbo.mean == "1" & DT_F[,..val.Tbo] == "1" ))/length(which(DT_F[,..val.Tbo] == "1" ))
	FLTur=length(which(DT_F$LTur.mean == "1" & DT_F[,..val.Tur] == "1" ))/length(which(DT_F[,..val.Tur] == "1" ))
	FRTur=length(which(DT_F$RTur.mean == "1" & DT_F[,..val.Tur] == "1" ))/length(which(DT_F[,..val.Tur] == "1" ))
	FRzTur=length(which(DT_F$RzTur.mean == "1" & DT_F[,..val.Tur] == "1" ))/length(which(DT_F[,..val.Tur] == "1" ))
###preparing plots
	disp.F<-rbind(FLTae, FLTbo,FLTur,FRTae, FRTbo,FRTur,FRzTae, FRzTbo,FRzTur)
	colnames(disp.F)<-"value"
	disp.F<-data.frame(disp.F)
	disp.F$concat<-c("LTae", "LTbo", "LTur", "RTae", "RTbo", "RTur", "RzTae", "RzTbo", "RzTur")
	disp.F$comp<-c("L","L","L","R","R","R","Rz","Rz","Rz")
	disp.F$kingdom<-c(rep.int("F",9))
	disp.F$geno<-c("Tae","Tbo","Tur","Tae","Tbo","Tur","Tae","Tbo","Tur")

###combining bacterial and fungal plots
	DT.disp<-data.table(rbind(disp.B, disp.F), keep.rownames=T, key="rn")
	DT.disp$concat1=paste(DT.disp$comp, DT.disp$kingdom, sep=".")
	pp1<-ggplot(DT.disp, aes(x=concat1, y=value))
	panel_a=pp1+geom_boxplot(colour="Black", width=1)+
	geom_jitter(aes(colour=geno),size=2)+ 
	scale_colour_manual(values=color_palette) +
	theme_bw() + theme_change 
#	pdf("Supplementary_Fig_14a.pdf",paper="A4" ,useDingbats=FALSE)
	gridExtra::grid.arrange(panel_a, nrow=2, ncol=2)
#	dev.off()

###testing significance
DT.L<-DT.disp[DT.disp$comp == "L",]
summary(out.B<-aov(data=DT.L, value ~ as.factor(kingdom)))

DT.R<-DT.disp[DT.disp$comp == "R",]
summary(out.B<-aov(data=DT.R, value ~ as.factor(kingdom)))

DT.Rz<-DT.disp[DT.disp$comp == "Rz",]
summary(out.B<-aov(data=DT.Rz, value ~ as.factor(kingdom)))









