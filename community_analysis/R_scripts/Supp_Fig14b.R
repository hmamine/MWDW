####################################################################
#R scripts written by M. Amine Hassani - hassani.medamine@gmail.com#
####################################################################
pkg=c("plyr", "ggplot2", "phyloseq", "data.table", "reshape2", "GUniFrac", "ape", "phytools", "metagenomeSeq","PMCMR","minpack.lm", "Hmisc", "reltools", "stats4")

lapply(pkg, library, character.only = TRUE)

rm(list = ls()[grep("*.*", ls())])

color_palette<-c("#cc2a36","#166590","#11902d","#000000")
B<-list()
F<-list()

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
	
	LIST<-list (L.Tae,R.Tae,Rz.Tae,L.Tbo,R.Tbo,Rz.Tbo,L.Tur,R.Tur,Rz.Tur)
	
	sd=read.table("sample_data_merge.txt", sep="\t", row.names=1, header=1)
	SD=sample_data(sd)

#computing within group dispersal
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
###subset to bacteria	
	physeq_B=subset_taxa( physeq_B, Kingdom == "Bacteria" )
###transform counts to rel. abundance
	physeq.B.ra = transform_sample_counts(physeq_B, function(x) x/sum(x))
###computing means
	DT.B.ra=data.table(otu_table(physeq.B.ra), keep.rownames=T, key="rn")
	DT.B=DT.B.ra
	COLS<-colnames(DT.B)[-1]
	for (col in COLS) set (DT.B, j=col, value = +(DT.B[[col]]>0))
###computing dispersal rates
	for( i in LIST ) 
	{
		print (paste0("computing within group dispesal for  ", i))
		DT_i<-DT.B[,c("rn",..i)]
		x<-length(colnames(DT_i))-1
		print(x)
		DT_i$Frq<-DT_i[,.(rowSums(.SD,na.rm = T)),.SDcols = i]/ x
		DT_i<-DT_i[which(!Frq == "0")]
		y<-sum(DT_i$Frq)/length(DT_i$rn)
		print (paste0 ("disparsal rate = ",y))
		B[i[[1]]]<- y
		print ("#########################################################")
	}
###preparing plots
	disp.B<-rbind(B$LTae1, B$LTbo1,B$LTur1,B$RTae1,B$RTbo1,B$RTur1,B$RzTae1,B$RzTbo1,B$RzTur1)
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
###transform counts to rel. abundance
	physeq.F.ra = transform_sample_counts(physeq_F, function(x) x/sum(x))
###computing means
	DT.F.ra=data.table(otu_table(physeq.F.ra), keep.rownames=T, key="rn")
	DT.F=DT.F.ra
	COLS<-colnames(DT.F)[-1]
	for (col in COLS) set (DT.F, j=col, value = +(DT.F[[col]]>0))
###computing dispersal rates
	for( i in LIST ) 
	{
		print (paste0("computing within group dispesal for  ", i))
		DT_i<-DT.F[,c("rn",..i)]
		x<-length(colnames(DT_i))-1
		print(x)
		DT_i$Frq<-DT_i[,.(rowSums(.SD,na.rm = T)),.SDcols = i]/ x
		DT_i<-DT_i[which(!Frq == "0")]
		y<-sum(DT_i$Frq)/length(DT_i$rn)
		print (paste0 ("disparsal rate = ", y))
		F[i[[1]]]<- y
		print ("#########################################################")
	}
###preparing plots
	disp.F<-rbind(F$LTae1, F$LTbo1,F$LTur1,F$RTae1,F$RTbo1,F$RTur1,F$RzTae1,F$RzTbo1,F$RzTur1)
	colnames(disp.F)<-"value"
	disp.F<-data.frame(disp.F)
	disp.F$concat<-c("LTae", "LTbo", "LTur", "RTae", "RTbo", "RTur", "RzTae", "RzTbo", "RzTur")
	disp.F$comp<-c("L","L","L","R","R","R","Rz","Rz","Rz")
	disp.F$kingdom<-c(rep.int("F",9))
	disp.F$geno<-c("Tae","Tbo","Tur","Tae","Tbo","Tur","Tae","Tbo","Tur")

###combiing both plots
	DT.disp<-data.table(rbind(disp.B, disp.F), keep.rownames=T, key="rn")
	DT.disp$concat1=paste(DT.disp$comp, DT.disp$kingdom, sep=".")
	pp1<-ggplot(DT.disp, aes(x=concat1, y=value))
	p1=pp1+geom_boxplot(colour="Black", width=1)+
	geom_jitter(aes(colour=geno),size=2)+ 
	scale_colour_manual(values=color_palette) +
	theme_bw() + theme_change 
#	pdf("Supplementary_Fig_14b.pdf",paper="A4" ,useDingbats=FALSE)
	gridExtra::grid.arrange(p1, nrow=2, ncol=2)
#	dev.off()

###testing significance
DT.L<-DT.disp[DT.disp$comp == "L",]
summary(out.B<-aov(data=DT.L, value ~ as.factor(kingdom)))

DT.R<-DT.disp[DT.disp$comp == "R",]
summary(out.B<-aov(data=DT.R, value ~ as.factor(kingdom)))

DT.Rz<-DT.disp[DT.disp$comp == "Rz",]
summary(out.B<-aov(data=DT.Rz, value ~ as.factor(kingdom)))









