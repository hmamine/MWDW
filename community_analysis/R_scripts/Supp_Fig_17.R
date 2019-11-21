####################################################################
#R scripts written by M. Amine Hassani - hassani.medamine@gmail.com#
####################################################################
 pkg=c("plyr", "ggplot2", "phyloseq", "data.table", "reshape2", "GUniFrac", "ape", "phytools", "metagenomeSeq","PMCMR","minpack.lm", "Hmisc", "reltools", "stats4","gridExtra")

lapply(pkg, library, character.only = TRUE)

rm(list = ls()[grep("*.*", ls())])

Palette.phylum <-c("#4e8e8e","#8e8e4e","#8e4e8e")
Palette.shape<-c(21, 21	, 16)
alpha=0.05
minZ=1.5
maxZ=8.5
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
	All<-list(L.Tae, L.Tbo, L.Tur,R.Tae, R.Tbo, R.Tur,Rz.Tae, Rz.Tbo, Rz.Tur, UnSoil)
#uploading data 
	mat=read.table( "otu_table_merge_fungi.txt", sep="\t", row.names=1, header=T)
	mat=as.matrix(mat)
	OTU=otu_table(mat, taxa_are_rows=T) 
	tax=read.table("taxonomy_fungi.txt", sep="\t", row.names=1, header=1)
	tax=as.matrix(tax)
	TAXA=tax_table(tax)
	sd=read.table("sample_data_merge.txt", sep="\t", row.names=1, header=1)
	SD=sample_data(sd)
	physeq= phyloseq(OTU, TAXA, SD)
#subset to fungi
	physeq_F=subset_taxa( physeq, Kingdom == "Fungi" )
#transform counts to rel. abund.
	physeq.ra = transform_sample_counts(physeq_F, function(x) x/sum(x))
	DT.ra=data.table(otu_table(physeq.ra), keep.rownames=T, key="rn")
	DT.taxa <- data.table(tax, keep.rownames=T, key="rn")
	ColPhylum<-c("Ascomycota","Basidiomycota")
	DT.taxa$ColPhylum <- ifelse(!DT.taxa$Phylum %in% ColPhylum, "Other", ifelse(DT.taxa$Phylum == "Ascomycota", "Ascomycota", "Basidiomycota"))
#computing occurance frequence and rel. abund.	
DT<-list()
p<-list()
for (i in All) 
{
	print (i)
	DT.I1<-DT.ra[,c("rn",..i)]
	for (col in i) set (DT.I1, j=col, value = +(DT.I1[[col]]>0))
	DT.I1$Frq<-DT.I1[,.(rowSums(.SD,na.rm = T)),.SDcols = i]
	DT.I1$Frq<-DT.I1$Frq/length(i)
	DT.I2<-DT.ra[,c("rn",..i)]
	DT.I2$mean<-DT.I2[,.(rowMeans(.SD,na.rm = T)),.SDcols = i]
	DT.I=DT.I1[DT.I2][DT.taxa]
	DT.I=DT.I[DT.I$mean!=0]
	DT[[i[[1]]]]<-DT.I
	pp1=ggplot(DT.I, aes(y=Frq, x=log2(mean), color=ColPhylum))
	p1=pp1+geom_jitter(position = position_jitter(width = 0.05, height = 0.05), size=1) +
	theme_bw()+ scale_color_manual(values=Palette.phylum)+
	geom_vline(xintercept = log2(0.01), linetype="longdash")+
	scale_y_continuous(limits=c(0, 1.1), breaks=c(0.00,.25,.5,.75,1.00))+
	theme_change
#	print(p1)
	p[[i[[1]]]]<-p1
	print ("#############################################")
}
#	pdf("Supplementary_Fig_17.pdf", paper="A4", useDingbats=FALSE)
	grid.arrange(grobs = p, ncol = 3)
#	dev.off()
