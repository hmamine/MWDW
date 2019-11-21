####################################################################
#R scripts written by M. Amine Hassani - hassani.medamine@gmail.com#
####################################################################
pkg=c("plyr", "ggplot2", "phyloseq", "data.table", "reshape2", "GUniFrac", "ape", "phytools", "metagenomeSeq","PMCMR", "venn", "magrittr", "pracma", "seqtime", "gridGraphics")

lapply(pkg, library, character.only = TRUE)

rm(list = ls()[grep("*.*", ls())])

alpha=0.05
theme_new <- theme (
	panel.border = element_rect(),
	panel.grid.major = element_blank(),
	panel.grid.minor = element_blank(),
	legend.position="none",
	axis.title.x=element_blank(),
	axis.title.y=element_blank(),
#       axis.text.x=element_blank(),
#       axis.ticks.x=element_blank(),
	panel.background = element_rect(fill="white", colour ="black", size = 1)
	)

Shape.1=c(21,19)
Tae<-list()
Tbo<-list()
Tur<-list()

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
#computing enriched OTUs
res<-vector( "list" )
list.Comp<-list ( "Leaf","Root", "Rhizo" )
for( i in list.Comp ) 
{	
	#subset to leaf samples
	print( i )
	subset.i=c( i )
	physeq.sub=subset_samples( physeq_F, Comp %in% subset.i )
	#filtering OTUs with by occurance
	otu.i=otu_table(physeq.sub)
	taxa.i=tax_table(physeq.sub)
	sd.i=sample_data(physeq.sub)
	flt=filterTaxonMatrix(otu.i,minocc=2, keepSum = TRUE, return.filtered.indices = TRUE)
	otu.i.filtered=flt$mat
	taxa.i.filtered=taxa.i[setdiff(1:nrow(taxa.i),flt$filtered.indices),]
	dummyTaxonomy=c("k__dummy","p__","c__","o__","f__","g__","s__", "sc__")
	TAXA.filt=rbind(taxa.i.filtered, dummyTaxonomy)
	rownames(TAXA.filt)[nrow(TAXA.filt)]="SUM"
	rownames(otu.i.filtered)[nrow(otu.i.filtered)]="SUM"
	NewOTU=otu_table(otu.i.filtered, taxa_are_rows = TRUE)
	NewTaxa=tax_table(TAXA.filt)
	physeq.filter=phyloseq(NewOTU, NewTaxa, sd.i)
	DT.tax=data.table(TAXA.filt, keep.rownames=T, key="rn")
	ColPhylum<-c("Ascomycota","Basidiomycota")
	DT.tax$ColPhylum <- ifelse(!DT.tax$Phylum %in% ColPhylum, "Other", ifelse(DT.tax$Phylum == "Ascomycota", "Ascomycota", "Basidiomycota"))
	#enrichment test
	m = as(otu_table(physeq.filter), "matrix") + 1L
	t = data.frame(as(tax_table(physeq.filter), "matrix"))
	T=AnnotatedDataFrame(t)
	s = as(sample_data(physeq.filter), "data.frame")
	S =AnnotatedDataFrame(s)
	obj = newMRexperiment(m,phenoData=S,featureData=T) 
	p=cumNormStatFast(obj)
	objTrim=cumNorm(obj, p=p)
	Geno = pData(obj)$Geno
	settings = zigControl(maxit = 20, verbose = TRUE)
	dsg1=model.matrix(~0+Geno, data =s)
	res1 = fitZig(obj = objTrim, mod = dsg1, control = settings)
	zigFit1 = res1$fit
	finalMod1 = res1$fit$design
	#Tae vs Tbo
	c.mat1 = makeContrasts ( GenoTae	-	GenoTbo, levels = finalMod1)
	fit_c.mat1 = contrasts.fit(zigFit1, c.mat1)
	fit1 = eBayes(fit_c.mat1)
	DT_1.fc=data.table(fit1$coefficients, keep.rownames=T, key="rn")
	DT_1.p=data.table(fit1$p.value, keep.rownames=T, key="rn")
	setnames(DT_1.fc, "GenoTae - GenoTbo", "fc.TaeVsTbo")
	setnames(DT_1.p, "GenoTae - GenoTbo", "p.TaeVsTbo")
	#Tae vs Tur
	c.mat2 = makeContrasts ( GenoTae	-	GenoTur, levels = finalMod1)
	fit_c.mat2 = contrasts.fit(zigFit1, c.mat2)
	fit2 = eBayes(fit_c.mat2)
	DT_2.fc=data.table(fit2$coefficients, keep.rownames=T, key="rn")
	DT_2.p=data.table(fit2$p.value, keep.rownames=T, key="rn")
	setnames(DT_2.fc, "GenoTae - GenoTur", "fc.TaeVsTur")
	setnames(DT_2.p, "GenoTae - GenoTur", "p.TaeVsTur")
	#Tbo vs Tur
	c.mat3 = makeContrasts ( GenoTbo	-	GenoTur, levels = finalMod1)
	fit_c.mat3 = contrasts.fit(zigFit1, c.mat3)
	fit3 = eBayes(fit_c.mat3)
	DT_3.fc=data.table(fit3$coefficients, keep.rownames=T, key="rn")
	DT_3.p=data.table(fit3$p.value, keep.rownames=T, key="rn")
	setnames(DT_3.fc, "GenoTbo - GenoTur", "fc.TboVsTur")
	setnames(DT_3.p, "GenoTbo - GenoTur", "p.TboVsTur")

	DT.m<-DT.tax[DT_1.fc,][DT_1.p,][DT_2.fc,][DT_2.p,][DT_3.fc,][DT_3.p,]
	DT.m$ptest <- ifelse(DT.m$p.TaeVsTbo < alpha | DT.m$p.TaeVsTur < alpha | DT.m$p.TboVsTur < alpha, "IOS", "IANS")

	DT.m$p.shape1 <- ifelse ( DT.m$fc.TaeVsTbo > 0, "Enriched", "Depleted")
	DT.m$p.shape2 <- ifelse ( DT.m$fc.TaeVsTur > 0, "Enriched", "Depleted")
	DT.m$p.shape3 <- ifelse ( DT.m$fc.TboVsTur > 0, "Enriched", "Depleted")
	res[[i]]<-DT.m
	Tae.1<-DT.m[which(DT.m$p.TaeVsTbo<.05 & DT.m$fc.TaeVsTbo>0),]$rn
	Tae.2<-DT.m[which(DT.m$p.TaeVsTur<.05 & DT.m$fc.TaeVsTur>0),]$rn
	V.Tae<- venn(list(Tae.1, Tae.2), zcolor="style")
	Tae[[i]]<-unique(c(Tae.1, Tae.2), fromLast = TRUE)
	Tbo.1<-DT.m[which(DT.m$p.TaeVsTbo<.05 & DT.m$fc.TaeVsTbo<0),]$rn
	Tbo.2<-DT.m[which(DT.m$p.TboVsTur<.05 & DT.m$fc.TboVsTur>0),]$rn
	V.Tbo<- venn(list(Tbo.1, Tbo.2), zcolor="style")
	Tbo[[i]]<-unique(c(Tbo.1, Tbo.2), fromLast = TRUE)
	Tur.1<-DT.m[which(DT.m$p.TaeVsTur<.05 & DT.m$fc.TaeVsTur<0),]$rn
	Tur.2<-DT.m[which(DT.m$p.TboVsTur<.05 & DT.m$fc.TboVsTur<0),]$rn
	V.Tur<- venn(list(Tur.1, Tur.2), zcolor="style")
	Tur[[i]]<-unique(c(Tur.1, Tur.2), fromLast = TRUE)
}
#preparing the plots
	venn(list(Tae$Leaf, Tae$Root, Tae$Rhizo), snames=c("Leaf", "Root", "Rhizo"), zcolor="style")
	grid.echo()
	panel_a<-grid.grab()
	venn(list(Tbo$Leaf, Tbo$Root, Tbo$Rhizo), snames=c("Leaf", "Root", "Rhizo"), zcolor="style")
	grid.echo()
	panel_b<-grid.grab()
	venn(list(Tur$Leaf, Tur$Root, Tur$Rhizo), snames=c("Leaf", "Root", "Rhizo"), zcolor="style")
	grid.echo()
	panel_c<-grid.grab()

#	pdf("Supplementary_Fig_13abc.pdf", useDingbats=FALSE, paper="A4")
	gridExtra::grid.arrange( panel_a, panel_b, panel_c, nrow=2, ncol=2)
#	dev.off()

#saving enriched OTUs in txt file for ploting using Graphlan
	DT.taxa=data.table(tax, keep.rownames=T, key="rn")
	Uniq.Tae<-unique(c(Tae$Leaf, Tae$Root, Tae$Rhizo), fromLast = TRUE)
	DT.uniq.Tae<-DT.taxa[DT.taxa$rn %in% Uniq.Tae,]
	DT.uniq.Tae$Leaf<- ifelse (DT.uniq.Tae$rn %in% Tae$Leaf, "1", "0")
	DT.uniq.Tae$Root<- ifelse (DT.uniq.Tae$rn %in% Tae$Root, "1", "0")
	DT.uniq.Tae$Rhizo<- ifelse (DT.uniq.Tae$rn %in% Tae$Rhizo, "1", "0")
#	write.table(DT.uniq.Tae, "OTUs_Unique_Tae.txt", sep="\t")
	Uniq.Tbo<-unique(c(Tbo$Leaf, Tbo$Root, Tbo$Rhizo), fromLast = TRUE)
	DT.uniq.Tbo<-DT.taxa[DT.taxa$rn %in% Uniq.Tbo,]
	DT.uniq.Tbo$Leaf<- ifelse (DT.uniq.Tbo$rn %in% Tbo$Leaf, "1", "0")
	DT.uniq.Tbo$Root<- ifelse (DT.uniq.Tbo$rn %in% Tbo$Root, "1", "0")
	DT.uniq.Tbo$Rhizo<- ifelse (DT.uniq.Tbo$rn %in% Tbo$Rhizo, "1", "0")
#	write.table(DT.uniq.Tbo, "OTUs_Unique_Tbo.txt", sep="\t")
	Uniq.Tur<-unique(c(Tur$Leaf, Tur$Root, Tur$Rhizo), fromLast = TRUE)
	DT.uniq.Tur<-DT.taxa[DT.taxa$rn %in% Uniq.Tur,]
	DT.uniq.Tur$Leaf<- ifelse (DT.uniq.Tur$rn %in% Tur$Leaf, "1", "0")
	DT.uniq.Tur$Root<- ifelse (DT.uniq.Tur$rn %in% Tur$Root, "1", "0")
	DT.uniq.Tur$Rhizo<- ifelse (DT.uniq.Tur$rn %in% Tur$Rhizo, "1", "0")
#	write.table(DT.uniq.Tur, "OTUs_Unique_Tur.txt", sep="\t")
