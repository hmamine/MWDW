####################################################################
#R scripts written by M. Amine Hassani - hassani.medamine@gmail.com#
####################################################################
#required packages
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
	panel.background = element_rect(fill="white", colour ="black", size = 1)
	)
Palette.phylum <-c("#4e8e8e","#8e8e4e","#8e4e8e")
#Palette.phylum <-c("#8e8e4e","#8e4e8e")
Shape.1=c(21,19)

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
#subset to Fungi	
	physeq_F=subset_taxa( physeq, Kingdom == "Fungi" )
#subset to root samples
	physeq.sub=subset_samples( physeq_F, Comp == "Root" )
#computing relative abundance
	physeq.sub.ra = transform_sample_counts(physeq.sub, function(x) x/sum(x))
#filtering OTUs by occurance
	otu.i=otu_table(physeq.sub)
	taxa.i=tax_table(physeq.sub)
	sd.i=sample_data(physeq.sub)
	flt=filterTaxonMatrix(otu.i,minocc=2, keepSum = TRUE, return.filtered.indices = TRUE)
	otu.i.filtered=flt$mat
	taxa.i.filtered=taxa.i[setdiff(1:nrow(taxa.i),flt$filtered.indices),]
	dummyTaxonomy=c("k__dummy","p__","c__","o__","f__","g__","s__","sc__")
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
#Tbo vs Tur
	c.mat2 = makeContrasts ( GenoTae	-	GenoTur, levels = finalMod1)
	fit_c.mat2 = contrasts.fit(zigFit1, c.mat2)
	fit2 = eBayes(fit_c.mat2)
	DT_2.fc=data.table(fit2$coefficients, keep.rownames=T, key="rn")
	DT_2.p=data.table(fit2$p.value, keep.rownames=T, key="rn")
	setnames(DT_2.fc, "GenoTae - GenoTur", "fc.TaeVsTur")
	setnames(DT_2.p, "GenoTae - GenoTur", "p.TaeVsTur")
#Tae vs Tur
	c.mat3 = makeContrasts ( GenoTbo	-	GenoTur, levels = finalMod1)
	fit_c.mat3 = contrasts.fit(zigFit1, c.mat3)
	fit3 = eBayes(fit_c.mat3)
	DT_3.fc=data.table(fit3$coefficients, keep.rownames=T, key="rn")
	DT_3.p=data.table(fit3$p.value, keep.rownames=T, key="rn")
	setnames(DT_3.fc, "GenoTbo - GenoTur", "fc.TboVsTur")
	setnames(DT_3.p, "GenoTbo - GenoTur", "p.TboVsTur")
#preparing plots
	DT.m<-DT.tax[DT_1.fc,][DT_1.p,][DT_2.fc,][DT_2.p,][DT_3.fc,][DT_3.p,]
	DT.m$ptest <- ifelse(DT.m$p.TaeVsTbo < alpha | DT.m$p.TaeVsTur < alpha | DT.m$p.TboVsTur < alpha, "IOS", "IANS")
	DT.m$p.shape1 <- ifelse ( DT.m$fc.TaeVsTbo > 0, "Enriched", "Depleted")
	DT.m$p.shape2 <- ifelse ( DT.m$fc.TaeVsTur > 0, "Enriched", "Depleted")
	DT.m$p.shape3 <- ifelse ( DT.m$fc.TboVsTur > 0, "Enriched", "Depleted")

##volcano plot
###Tae vs Tbo - Panel a
	pp4=ggplot(DT.m, aes(x= fc.TaeVsTbo, y=-log2(p.TaeVsTbo), color = ColPhylum))
	panel_a=pp4+geom_jitter( size=1.5 )+theme_bw( )+
	scale_x_continuous( breaks=c( -10,-8,-6,-4,-2,0,2,4,6,8,10 ) ) +
	geom_hline(yintercept = -log2(0.05), linetype="longdash")+
	geom_vline(xintercept = 0, linetype="longdash")+
	scale_color_manual(values=Palette.phylum)+
	ggtitle("Tae-Tbo") +
	theme_new 
###Tae vs Tur - Panel b
	pp5=ggplot(DT.m, aes(x= fc.TaeVsTur, y=-log2(p.TaeVsTur), color = ColPhylum))
	panel_b=pp5+geom_jitter( size=1.5 )+theme_bw( )+
	scale_x_continuous( breaks=c( -14,-12,-10,-8,-6,-4,-2,0,2,4,6,8,10,12,14 ) ) +
	geom_hline(yintercept = -log2(0.05), linetype="longdash")+
	geom_vline(xintercept = 0, linetype="longdash")+
	scale_color_manual(values=Palette.phylum)+
	ggtitle("Tae-Tur") +
	theme_new
###Tbo vs Tur - Panel c
	pp6=ggplot(DT.m, aes(x= fc.TboVsTur, y=-log2(p.TboVsTur), color = ColPhylum))
	panel_c=pp6+geom_jitter( size=1.5 )+theme_bw( ) +
	scale_x_continuous( breaks=c( -10,-8,-6,-4,-2,0,2,4,6,8,10 ) ) +
	geom_hline(yintercept = -log2(0.05), linetype="longdash")+
	geom_vline(xintercept = 0, linetype="longdash")+
	scale_color_manual(values=Palette.phylum)+
	ggtitle("Tbo-Tur") +	
	theme_new 
#pdf("volcano_plot.pdf", useDingbats=FALSE, paper="A4")
gridExtra::grid.arrange(panel_a, panel_b, panel_c, nrow=2, ncol=2)
#dev.off()

##venn plots - Panel d
###Tae
	Tae.1<-DT.m[which(DT.m$p.TaeVsTbo<.05 & DT.m$fc.TaeVsTbo>0),]$rn
	Tae.2<-DT.m[which(DT.m$p.TaeVsTur<.05 & DT.m$fc.TaeVsTur>0),]$rn
	V.Tae<- venn(list(Tae.1, Tae.2), zcolor="style")
	grid.echo()
	TAE<-grid.grab()
###Tbo
	Tbo.1<-DT.m[which(DT.m$p.TaeVsTbo<.05 & DT.m$fc.TaeVsTbo<0),]$rn
	Tbo.2<-DT.m[which(DT.m$p.TboVsTur<.05 & DT.m$fc.TboVsTur>0),]$rn
	V.Tbo<- venn(list(Tbo.1, Tbo.2), zcolor="style")
	grid.echo()
	TBO<-grid.grab()
###Tur
	Tur.1<-DT.m[which(DT.m$p.TaeVsTur<.05 & DT.m$fc.TaeVsTur<0),]$rn
	Tur.2<-DT.m[which(DT.m$p.TboVsTur<.05 & DT.m$fc.TboVsTur<0),]$rn
	V.Tur<- venn(list(Tur.1, Tur.2), zcolor="style")
	grid.echo()
	TUR<-grid.grab()
#pdf("Supplementary_Fig_11.pdf", useDingbats=FALSE, paper="A4")
gridExtra::grid.arrange( TAE, TBO, TUR, nrow=2, ncol=2)
#dev.off()

