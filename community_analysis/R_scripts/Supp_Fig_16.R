####################################################################
#R scripts written by M. Amine Hassani - hassani.medamine@gmail.com#
####################################################################

#required R packages
pkg=c("plyr", "ggplot2", "phyloseq", "data.table", "reshape2", "ape", "PMCMR", "reldist","vegan","bipartite","phangorn","metagenomeSeq")

lapply(pkg, library, character.only = TRUE)

rm(list = ls()[grep("*.*", ls())])

#required source scripts
source("/MWDW/source_scripts/MetacommunityDynamicsFctsOikos.R") #Tucker et al 2016
source("/MWDW/source_scripts/PANullDevFctsOikos.R") #Tucker et al 2016
set.seed(131)

color_palette<-c("#cc2a36","#166590","#11902d","#999999")
theme_new <- theme (
	panel.grid.major = element_blank(),
	panel.grid.minor = element_blank(),
	legend.position="none",
	axis.text.x = element_text(angle=90, vjust=1)
	)
Bnull<-list()
Bnull.out<-data.frame()
LIST<-list ( "Leaf", "Root", "Rhizo", "UnSoil" )

# uploading and preparing phyloseq objects
	mat=read.table("otu_table_merge_bacteria.txt", sep="\t", row.names=1, header=T)
	mat=as.matrix(mat)
	OTU=otu_table(mat, taxa_are_rows=T) 
	tax=read.table("taxonomy_bacteria.txt", sep="\t", row.names=1, header=1)
	tax=as.matrix(tax)
	TAXA=tax_table(tax)
	sd=read.table("sample_data_merge.txt", sep="\t", row.names=1, header=1)
	SD=sample_data(sd)
	physeq=phyloseq(OTU, TAXA, SD)
	physeq_B=subset_taxa( physeq, Kingdom  == "Bacteria" )

###############################################
#This script took 07:35:46 to run in HPCluster#
###############################################

for( i in LIST ) 
{
	print(paste0("Computing null deviation for ", i, " samples " ))
#	print( paste0("subset ", i, " samples" ))
	SUB=c( i )
	physeq_SUB=subset_samples( physeq_B, Comp %in% SUB )
	comm=data.frame(otu_table(physeq_SUB))

	map=sd[colnames(comm),]
	rdp=tax[rownames(comm),]

	list.OTUs=row.names(comm)[rowSums(comm)>1]
        COMM=comm[rownames(comm) %in% list.OTUs, ]

#	print("preparing the data and setting parameters ... ")
	comm.t=t(COMM)
	bbs.sp.site <- comm.t
	patches=nrow(bbs.sp.site)
	rand <- 1000

	null.alphas <- matrix(NA, ncol(comm.t), rand)
	null.alpha <- matrix(NA, ncol(comm.t), rand)
	expected_beta <- matrix(NA, 1, rand)
	null.gamma <- matrix(NA, 1, rand)
	null.alpha.comp <- numeric()
	bucket_bray_res <- matrix(NA, patches, rand)

	bbs.sp.site = ceiling(bbs.sp.site/max(bbs.sp.site)) 
	mean.alpha = sum(bbs.sp.site)/nrow(bbs.sp.site) #mean.alpha
	gamma <- ncol(bbs.sp.site) #gamma
	obs_beta <- 1-mean.alpha/gamma
	obs_beta_all <- 1-rowSums(bbs.sp.site)/gamma

#	print (" Generating null patches ...")
	for (randomize in 1:rand) 
	{  
		null.dist = comm.t
		for (species in 1:ncol(null.dist)) 
		{
			tot.abund = sum(null.dist[,species])
			null.dist[,species] = 0
			for (individual in 1:tot.abund) 
			{
			sampled.site = sample(c(1:nrow(bbs.sp.site)), 1)
			null.dist[sampled.site, species] = null.dist[sampled.site, species] + 1
			}
		}
#		print ("Calculating null deviation for null patches ...")
		null.alphas[,randomize] <- apply(null.dist, 2, function(x){sum(ifelse(x > 0, 1, 0))})
		null.gamma[1, randomize] <- sum(ifelse(rowSums(null.dist)>0, 1, 0))
		expected_beta[1, randomize] <- 1 - mean(null.alphas[,randomize]/null.gamma[,randomize])
		null.alpha <- mean(null.alphas[,randomize])
		null.alpha.comp <- c(null.alpha.comp, null.alpha)
  
		bucket_bray <- as.matrix(vegdist(null.dist, "bray"))
		diag(bucket_bray) <- NA
		bucket_bray_res[,randomize] <- apply(bucket_bray, 2, FUN="mean", na.rm=TRUE)
	}
	## Computing beta-diversity for observed communities
	beta_comm_abund <- vegdist(comm.t, "bray")
	res_beta_comm_abund <- as.matrix(as.dist(beta_comm_abund))
	diag(res_beta_comm_abund) <- NA
	# output beta diversity (Bray)
	beta_div_abund_stoch <- apply(res_beta_comm_abund, 2, FUN="mean", na.rm=TRUE)
	# output abundance beta-null deviation
	bray_abund_null_dev <- beta_div_abund_stoch - mean(bucket_bray_res)

	betanull.out=data.frame(I(beta_div_abund_stoch),I(bray_abund_null_dev),stringsAsFactors=FALSE)
	colnames(betanull.out)=c("BRAY_beta_div_abund_stoch", "BRAY_AbundanceNullDeviation")

	Bnull[[i]]<-betanull.out
	Bnull.out<-rbind(Bnull.out, betanull.out)

#	print( paste0("Computing null deviation for ", i, " samples done" ))
	print("##############################################################")
}
###saving the data as txt
print("Computing BC null deviations completed")
#write.table(Bnull.out, "bnullout_bacteria.txt", quote=FALSE, sep="\t")
###uploading data - Bacteria 
	Bnull=read.table("bnullout_bacteria.txt", sep="\t", row.names=1, header=T)
	DT.B=data.table(Bnull, keep.rownames=T, key="rn")
	setnames(DT.B, c("BRAY_beta_div_abund_stoch", "BRAY_AbundanceNullDeviation"), c("BCB","NullB"))	
###uploading data - Fungi
	Fnull=read.table("bnullout_fungi.txt", sep="\t", row.names=1, header=T)
	DT.F=data.table(Fnull, keep.rownames=T, key="rn")
	setnames(DT.F, c("BRAY_beta_div_abund_stoch", "BRAY_AbundanceNullDeviation"), c("BCF","NullF"))
###preparing the plot
	sd=read.table("sample_data_merge.txt", sep="\t", row.names=1, header=1)
	DT.sd=data.table(sd, keep.rownames=T, key="rn")
	DT=DT.sd[DT.B][DT.F]
	DT<-DT[,-c("BCB","BCF")]
	DT.m<-melt(DT, measure.vars =c("NullB","NullF"))
	DT.m$variable<-as.factor(DT.m$variable)
	DT.m$VAR<-paste0(DT.m$Concat,DT.m$variable)
	DT.m$kingdom<-ifelse(DT.m$variable=="NullB","B","F")
	pp3=ggplot(DT.m, aes(x=VAR, y=as.numeric(value), color=kingdom))
	p3=pp3+geom_smooth(aes(group = kingdom))+ theme_bw(base_size=10)+theme_new +
	geom_point(size=1, color="Black")

#	pdf("Supplementary_Fig_16.pdf",paper="A4" ,useDingbats=FALSE)
	gridExtra::grid.arrange(p3)
#	dev.off()

