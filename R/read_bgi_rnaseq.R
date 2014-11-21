
# @importFrom colbycol cbc.read.table
bgi_readBasicExpressionData <- function(filename='all.gene.rpkm.xls') {
  require(colbycol)
  file <- filename
  if (grepl('gz$',filename)) {
  	file <- gzfile(filename)
  	dats <- read.table(file,sep="\t",header=T)
  } else {
	dats <- read.table( file,sep="\t",header=TRUE,nrows=1 )
	maxcol <- match('Symbol',names(dats))
	dats <- as.data.frame ( colbycol::cbc.read.table(file,sep='\t',header=T,just.read=1:maxcol) )
  }
  m <- as.matrix(dats)
  m[m=='-'] <- '0'
  dats <- as.data.frame(m)
  for( reads_col in grep("_Uniq",names(dats),value=T) ) {
  	dats[[reads_col]] <- as.numeric( dats[[ reads_col ]] )
  }
  names(dats) <- gsub("_Uniq.*","",names(dats))
  rownames(dats) <- dats$GeneID
  dats
}

get_reads_for_design <- function(all.reads,...) {
	design<-list(...)
	design.reads <- all.reads[,unlist(design)]
	attributes(design.reads)$length <- all.reads$Length
	attributes(design.reads)$conditions <- unlist( sapply(names(design),function(cond) { rep(cond, length(design[[cond]])) }) )
	attributes(design.reads)$design <-  design
	design.reads
}

rnaseq.prepareDifferential.DESeq <- function(all.reads,filtered=T,...) {
	getBiocLiteLib('DESeq')
	getBiocLiteLib('edgeR')

	require(edgeR)
	require(DESeq)
	design_reads <- get_reads_for_design(all.reads,...)
	y <- edgeR::DGEList(counts=design_reads,group=attributes(design_reads)$conditions,genes=data.frame(Length=attributes(design_reads)$length))
	if (filtered) {
		keep <- rowSums(edgeR::cpm(y)>1) >= min ( sapply(attributes(design_reads)$design, length ) )
	} else {
		keep <- rowSums(edgeR::cpm(y)>1) >= 0
	}
	y <- y[keep,]

	cds = newCountDataSet( y$counts, factor( attributes(design_reads)$conditions ) )
	cds = estimateSizeFactors( cds )
	cds = estimateDispersions( cds )
	cds
}

rnaseq.prepareDifferential.EdgeR <- function(all.reads,filtered=T,...) {
	getBiocLiteLib('edgeR')
	require(edgeR)
	design_reads <- get_reads_for_design(all.reads,...)
	y <- edgeR::DGEList(counts=design_reads,group=attributes(design_reads)$conditions,genes=data.frame(Length=attributes(design_reads)$length))
	if (filtered) {
		keep <- rowSums(edgeR::cpm(y)>1) >= min ( sapply(attributes(design_reads)$design, length ) )
	} else {
		keep <- rowSums(edgeR::cpm(y)>1) >= 0
	}
	y <- y[keep,]
	y$samples$lib.size <- colSums(y$counts)
	y <- edgeR::calcNormFactors(y)
	y <- edgeR::estimateCommonDisp(y,verbose=FALSE)
	y <- edgeR::estimateTagwiseDisp(y)
	y
}

rnaseq.getDifferentialGenes.DESeq <- function(deseq,up=T,adj.pval=0.1,taxonomy=9606) {
	conds = as.character(unique(conditions(deseq)))
	res = nbinomTest( deseq, conds[1], conds[2] )
	res = res[ res$padj <= adj.pval, ]
	names(res)[1] <- 'geneid'
	metadata<- convertEntrezIds(taxonomy,res$geneid)
	retvals <- subset( merge(res,metadata,all=T,by='geneid'),select=c('geneid','log2FoldChange','pval','uniprot','genename'))
	names(retvals) <- c('geneid','logFC','PValue','uniprot','genename')
	retvals
}

rnaseq.getDifferentialGenes.EdgeR <- function(edgeR,up=T,adj.pval=0.05,taxonomy=9606) {
	et <- edgeR::exactTest(edgeR,pair=as.character(unique(edgeR$samples$group)))
	if (up) {
		diffs <- et$table[which(decideTestsDGE(et,p=adj.pval,adjust="BH") > 0),]
	} else {
		diffs <- et$table[which(decideTestsDGE(et,p=adj.pval,adjust="BH") < 0),]
	}
	diffs$geneid <- rownames(diffs)
	metadata<- convertEntrezIds(taxonomy,rownames(diffs))
	subset(merge(diffs,metadata,all=T,by='geneid'),select=c('geneid','logFC','PValue','uniprot','genename'))
}

#' Get the biological coefficient of variation for a set of reads
#'
#' @param	all.reads	Reads as read in by \code{\link{rnaseq.readExpression}}
#' @param	filtered	Toggle if there should be a filtering for minimum number of reads (default TRUE)
#' @return	The biological coefficient of variation
#' @export
rnaseq.edgeR.getBCV <- function(all.reads,filtered=T,...) {
	dge_analysis <- rnaseq.prepareDifferential.EdgeR(all.reads,filtered,...)
	sqrt(dge_analysis$common.dispersion)
}

#' @export
rnaseq.edgeR.getDifferential <- function(all.reads,filtered=T,pval=0.05,taxonomy=9606,...) {
	dge_analysis <- rnaseq.prepareDifferential.EdgeR(all.reads,filtered,...)
	rbind( rnaseq.getDifferentialGenes.EdgeR(dge_analysis,up=T,pval,taxonomy), rnaseq.getDifferentialGenes.EdgeR(dge_analysis,up=F,pval,taxonomy))
}


#' @export
rnaseq.DESeq.getDifferential <- function(all.reads,filtered=T,pval=0.1,taxonomy=9606,...) {
	deseq_analysis <- rnaseq.prepareDifferential.DESeq(all.reads,filtered,...)
	rbind( rnaseq.getDifferentialGenes.DESeq(deseq_analysis,up=T,pval,taxonomy), rnaseq.getDifferentialGenes.DESeq(deseq_analysis,up=F,pval,taxonomy))
}

max_values <- function(df,fields,distance) {
	expand.grid(fields,fields)
}

#' Perform a differential analysis on a set of expression values from RNAseq experiments
#'
#' @param	all.reads	Reads as read in by \code{\link{rnaseq.readExpression}}
#' @param	filtered	Toggle if there should be a filtering for minimum number of reads (default TRUE)
#' @param	clonal.variation.threshold	Threshold to drop genes if there is a log2 variation between clones. Set to NA to remove this filtering
#' @param	pval.edger	EdgeR adjusted p-value cutoff for differential expression
#' @param	pval.deseq	DESeq adjusted p-value cutoff for differential expression
#' @return	Data.frame with geneid logFC PValue uniprot and genename
#' @export
rnaseq.getDifferential <- function(all.reads,filtered=T,clonal.variation.threshold=1,pval.edger=0.05,pval.deseq=0.1,taxonomy=9606,...) {
	dge_analysis <- rnaseq.prepareDifferential.EdgeR(all.reads,filtered,...)
	deseq_analysis <- rnaseq.prepareDifferential.DESeq(all.reads,filtered,...)
	logrpkms <- rpkm(dge_analysis,as.numeric(dge_analysis$genes$Length),log=T,prior.count=2)
	if ( !is.na(clonal.variation.threshold)) {
		samples <- as.data.frame(dge_analysis$samples)
		samples$sample <- rownames(samples)
		groups <- plyr::dlply( samples, plyr::.(group), function(df) { df$sample } )
		for(group in names(groups)) {
			combinations <- expand.grid(groups[[group]],groups[[group]])
			for (rownum in 1:nrow(combinations)) {
				logrpkms <- logrpkms[which( abs( logrpkms[,as.character(combinations[rownum,1]) ] - logrpkms[,as.character(combinations[rownum,2])] ) <= clonal.variation.threshold),]
			}
		}
	}
	rpkm_wanted <- rownames(logrpkms)
	edger.up <- rnaseq.getDifferentialGenes.EdgeR(dge_analysis,up=T,pval.edger,taxonomy)
	edger.down <- rnaseq.getDifferentialGenes.EdgeR(dge_analysis,up=F,pval.edger,taxonomy)
	deseq.up <- rnaseq.getDifferentialGenes.DESeq(deseq_analysis,up=T,pval.deseq,taxonomy)
	deseq.down <- rnaseq.getDifferentialGenes.DESeq(deseq_analysis,up=F,pval.deseq,taxonomy)

	wanted <- rbind ( subset(edger.up,geneid %in% intersect(rpkm_wanted, deseq.up$geneid) ) , subset(edger.down,geneid %in% intersect(rpkm_wanted, deseq.down$geneid ) ) )
	wanted
}

#' Read an expression file obtained from BGI
#' @return Data frame with expression data
#' @export
rnaseq.readExpression <- function(filename) {
	bgi_readBasicExpressionData(filename)
}