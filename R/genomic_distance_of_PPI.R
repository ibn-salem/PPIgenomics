#'
#'
#' A script to analyse genomic distance distribution of gene pairs that encode 
#' for proteins that have direct protein-protein interactinos (PPI).
#' 
#'

require(stringr)        # for some string functionality
require(biomaRt)        # to retrieve human paralogs from Ensembl
require(BSgenome.Hsapiens.UCSC.hg19)
require(ggplot2)

# set some parameters:

# to download thies files, run the script data/download.sh
HIPPIE_SCORE_TH <- 0.72
HIPPIE_FILE <- "data/HIPPIE/hippie_current.txt"

TAD_FILE <- "data/Rao2014/GSE63525_IMR90_Arrowhead_domainlist.txt.bed"

outPrefix <- "resutls/PPI_genomics"

# create directory, if not exist
dir.create(dirname(outPrefix), recursive=TRUE, showWarnings = FALSE)


#' Add linear distance between genes.
#' 
#' Distance is measured from start of each region and reproted in kilobaes. If
#' the genes are on different chromosome, NA is reported.
#' 
#' @param genePair a \code{data.frames} where each row is a gene pair with the 
#'   first columns holding gnee IDs
#' @param tssGR a \code{\link{GRanges}} object with genes. The names should 
#'   match the gene ids in \code{genePairs}.
#' @return a \code{data.frame} with with the same data as \code{genePair} but 
#'   with an additional column \code{dist} holding the pairwise distances in kb.
addPairDistKb <- function(genePairs, tssGR){
  
  # get chromosomes of gene pairs
  chr1 <- as.character(seqnames(tssGR[genePairs[,1]]))
  chr2 <- as.character(seqnames(tssGR[genePairs[,2]]))
  
  sameChrom <-  chr1 == chr2
  s1 = start(tssGR[genePairs[,1]])
  s2 = start(tssGR[genePairs[,2]])
  
  # add a new column "dist" to the data.frame
  genePairs[, "dist"] = ifelse(sameChrom, abs(s2-s1)/1000, NA)
  return(genePairs)
}


#=======================================================================
# Analyse genomic distance distribution of PPI and non PPI gene pairs
#=======================================================================

#-------------------------------------------------------------------
# get tssGR for ENSG
#-------------------------------------------------------------------
seqInfo <- seqinfo(BSgenome.Hsapiens.UCSC.hg19)

ensemblGRCh37 <- useMart(host="grch37.ensembl.org", biomart="ENSEMBL_MART_ENSEMBL",dataset="hsapiens_gene_ensembl", verbose=FALSE)

geneAttributes = c("ensembl_gene_id", "hgnc_symbol", "chromosome_name", "start_position", "end_position", "strand", "status", "gene_biotype")
geneFilters="chromosome_name"
# read "normal" human chromosome names (without fixes and patches)
geneValues=c(1:22, "X", "Y")
allGenes = getBM(attributes=geneAttributes, mart=ensemblGRCh37, filters=geneFilters, values=geneValues)

# unique gene entry by ENSG ID symbol:
genes = allGenes[!duplicated(allGenes$ensembl_gene_id),]

# make GRanges object for all known prot coding genes

tssGR = GRanges(
        paste0("chr", genes$chromosome_name),
        IRanges(genes$start_position, genes$start_position),
        strand = ifelse(genes$strand == 1, '+', '-'), 
        names = genes$ensembl_gene_id, 
        genes[,c("hgnc_symbol", "status", "gene_biotype")],
        seqinfo=seqInfo
        )
names(tssGR) = genes$ensembl_gene_id
tssGR <- sort(tssGR)


#-------------------------------------------------------------------
# get mapping of entrez IDs to ENGS from ENSEMBL
#-------------------------------------------------------------------

entrezAttributes <- c("entrezgene", "ensembl_gene_id")
entrezFilters <- c("chromosome_name", "with_entrezgene")
entrezValues <- list("chromosome_name"=c(1:22, "X", "Y"), with_entrezgene=TRUE) 
entrezToEnsgDF = getBM(attributes=entrezAttributes, mart=ensemblGRCh37, filters=entrezFilters, values=entrezValues)

# take only unique entrez IDs 
entrezToEnsgDF <- entrezToEnsgDF[!duplicated(entrezToEnsgDF$entrezgene),]


#-----------------------------------------------------------------------
# Parse HIPPIE 
#-----------------------------------------------------------------------

hippieDF <- read.table(HIPPIE_FILE, header=FALSE, sep="\t", quote="")

# get index in mapping table for each entrez gene ID in HIPPIE
idxG1 <- match(as.character(hippieDF[,2]), entrezToEnsgDF$entrezgene)
idxG2 <- match(as.character(hippieDF[,4]), entrezToEnsgDF$entrezgene)

hippie <- data.frame(
	g1_ENSG = entrezToEnsgDF$ensembl_gene_id[idxG1],
	g2_ENSG = entrezToEnsgDF$ensembl_gene_id[idxG2],
	symbol1 = str_split_fixed(as.character(hippieDF[,1]), "_", 2)[,1],
	symbol2 = str_split_fixed(as.character(hippieDF[,3]), "_", 2)[,1],
	score = hippieDF[,5],
	stringsAsFactors=FALSE)

message("INFO: After parsing: ", nrow(hippie))

# filter out interactions that could not be mapped to ENSG
hippie <- hippie[!is.na(hippie[,1]) & !is.na(hippie[,2]),]
message("INFO: After ENSG mapping: ", nrow(hippie))

# filter out interaction bellow score threshold
hippie <- hippie[hippie$score >= HIPPIE_SCORE_TH,]


#-----------------------------------------------------------------------
# generate random interaction network
#-----------------------------------------------------------------------
randNet <- hippie[rep(1:nrow(hippie), 10) ,c("g1_ENSG", "g2_ENSG", "score")]
randNet[,2] <- sample(randNet[,2])


#-----------------------------------------------------------------------
# combine HIPPIE and random interactions
#-----------------------------------------------------------------------
pairsDF <- rbind(
	hippie[,c("g1_ENSG", "g2_ENSG", "score")],
	randNet
	)
pairsDF$group <- rep(c("PPI", "shuffled"), c(nrow(hippie), nrow(randNet)))

message("INFO: After filtering score >= ", HIPPIE_SCORE_TH, " : ", sum(pairsDF$group == "PPI"), " and shuffled: ",sum(pairsDF$group == "shuffled"))


#-----------------------------------------------------------------------
# Annotate gene pairs with genomic distance and filter for same chrom.
#-----------------------------------------------------------------------

# add distance 
pairsDF <- addPairDistKb(pairsDF, tssGR)

# filter for pairs on same chromosome (with dist != NA)
pairsDF <- pairsDF[!is.na(pairsDF$dist),]

message("INFO: After filtering out different chromosomes : ", 
        sum(pairsDF$group == "PPI"), " and shuffled: ", 
        sum(pairsDF$group == "shuffled"))


message("INFO: PPI pairs with dist==0: ", sum(pairsDF$group == "PPI" & pairsDF$dist == 0))
message("INFO: PPI pairs with same ID: ", sum(pairsDF$group == "PPI" & pairsDF[,1] == pairsDF[,2]))

# filter out pairs with same ID
pairsDF <- pairsDF[!pairsDF[,1] == pairsDF[,2],]
message("INFO: After filtering out homo-dimers (pairs with same ID): ", sum(pairsDF$group == "PPI"), " and shuffled: ",sum(pairsDF$group == "shuffled"))

pairsDF <- pairsDF[pairsDF$dist <= 1000,]
message("INFO: After filtering distance <= 1000kb: ", sum(pairsDF$group == "PPI"), " and shuffled: ",sum(pairsDF$group == "shuffled"))



#-----------------------------------------------------------------------
# annotate to be in same TAD
#-----------------------------------------------------------------------
getPairAsGR <- function(genePairs, tssGR){
  # get chromosomes of gene pairs
  chrom = seqnames(tssGR[genePairs[,1]])
  
  s1 = start(tssGR[genePairs[,1]])
  s2 = start(tssGR[genePairs[,2]])
  up = apply(cbind(s1, s2), 1, min)
  down = apply(cbind(s1, s2), 1, max)
  GR = GRanges(chrom, IRanges(up, down))
  # add gene IDs and other annotations
  mcols(GR) = genePairs
  return(GR)
}

#-----------------------------------------------------------------------
# add column to indicate that query lies within at least one subject object
#-----------------------------------------------------------------------
addWithinSubject <- function(query, subject, colName="inRegion"){
    mcols(query)[, colName] = countOverlaps(query, subject, type="within") >= 1
    return(query)
}

#-----------------------------------------------------------------------
# parse TADs from Rao et al. 
#-----------------------------------------------------------------------

# parse TADs from bed file
tadGR <- import(TAD_FILE, seqinfo=seqInfo)

# get gene-pair spanning regions as GRanges
pairsGR <- getPairAsGR(pairsDF, tssGR)

# check overlap of gnee-pair spanning region is within any TAD
pairsDF$inTAD <- countOverlaps(pairsGR, tadGR, type="within") >= 1
pairsDF$inTAD <- factor(pairsDF$inTAD, c(TRUE, FALSE), c("Same TAD", "Not same TAD"))



#===============================================================================
# plot geomic distance distribution
#===============================================================================

# compute p-value for distance difference between HIPPIE and shuffled
pVal <- wilcox.test(dist ~ group, data=pairsDF)$p.value


p <- ggplot(pairsDF, aes(dist, ..density.., fill=group, color=group)) + 
	geom_histogram(binwidth=50, alpha=.5, position="identity") + 
	labs(title=paste("p =", signif(pVal, 3)), x="Genomic distance [kb]")  + 
	theme_bw()	 
ggsave(p, file=paste0(outPrefix, ".hippie_genomic_distance.v03.hist.pdf"), w=7, h=3.5)

p <- p + facet_grid(inTAD~., margins=TRUE, scales="free_y")
ggsave(p, file=paste0(outPrefix, ".hippie_genomic_distance.v03.hist.byTAD.pdf"), w=7, h=7)


