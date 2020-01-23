if (getRversion() >= "2.15.1") {
    utils::globalVariables(c("BSgenome.Hsapiens.UCSC.hg19",
                            "BSgenome.Hsapiens.UCSC.hg38", 
                            "seqlevels<-", "seqlevels"))
}
#' @title Get bam file directories, sample names, and whole genomic bins
#'
#' @description Get bam file directories, sample names, and whole genomic
#' bins from .bed file
#'
#' @param bamdir vector of the directory of a bam file. Should be in the same
#'  order as sample names in \code{sampname}.
#' @param sampname vector of sample names. Should be in the same order as bam
#'  directories in \code{bamdir}.
#' @param hgref reference genome. This should be either 'hg19' or 'hg38'.
#' @param resolution numeric value of fixed bin-length. Default is \code{500}.
#'  Unit is "kb". 
#'
#' @return A list with components
#'     \item{bamdir}{A vector of bam directories}
#'     \item{sampname}{A vector of sample names}
#'     \item{ref}{A GRanges object specifying whole genomic bin positions}
#'
#' @examples
#' library(WGSmapp)
#' library(BSgenome.Hsapiens.UCSC.hg38)
#' bamfolder <- system.file('extdata', package = 'WGSmapp')
#' bamFile <- list.files(bamfolder, pattern = '*.dedup.bam$')
#' bamdir <- file.path(bamfolder, bamFile)
#' sampname_raw <- sapply(strsplit(bamFile, '.', fixed = TRUE), '[', 1)
#' bambedObj <- get_bam_bed(bamdir = bamdir, sampname = sampname_raw, 
#'                         hgref = "hg38")
#' bamdir <- bambedObj$bamdir
#' sampname_raw <- bambedObj$sampname
#' ref_raw <- bambedObj$ref
#'
#' @author Rujin Wang \email{rujin@email.unc.edu}
#' @import utils
#' @import BSgenome.Hsapiens.UCSC.hg19
#' @importFrom GenomicRanges GRanges tileGenome
#' @importFrom IRanges IRanges
#' @importFrom GenomeInfoDb seqnames seqinfo seqlevels
#' @export
get_bam_bed <- function(bamdir, sampname, hgref = "hg19", resolution = 500) {
    if(!hgref %in% c("hg19", "hg38")){
        stop("Reference genome should be either hg19 or hg38. ")
    }
    if(hgref == "hg19") {
        genome <- BSgenome.Hsapiens.UCSC.hg19
    }else if(hgref == "hg38") {
        genome <- BSgenome.Hsapiens.UCSC.hg38
    }
    if(resolution <= 0){
        stop("Invalid fixed bin length. ")
    }
    bins <- tileGenome(seqinfo(genome), 
                    tilewidth = resolution * 1000, 
                    cut.last.tile.in.chrom = TRUE)
    ref <- bins[which(as.character(seqnames(bins)) %in% paste0("chr", 
                                                        seq_len(22)))]
    if (!any(grepl("chr", seqlevels(ref)))) {
        seqlevels(ref) <- paste(c(seq_len(22), "X", "Y"), sep = "")
        ref <- sort(ref)
    } else {
        seqlevels(ref) <- paste("chr", c(seq_len(22), "X", "Y"), sep = "")
        ref <- sort(ref)
    }
    list(bamdir = bamdir, sampname = sampname, ref = ref)
}

