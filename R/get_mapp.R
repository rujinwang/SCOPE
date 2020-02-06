if (getRversion() >= "2.15.1") {
    utils::globalVariables(c("mapp_hg19", 
                            "mapp_hg38", "seqlevelsStyle<-"))
}
#' @title Compute mappability
#' @name get_mapp
#'
#' @description Compute mappability for each bin. Note that scDNA
#' sequencing is whole-genome amplification and the mappability
#' score is essential to determine variable binning method.
#' Mappability track for 100-mers on the GRCh37/hg19 human
#' reference genome from ENCODE is pre-saved. Compute the mean
#' of mappability scores that overlapped reads map to bins,
#' weighted by the width of mappability tracks on the genome
#' reference. Use liftOver utility to calculate mappability
#' for hg38, which is pre-saved as well.
#'
#' @param ref GRanges object returned from \code{get_bam_bed}
#' @param hgref reference genome. This should be either 'hg19' or 'hg38'.
#'
#' @return
#'   \item{mapp}{Vector of mappability for each bin/target}
#'
#' @examples
#' \dontrun{
#' library(WGSmapp)
#' library(BSgenome.Hsapiens.UCSC.hg38)
#' bamfolder <- system.file('extdata', package = 'WGSmapp')
#' bamFile <- list.files(bamfolder, pattern = '*.dedup.bam$')
#' bamdir <- file.path(bamfolder, bamFile)
#' sampname_raw <- sapply(strsplit(bamFile, '.', fixed = TRUE), '[', 1)
#' bambedObj <- get_bam_bed(bamdir = bamdir,
#'                             sampname = sampname_raw, 
#'                             hgref = "hg38")
#' bamdir <- bambedObj$bamdir
#' sampname_raw <- bambedObj$sampname
#' ref_raw <- bambedObj$ref
#' 
#' mapp <- get_mapp(ref_raw, hgref = "hg38")
#' }
#'
#' @author Rujin Wang \email{rujin@email.unc.edu}
#' @import utils
#' @importFrom GenomicRanges GRanges pintersect
#' @importFrom IRanges IRanges RangesList Views countOverlaps findOverlaps width
#' @importFrom GenomeInfoDb mapSeqlevels seqlevelsStyle seqnames
#' @importFrom S4Vectors queryHits subjectHits
#' @export
get_mapp <- function(ref, hgref = "hg19") {
    if(!hgref %in% c("hg19", "hg38")){
        stop("Reference genome should be either hg19 or hg38. ")
    }
    if(hgref == "hg19") {
        mapp_gref <- mapp_hg19
    }else {
        mapp_gref <- mapp_hg38
    }
    mapp <- rep(1, length(ref))
    seqlevelsStyle(ref) <- "UCSC"
    for (chr in as.character(unique(seqnames(ref)))) {
        message("Getting mappability for ", chr, sep = "")
        chr.index <- which(as.matrix(seqnames(ref)) == chr)
        ref.chr <- ref[which(as.character(seqnames(ref)) == chr)]
        mapp.chr <- rep(1, length(ref.chr))
        overlap <- as.matrix(findOverlaps(ref.chr, mapp_gref))
        for (i in unique(overlap[, 1])) {
            index.temp <- overlap[which(overlap[, 1] == i), 2]
            overlap.sub <- findOverlaps(ref.chr[i], mapp_gref[index.temp])
            overlap.intersect <- pintersect(ref.chr[i][queryHits(
                                overlap.sub)], mapp_gref[index.temp][
                                subjectHits(overlap.sub)])
            mapp.chr[i] <- sum((mapp_gref$score[index.temp]) * 
                                (width(overlap.intersect)))/sum(
                                width(overlap.intersect))
        }
        mapp[chr.index] <- mapp.chr
    }
    mapp
}

