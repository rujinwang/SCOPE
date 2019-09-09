#' @title Get bam file directories, sample names, and whole genomic bins
#'
#' @description Get bam file directories, sample names, and whole genomic
#' bins from .bed file
#'
#' @param bamdir vector of the directory of a bam file. Should be in the same
#'  order as sample names in \code{sampname}.
#' @param bedFile path to the whole genome bed file
#' @param sampname vector of sample names. Should be in the same order as bam
#'  directories in \code{bamdir}.
#'
#' @return A list with components
#'     \item{bamdir}{A vector of bam directories}
#'     \item{sampname}{A vector of sample names}
#'     \item{ref}{A GRanges object specifying whole genomic bin positions}
#'
#' @examples
#' library(WGSmapp)
#' bedFile <- system.file('extdata', 'scWGA500kbsort.bed', package = 'SCOPE')
#' bamfolder <- system.file('extdata', package = 'WGSmapp')
#' bamFile <- list.files(bamfolder, pattern = '*.dedup.bam$')
#' bamdir <- file.path(bamfolder, bamFile)
#' sampname_raw = sapply(strsplit(bamFile, '.', fixed = TRUE), '[', 1)
#' bambedObj <- getbambed_scope(bamdir = bamdir, bedFile = bedFile,
#'                             sampname = sampname_raw)
#' bamdir <- bambedObj$bamdir
#' sampname_raw <- bambedObj$sampname
#' ref_raw <- bambedObj$ref
#'
#' @author Rujin Wang \email{rujin@email.unc.edu}
#' @import utils
#' @import GenomeInfoDb
#' @importFrom GenomicRanges GRanges
#' @importFrom IRanges IRanges
#' @export
getbambed_scope = function(bamdir, bedFile, sampname) {
    if (!file.exists(bedFile))
        stop("Please check the bed file directory provided.
            File could not be \nfound!")
    exomtarg <- read.table(bedFile, sep = "\t")
    ref <- GRanges(seqnames = exomtarg[, 1], ranges =
        IRanges(start = exomtarg[, 2], end = exomtarg[, 3]))
    if (!any(grepl("chr", seqlevels(ref)))) {
        seqlevels(ref) = paste(c(seq_len(22), "X", "Y"), sep = "")
        ref <- sort(ref)
    } else {
        seqlevels(ref) = paste("chr", c(seq_len(22), "X", "Y"), sep = "")
        ref <- sort(ref)
    }
    list(bamdir = bamdir, sampname = sampname, ref = ref)
}

