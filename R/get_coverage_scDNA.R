#' @title Get read coverage from single-cell DNA sequencing
#'
#' @description Get read coverage for each genomic bin across all single
#'     cells from scDNA-seq. Blacklist regions, such as segmental duplication 
#'     regions and gaps near telomeres/centromeres will be masked prior to
#'     getting coverage. 
#'
#' @param bambedObj object returned from \code{get_bam_bed}
#' @param mapqthres mapping quality threshold of reads
#' @param seq the sequencing method to be used. This should be either
#' 'paired-end' or 'single-end'
#' @param hgref human reference genome. This should be either 'hg19' or 'hg38'.
#'
#' @return
#'   \item{Y}{Read depth matrix}
#'
#' @examples
#' library(WGSmapp)
#' bamfolder <- system.file('extdata', package = 'WGSmapp')
#' bamFile <- list.files(bamfolder, pattern = '*.dedup.bam$')
#' bamdir <- file.path(bamfolder, bamFile)
#' sampname_raw <- sapply(strsplit(bamFile, '.', fixed = TRUE), '[', 1)
#' bambedObj <- get_bam_bed(bamdir = bamdir,
#'                             sampname = sampname_raw)
#'
#' # Getting raw read depth
#' coverageObj <- get_coverage_scDNA(bambedObj,
#'                                 mapqthres = 40,
#'                                 seq = 'paired-end')
#' Y_raw <- coverageObj$Y
#'
#' @author Rujin Wang \email{rujin@email.unc.edu}
#' @import Rsamtools
#' @importFrom GenomicRanges GRanges
#' @importFrom IRanges IRanges RangesList Views countOverlaps
#' @importFrom GenomeInfoDb seqnames
#' @export
get_coverage_scDNA <- function(bambedObj, mapqthres, seq, hgref = "hg19") {
    if(!hgref %in% c("hg19", "hg38")){
        stop("Reference genome should be either hg19 or hg38!")
    } 
    ref <- bambedObj$ref
    bamdir <- bambedObj$bamdir
    sampname <- bambedObj$sampname

    if(hgref == "hg19"){
        # Get segmental duplication regions
        seg.dup <- read.table(system.file("extdata", 
                                        "GRCh37GenomicSuperDup.tab", 
                                        package = "WGSmapp"), header = TRUE)
        # Get hg19 gaps
        gaps <- read.table(system.file("extdata", "hg19gaps.txt", 
                                        package = "WGSmapp"), 
                                        header = TRUE)
    } else{
        # Get segmental duplication regions
        seg.dup <- read.table(system.file("extdata", 
                                        "GRCh38GenomicSuperDup.tab", 
                                        package = "WGSmapp"))
        # Get hg19 gaps
        gaps <- read.table(system.file("extdata", "hg38gaps.txt", 
                                        package = "WGSmapp"))
    }

    seg.dup <- seg.dup[!is.na(match(seg.dup[,1], 
                        paste('chr', c(seq_len(22), 'X', 'Y'), 
                        sep = ''))),]
    seg.dup <- GRanges(seqnames = seg.dup[,1], 
                        ranges = IRanges(start = seg.dup[,2], 
                                        end = seg.dup[,3]))
    gaps <- gaps[!is.na(match(gaps[,2], 
                        paste('chr', c(seq_len(22), 'X', 'Y'), 
                        sep = ''))),]
    gaps <- GRanges(seqnames = gaps[,2], 
                        ranges = IRanges(start = gaps[,3], 
                        end = gaps[,4]))
    # Generate mask region
    mask.ref <- sort(c(seg.dup, gaps))
    
    
    Y <- matrix(nrow = length(ref), ncol = length(sampname))
    rownames(Y) <- paste(seqnames(ref), ":", start(ref), "-", 
                        end(ref), sep = "")
    colnames(Y) <- sampname
    for (i in seq_len(length(sampname))) {
        bamurl <- bamdir[i]
        what <- c("rname", "pos", "mapq", "qwidth")
        if (seq == "paired-end") {
            flag <- scanBamFlag(isPaired = TRUE, isDuplicate = FALSE,
                isUnmappedQuery = FALSE, isNotPassingQualityControls = FALSE,
                isFirstMateRead = TRUE)
            param <- ScanBamParam(what = what, flag = flag)
            bam <- scanBam(bamurl, param = param)[[1]]
        } else if (seq == "single-end") {
            flag <- scanBamFlag(isPaired = FALSE, isDuplicate = FALSE,
                isUnmappedQuery = FALSE, isNotPassingQualityControls = FALSE)
            param <- ScanBamParam(what = what, flag = flag)
            bam <- scanBam(bamurl, param = param)[[1]]
        }
        message("Getting coverage for sample ", i, ": ",
            sampname[i], "...", sep = "")
        if (length(bam$rname) == 0) {
            Y[, i] <- 0  # Failed library preparation
        } else {
            if (any(grepl("chr", bam$rname) == TRUE)) {
                bam.ref <- GRanges(seqnames = bam$rname, ranges =
                    IRanges(start = bam[["pos"]], width = bam[["qwidth"]]))
            } else {
                bam.ref <- GRanges(seqnames = paste0("chr", bam$rname), ranges =
                    IRanges(start = bam[["pos"]], width = bam[["qwidth"]]))
            }
            bam.ref <- bam.ref[bam$mapq >= mapqthres]
            bam.ref <- suppressWarnings(bam.ref[countOverlaps(bam.ref,
                mask.ref) == 0])
            Y[, i] <- countOverlaps(ref, bam.ref)
        }
    }
    list(Y = Y)
}