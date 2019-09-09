#' @title Get read coverage from single-cell DNA sequencing
#'
#' @description Get read coverage for each genomic bin across all single
#'     cells from scDNA-seq.
#'
#' @param bambedObj object returned from \code{getbambed_scope}
#' @param mapqthres mapping quality threshold of reads
#' @param mask.ref a GRanges object indicating bad regions/bins,
#' such as segmental duplication regions and gaps near
#' telomeres/centromeres, which need to be masked prior to
#' getting coverage
#' @param seq the sequencing method to be used. This should be either
#' 'paired-end' or 'single-end'
#'
#' @return
#'   \item{Y}{Read depth matrix}
#'
#' @examples
#' library(WGSmapp)
#' bedFile <- system.file('extdata', 'scWGA500kbsort.bed',
#'                         package = 'SCOPE')
#' bamfolder <- system.file('extdata', package = 'WGSmapp')
#' bamFile <- list.files(bamfolder, pattern = '*.dedup.bam$')
#' bamdir <- file.path(bamfolder, bamFile)
#' sampname_raw = sapply(strsplit(bamFile, '.', fixed = TRUE), '[', 1)
#' bambedObj <- getbambed_scope(bamdir = bamdir,
#'                             bedFile = bedFile,
#'                             sampname = sampname_raw)
#' # Get segmental duplication regions
#' seg.dup = read.table(system.file('extdata',
#'                     'GRCh37GenomicSuperDup.tab',
#'                     package = 'WGSmapp'),
#'                     head = TRUE)
#' seg.dup = seg.dup[!is.na(match(seg.dup[,1],
#'                     paste('chr', c(1:22, 'X', 'Y'), sep = ''))),]
#' seg.dup = GRanges(seqnames = seg.dup[,1],
#'                     ranges = IRanges(start=seg.dup[,2], end = seg.dup[,3]))
#' # Get hg19 gaps
#' gaps = read.table(system.file('extdata', 'hg19gaps.txt',
#'                     package = 'WGSmapp'), head = TRUE)
#' gaps = gaps[!is.na(match(gaps[,2],
#'                     paste('chr', c(1:22, 'X', 'Y'), sep=''))),]
#' gaps = GRanges(seqnames = gaps[,2],
#'                 ranges = IRanges(start = gaps[,3], end = gaps[,4]))
#' # Generate mask region
#' mask.ref = sort(c(seg.dup, gaps))
#'
#' # Getting raw read depth
#' coverageObj <- getcoverage.scDNA(bambedObj,
#'                                 mapqthres = 40,
#'                                 mask.ref,
#'                                 seq='paired-end')
#' Y_raw <- coverageObj$Y
#'
#' @author Rujin Wang \email{rujin@email.unc.edu}
#' @import Rsamtools
#' @importFrom GenomicRanges GRanges
#' @importFrom IRanges IRanges RangesList Views countOverlaps
#' @importFrom GenomeInfoDb seqnames
#' @export
getcoverage.scDNA = function(bambedObj, mapqthres, mask.ref, seq) {
    ref <- bambedObj$ref
    bamdir <- bambedObj$bamdir
    sampname <- bambedObj$sampname

    Y = matrix(nrow = length(ref), ncol = length(sampname))
    rownames(Y) = paste(seqnames(ref), ":", start(ref), "-", end(ref), sep = "")
    colnames(Y) = sampname
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
                bam.ref = GRanges(seqnames = bam$rname, ranges =
                    IRanges(start = bam[["pos"]], width = bam[["qwidth"]]))
            } else {
                bam.ref = GRanges(seqnames = paste0("chr", bam$rname), ranges =
                    IRanges(start = bam[["pos"]], width = bam[["qwidth"]]))
            }
            bam.ref = bam.ref[bam$mapq >= mapqthres]
            bam.ref = suppressWarnings(bam.ref[countOverlaps(bam.ref,
                mask.ref) == 0])
            Y[, i] <- countOverlaps(ref, bam.ref)
        }
    }
    list(Y = Y)
}
