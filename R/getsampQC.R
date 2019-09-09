#' @title Get QC metrics for single cells
#'
#' @description Perform QC step on single cells.
#'
#' @param bambedObj object returned from \code{getbambed_scope}
#'
#' @return
#'   \item{QCmetric}{A matrix containing total number/proportion of reads,
#'     total number/proportion of mapped reads, total number/proportion
#'     of mapped non-duplicate reads, and number/proportion of reads with
#'     mapping quality greater than 20}
#'
#' @examples
#' library(WGSmapp)
#' bedFile <- system.file('extdata', 'scWGA500kbsort.bed', package = 'SCOPE')
#' bamfolder <- system.file('extdata', package = 'WGSmapp')
#' bamFile <- list.files(bamfolder, pattern = '*.dedup.bam$')
#' bamdir <- file.path(bamfolder, bamFile)
#' sampname_raw = sapply(strsplit(bamFile, '.', fixed = TRUE), '[', 1)
#' bambedObj <- getbambed_scope(bamdir = bamdir,
#'                             bedFile = bedFile,
#'                             sampname = sampname_raw)
#' QCmetric_raw = getsampQC(bambedObj)
#'
#' @author Rujin Wang \email{rujin@email.unc.edu}
#' @import Rsamtools
#' @export
getsampQC = function(bambedObj) {
    ref <- bambedObj$ref
    bamdir <- bambedObj$bamdir
    sampname <- bambedObj$sampname
    QCmetric = matrix(ncol = 8, nrow = length(sampname))
    rownames(QCmetric) = sampname
    colnames(QCmetric) = c("readlength", "total", "mapped",
        "mapped_prop", "non_dup", "non_dup_prop", "mapq20",
        "mapq20_prop")
    for (i in seq_len(length(sampname))) {
        cat("Getting sample QC metric for sample", i, "\n")
        what <- c("rname", "pos", "strand", "mapq", "qwidth", "flag")
        param <- ScanBamParam(what = what)
        aln <- scanBam(bamdir[i], param = param)
        aln <- aln[[1]]
        temp0 = round(mean(aln$qwidth, na.rm = TRUE))
        temp1 = length(aln$mapq)
        temp2 = sum(!is.na(aln$rname))
        temp3 = sum(!is.na(aln$rname) & aln$flag < 1024)
        temp4 = sum(aln$flag < 1024 & aln$mapq >= 20, na.rm = TRUE)
        QCmetric[i, ] = c(temp0, temp1, temp2, round(temp2/temp1, 3),
            temp3, round(temp3/temp1, 3), temp4, round(temp4/temp1,
            3))
    }
    return(QCmetric)
}
