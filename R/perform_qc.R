#' @title Quality control for cells and bins
#'
#' @description Perform QC step on single cells and bins.
#'
#' @usage
#' perform_qc(Y_raw, sampname_raw, ref_raw, QCmetric_raw,
#'         cov_thresh = 0, minCountQC = 20, 
#'         mapq20_thresh = 0.3, mapp_thresh = 0.9,
#'         gc_thresh = c(20, 80), nMAD = 3)
#'
#' @param Y_raw raw read count matrix returned
#'  from \code{\link{get_coverage_scDNA}}
#' @param sampname_raw sample names for quality control returned
#'  from \code{\link{get_bam_bed}}
#' @param ref_raw raw GRanges object with corresponding GC content
#'  and mappability for quality control returned from
#'  \code{\link{get_bam_bed}}
#' @param QCmetric_raw a QC metric for single cells returned from
#'  \code{\link{get_samp_QC}}
#' @param cov_thresh scalar variable specifying the lower bound of read count
#'  summation of each cell. Default is \code{0}
#' @param minCountQC the minimum read coverage required for
#'  normalization and EM fitting. Defalut is \code{20}
#' @param mapq20_thresh scalar variable specifying the lower threshold
#'  of proportion of reads with mapping quality greater than 20.
#'  Default is \code{0.3}
#' @param mapp_thresh scalar variable specifying mappability of
#'  each genomic bin. Default is \code{0.9}
#' @param gc_thresh vector specifying the lower and upper bound of
#'  GC content threshold for quality control. Default is \code{20-80}
#' @param nMAD scalar variable specifying the number of MAD from the median
#'  of total read counts adjusted by library size for each cell.
#'  Default is \code{3}
#'
#' @return A list with components
#'     \item{Y}{read depth matrix after quality control}
#'     \item{sampname}{sample names after quality control}
#'     \item{ref}{A GRanges object specifying whole genomic
#'         bin positions after quality control}
#'     \item{QCmetric}{A data frame of QC metric for single cells
#'         after quality control}
#'
#' @examples
#' Y_raw <- coverageObj.scopeDemo$Y
#' sampname_raw <- rownames(QCmetric.scopeDemo)
#' ref_raw <- ref.scopeDemo
#' QCmetric_raw <- QCmetric.scopeDemo
#' qcObj <- perform_qc(Y_raw = Y_raw, sampname_raw = sampname_raw,
#'                 ref_raw = ref_raw, QCmetric_raw = QCmetric_raw)
#'
#' @author Rujin Wang \email{rujin@email.unc.edu}
#' @import stats
#' @export
perform_qc <- function(Y_raw, sampname_raw, ref_raw, QCmetric_raw,
    cov_thresh = 0, minCountQC = 20, 
    mapq20_thresh = 0.3, mapp_thresh = 0.9,
    gc_thresh = c(20, 80), nMAD = 3) {
    if (length(ref_raw) != nrow(Y_raw)) {
        stop("Invalid inputs: length of ref and # of rows
            in read count matrix must be the same")
    }
    if (length(sampname_raw) != ncol(Y_raw)) {
        stop("Invalid inputs: length of sample names and # of cols in
            read count matrix must be the same")
    }
    if (nrow(QCmetric_raw) != ncol(Y_raw)) {
        stop("Invalid inputs: # of rows in QC metric and # of cols in
            read count matrix must be the same")
    }
    mapp <- ref_raw$mapp
    gc <- ref_raw$gc
    sampfilter1 <- (apply(Y_raw, 2, sum) <= cov_thresh)
    message("Removed ", sum(sampfilter1),
        " samples due to failed library preparation.")
    sampfilter2 <- (apply(Y_raw, 2, mean) <= minCountQC)
    message("Removed ", sum(sampfilter2),
        " samples due to failure to meet min coverage requirement.")    
    sampfilter3 <- (QCmetric_raw[, "mapq20_prop"] < mapq20_thresh)
    message("Removed ", sum(sampfilter3),
        " samples due to low proportion of mapped reads.")
    if (sum(sampfilter1 | sampfilter2 | sampfilter3) != 0) {
        Y <- Y_raw[, !(sampfilter1 | sampfilter2 | sampfilter3)]
        sampname <- sampname_raw[!(sampfilter1 | sampfilter2 
                                | sampfilter3)]
        QCmetric <- QCmetric_raw[!(sampfilter1 | sampfilter2 
                                | sampfilter3), ]
    } else {
        Y <- Y_raw
        sampname <- sampname_raw
        QCmetric <- QCmetric_raw
    }
    binfilter1 <- (gc < gc_thresh[1] | gc > gc_thresh[2])
    message("Excluded ", sum(binfilter1),
        " bins due to extreme GC content.")
    binfilter2 <- (mapp < mapp_thresh)
    message("Excluded ", sum(binfilter2),
        " bins due to low mappability.")
    if (sum(binfilter1 | binfilter2) != 0) {
        ref <- ref_raw[!(binfilter1 | binfilter2)]
        Y <- Y[!(binfilter1 | binfilter2), ]
    } else {
        ref <- ref_raw
        Y <- Y
    }
    Y.nonzero <- Y[apply(Y, 1, function(x) {
        !any(x == 0)
        }), , drop = FALSE]
    if(dim(Y.nonzero)[1] <= 10){
        message("Adopt arithmetic mean instead of geometric mean")
        pseudo.sample <- apply(Y, 1, mean)
        N <- apply(apply(Y, 2, function(x) {
            x/pseudo.sample
        }), 2, median, na.rm = TRUE)
    } else{
        pseudo.sample <- apply(Y.nonzero, 1, function(x) {
            exp(sum(log(x))/length(x))
        })
        N <- apply(apply(Y.nonzero, 2, function(x) {
            x/pseudo.sample
        }), 2, median)
    }
    sampfilter3 <- (N == 0)
    message("Removed ", sum(sampfilter3),
            " samples due to excessive zero read counts in 
            library size calculation.")
    if (sum(sampfilter3) != 0) {
        Y <- Y[, !(sampfilter3)]
        sampname <- sampname[!(sampfilter3)]
        QCmetric <- QCmetric[!(sampfilter3), ]
        N <- N[!(sampfilter3)]
    }
    Nmat <- matrix(nrow = nrow(Y), ncol = ncol(Y), data = N,
        byrow = TRUE)
    bin.sum <- apply(Y/Nmat, 1, sum)
    binfilter3 <- (bin.sum >= (median(bin.sum) -
        nMAD * mad(bin.sum))) & (bin.sum <= (median(bin.sum) +
        nMAD * mad(bin.sum)))
    Y <- Y[binfilter3, ]
    ref <- ref[binfilter3]
    QCmetric <- as.data.frame(QCmetric)
    message("There are ", ncol(Y), " samples and ",
        nrow(Y), " bins after QC step. ")
    list(Y = Y, sampname = sampname, ref = ref, QCmetric = QCmetric)
}
