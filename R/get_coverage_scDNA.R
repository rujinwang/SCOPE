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
#' @param hgref reference genome. This should be 'hg19', 'hg38' or 'mm10'.
#' Default is human genome \code{hg19}. 
#' @return
#'   \item{Y}{Read depth matrix}
#'
#' @examples
#' library(WGSmapp)
#' library(BSgenome.Hsapiens.UCSC.hg38)
#' bamfolder <- system.file('extdata', package = 'WGSmapp')
#' bamFile <- list.files(bamfolder, pattern = '*.dedup.bam$')
#' bamdir <- file.path(bamfolder, bamFile)
#' sampname_raw <- sapply(strsplit(bamFile, '.', fixed = TRUE), '[', 1)
#' bambedObj <- get_bam_bed(bamdir = bamdir,
#'                             sampname = sampname_raw, 
#'                             hgref = "hg38")
#'
#' # Getting raw read depth
#' coverageObj <- get_coverage_scDNA(bambedObj,
#'                                 mapqthres = 40,
#'                                 seq = 'paired-end', 
#'                                 hgref = "hg38")
#' Y_raw <- coverageObj$Y
#'
#' @author Rujin Wang \email{rujin@email.unc.edu}
#' @import Rsamtools
#' @importFrom GenomicRanges GRanges
#' @importFrom IRanges IRanges RangesList Views countOverlaps
#' @importFrom GenomeInfoDb seqnames
#' @export
get_coverage_scDNA <- function(bambedObj, mapqthres, seq, hgref = "hg19") {
    if(!hgref %in% c("hg19", "hg38", "mm10")){
        stop("Reference genome should be hg19, hg38 or mm10.")
    } 
    ref <- bambedObj$ref
    bamdir <- bambedObj$bamdir
    sampname <- bambedObj$sampname
    
    mask.ref <- get_masked_ref(hgref = hgref)
    
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


get_masked_ref <- function(hgref){
    if(hgref == "hg19"){
        # Get segmental duplication regions
        seg.dup <- read.table(system.file("extdata", 
                                "GRCh37GenomicSuperDup.tab", 
                                package = "WGSmapp"), header = TRUE)
        # Get hg19 gaps
        gaps <- read.table(system.file("extdata", "hg19gaps.txt", 
                                package = "WGSmapp"), 
                                header = TRUE)
    } else if(hgref == "hg38"){
        # Get segmental duplication regions
        seg.dup <- read.table(system.file("extdata", 
                                "GRCh38GenomicSuperDup.tab", 
                                package = "WGSmapp"))
        # Get hg19 gaps
        gaps <- read.table(system.file("extdata", "hg38gaps.txt", 
                                package = "WGSmapp"))
    } else if (hgref == "mm10"){
        black.list <- read.table(system.file("extdata", 
                                "mm10-blacklist.v2.bed", 
                                package = "WGSmapp"), header = FALSE, 
                                sep = '\t')
    }
    if(hgref != "mm10"){
        seg.dup <- seg.dup[!is.na(match(seg.dup[,1], 
                            paste('chr', c(seq_len(22), 'X', 'Y'), 
                            sep = ''))),]
        seg.dup <- GRanges(seqnames = as.character(seg.dup[,1]), 
                            ranges = IRanges(start = seg.dup[,2], 
                            end = seg.dup[,3]))
        gaps <- gaps[!is.na(match(gaps[,2], 
                            paste('chr', c(seq_len(22), 'X', 'Y'), 
                            sep = ''))),]
        gaps <- GRanges(seqnames = as.character(gaps[,2]), 
                            ranges = IRanges(start = gaps[,3], 
                            end = gaps[,4]))
        # Generate mask region
        mask.ref <- sort(c(seg.dup, gaps))
    } else{
        black.list <- black.list[!is.na(match(black.list[,1], 
                paste('chr', c(seq_len(19), 'X', 'Y'), sep = ''))),]
        black.list <- GRanges(seqnames = black.list[,1], 
                ranges = IRanges(start = black.list[,2], 
                end = black.list[,3]))
        mask.ref <- black.list
    }
    return(mask.ref)
}


