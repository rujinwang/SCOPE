#' @title Get bam file directories, sample names, and whole genomic bins
#'
#' @description Get bam file directories, sample names, and whole genomic
#' bins from .bed file
#'
#' @param bamdir vector of the directory of a bam file. Should be in the same
#'  order as sample names in \code{sampname}.
#' @param sampname vector of sample names. Should be in the same order as bam
#'  directories in \code{bamdir}.
#' @param hgref human reference genome. This should be either 'hg19' or 'hg38'.
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
#' bamfolder <- system.file('extdata', package = 'WGSmapp')
#' bamFile <- list.files(bamfolder, pattern = '*.dedup.bam$')
#' bamdir <- file.path(bamfolder, bamFile)
#' sampname_raw = sapply(strsplit(bamFile, '.', fixed = TRUE), '[', 1)
#' bambedObj <- getbambed_scope(bamdir = bamdir, sampname = sampname_raw)
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
getbambed_scope = function(bamdir, sampname, hgref = "hg19", resolution = 500) {
    if(!hgref %in% c("hg19", "hg38")){
        stop("Reference genome should be either hg19 or hg38!")
    }
    if(resolution <= 0){
        stop("Invalid fixed bin length!")
    }
    if(hgref == "hg19"){
        coordinates = system.file('extdata', 'scWGA500kbsort.bed', 
                                package = 'SCOPE')
    }else{
        coordinates = system.file('extdata', 'scWGA500kb.hg38.bed', 
                                package = 'SCOPE')
    }
    if(resolution == 500){
        bedFile <- read.table(coordinates, sep = "\t", stringsAsFactors = FALSE,
                            header = FALSE)
    } else{
        bed = read.table(coordinates, sep = "\t", stringsAsFactors = FALSE, 
                            header = FALSE)
        names(bed) = c("Chr","Start","End")
        bedFile = NULL
        window_width = resolution * 1000
        for(chrnum in 1:22) {
            bed2 = c()
            chrnum = paste0("chr", chrnum)
            tmp_bed = bed[which(bed$Chr == chrnum),]
            chr_start = min(tmp_bed$Start)
            chr_end = max(tmp_bed$End)
            message("Creating bed file with resolution = ", resolution, "kb for ", 
                    chrnum)
            while(chr_start <= chr_end){
                tmp_end = min(c(chr_start + window_width - 1,chr_end))
                bed2 = rbind(bed2, data.frame(Chr = chrnum, Start = chr_start, 
                                            End = tmp_end, stringsAsFactors = FALSE))
                chr_start = tmp_end + 1
            }
            bedFile = rbind(bedFile, bed2)
        }
    }
    ref <- GRanges(seqnames = bedFile[, 1], ranges =
        IRanges(start = bedFile[, 2], end = bedFile[, 3]))
    if (!any(grepl("chr", seqlevels(ref)))) {
        seqlevels(ref) = paste(c(seq_len(22), "X", "Y"), sep = "")
        ref <- sort(ref)
    } else {
        seqlevels(ref) = paste("chr", c(seq_len(22), "X", "Y"), sep = "")
        ref <- sort(ref)
    }
    list(bamdir = bamdir, sampname = sampname, ref = ref)
}

