#' A read count matrix in the toy dataset
#'
#' @docType data
#'
#' @format A read count matrix with 1544 bins and 39 cells
#'
#' @keywords datasets
"Y_sim"


#' A reference genome in the toy dataset
#'
#' @docType data
#'
#' @format A GRanges object with 1544 bins and 1 metadata column of GC content
#'
#' @keywords datasets
"ref_sim"


#' A post cross-sample segmentation integer copy number matrix returned by
#' SCOPE in the demo
#'
#' @docType data
#'
#' @format A post cross-sample segmentation integer copy number matrix of
#' five toy cells returned by SCOPE
#'
#' @keywords datasets
"iCN_sim"


#' Pre-stored normObj.scope data for demonstration purposes
#'
#' @docType data
#'
#' @format Pre-computed by SCOPE using pre-stored data \code{Y_sim}
#'
#' @keywords datasets
"normObj.scopeDemo"


#' Pre-stored coverageObj.scope data for demonstration purposes
#'
#' @docType data
#'
#' @format Pre-computed using whole genome sequencing data of
#' three single cells from 10X Genomics Single-Cell CNV solution
#'
#' @keywords datasets
"coverageObj.scopeDemo"


#' Pre-stored QCmetric data for demonstration purposes
#'
#' @docType data
#'
#' @format Pre-computed using whole genome sequencing data of
#' three single cells from 10X Genomics Single-Cell CNV solution
#'
#' @keywords datasets
"QCmetric.scopeDemo"




#' Pre-stored 500kb-size reference genome for demonstration purposes
#'
#' @docType data
#'
#' @format Pre-computed using whole genome sequencing data
#' with GC content and mappability scores
#'
#' @keywords datasets
"ref.scopeDemo"
