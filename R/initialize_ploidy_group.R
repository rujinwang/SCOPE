#' @title Group-wise ploidy pre-initialization
#'
#' @description Pre-estimate ploidies across cells with shared clonal 
#' memberships
#'
#' @usage
#' initialize_ploidy_group(Y, Yhat, ref, groups, 
#'                         maxPloidy = 6, minPloidy = 1.5,
#'                         minBinWidth = 5, SoS.plot = FALSE)
#' @param Y raw read depth matrix after quality control procedure
#' @param Yhat normalized read depth matrix
#' @param ref GRanges object after quality control procedure
#' @param groups clonal membership labels for each cell
#' @param maxPloidy maximum ploidy candidate. Defalut is \code{6}
#' @param minPloidy minimum ploidy candidate. Defalut is \code{1.5}
#' @param minBinWidth the minimum number of bins for a changed segment.
#'  Defalut is \code{5}
#' @param SoS.plot logical, whether to generate ploidy pre-estimation 
#'  plots. Default is \code{FALSE}. 
#'
#' @return
#'     \item{ploidy.SoS}{Vector of group-wise pre-estimated ploidies 
#'             for each cell}
#'
#' @examples
#' Gini <- get_gini(Y_sim)
#'
#' # first-pass CODEX2 run with no latent factors
#' normObj.sim <- normalize_codex2_ns_noK(Y_qc = Y_sim,
#'                                         gc_qc = ref_sim$gc,
#'                                         norm_index = which(Gini<=0.12))
#' Yhat.noK.sim <- normObj.sim$Yhat
#' beta.hat.noK.sim <- normObj.sim$beta.hat
#' fGC.hat.noK.sim <- normObj.sim$fGC.hat
#' N.sim <- normObj.sim$N
#'
#' # Group-wise ploidy initialization
#' clones <- c("normal", "tumor1", "normal", "tumor1", "tumor1")
#' ploidy.sim.group <- initialize_ploidy_group(Y = Y_sim, Yhat = Yhat.noK.sim, 
#'                                 ref = ref_sim, groups = clones)
#' ploidy.sim.group
#'
#' @author Rujin Wang \email{rujin@email.unc.edu}
#' @importFrom DNAcopy CNA smooth.CNA segment
#' @importFrom GenomeInfoDb seqnames
#' @importFrom BSgenome start
#' @importFrom IRanges end
#' @export
initialize_ploidy_group <- function(Y, Yhat, ref, groups, maxPloidy = 6,
                                    minPloidy = 1.5, minBinWidth = 5, 
                                    SoS.plot = FALSE) {
    ploidy.SoS <- rep(NA, ncol(Y))
    breaks <- matrix(0, nrow(Y), ncol(Y))
    RCNP <- matrix(0, nrow(Y), ncol(Y))
    final <- matrix(0, nrow(Y), ncol(Y))
    X <- seq(minPloidy, maxPloidy, by = 0.05)
    n_ploidy <- length(X)
    SoS <- matrix(0, n_ploidy, ncol(Y))
    normal <- (Y + 1)/(Yhat + 1)
    for (k in seq_len(ncol(Y))) {
        if (k%%5 == 1) {
            cat(k, "\t")
        }
    
        lr <- log(normal[, k])
        loc <- data.frame(seq = as.character(seqnames(ref)), 
                        start = start(ref), end = end(ref))
        CNA.object <- CNA(genomdat = lr, chrom = loc[, 1],
                    maploc = as.numeric(loc[, 2]), data.type = "logratio")
        CNA.smoothed <- smooth.CNA(CNA.object)
        segs <- segment(CNA.smoothed, verbose = 0, min.width = minBinWidth)
        frag <- segs$output[, 2:3]
        len <- dim(frag)[1]
        bps <- array(0, len)
        for (j in seq_len(len)) {
            bps[j] <- which((loc[, 1] == frag[j, 1]) &
                            (as.numeric(loc[, 2]) == frag[j, 2]))
        }
        bps <- sort(bps)
        bps[(len = len + 1)] <- nrow(Y)
        breaks[bps, k] <- 1
        RCNP[, k][seq_len(bps[2])] <- median(normal[,
                            k][seq_len(bps[2])])
        for (i in 2:(len - 1)) {
            RCNP[, k][bps[i]:(bps[i + 1] - 1)] <- median(normal[,
                                    k][bps[i]:(bps[i + 1] - 1)])
        }
        RCNP[, k] <- RCNP[, k]/mean(RCNP[, k])
    
        SCNP <- RCNP[, k] %o% X
        FSCP <- round(SCNP)
        Diff2 <- (SCNP - FSCP)^2
        SoS[, k] <- colSums(Diff2, na.rm = FALSE, dims = 1)
    }

    cat("\n", "Initialize ploidy by group", "\n")
    for (G in unique(groups)) {
        cat(G, "\t")
        g <- which(G == groups)
        ploidy.SoS[g] <- X[which.min(apply(SoS[, g, drop = FALSE], 1, sum))]
    
        if(SoS.plot){
            par(mfrow = c(1,2))
            par(mar = c(5,4,4,2))
            hist(apply(Y[, g, drop = FALSE], 1, median), 100, 
                main = 'Read depth distribution', 
                xlab = 'Group-wise Median Coverage per bin')
            plot(X, apply(SoS[, g, drop = FALSE], 1, sum), 
                xlab = "ploidy", ylab = "Group-wise sum of squared errors", 
                main = "First-pass estimation of ploidy", pch = 16)
            abline(v = X[which.min(apply(SoS[, g, drop = FALSE], 1, sum))], 
                lty = 2)
        }
    }
    return(ploidy.SoS)
}