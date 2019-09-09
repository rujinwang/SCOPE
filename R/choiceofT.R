#' @title Visualize EM fitting for each cell.
#'
#' @description A pdf file containing EM fitting results and plots is generated.
#'
#' @usage choiceofT(Y_qc, gc_qc, norm_index, T, ploidyInt, beta0,
#'                 minCountQC = 20, filename)
#' @param Y_qc read depth matrix across all cells after quality control
#' @param gc_qc vector of GC content for each bin after quality control
#' @param norm_index indices of normal/diploid cells
#' @param T a vector of integers indicating number of CNV groups.
#'  Use BIC to select optimal number of CNV groups.
#'  If \code{T = 1}, assume all reads are from normal regions
#'  so that EM algorithm is not implemented. Otherwise,
#'  we assume there is always a CNV group of heterozygous deletion
#'  and a group of null region. The rest groups are representative
#'  of different duplication states.
#' @param ploidyInt a vector of initialized ploidy return from
#'  \code{PreEst_Ploidy}
#' @param beta0 a vector of initialized bin-specific biases returned
#'  from CODEX2 without latent factors
#' @param minCountQC the minimum read coverage required for EM fitting.
#'  Defalut is \code{20}
#' @param filename the name of output pdf file
#'
#' @return pdf file with EM fitting results and two plots:
#' log likelihood, and BIC versus the number of CNV groups.
#'
#' @examples
#' Gini = getGini(Y_sim)
#' # first-pass CODEX2 run with no latent factors
#' normObj.sim <- normalize_codex2_ns_noK(Y_qc =Y_sim,
#'                                         gc_qc = ref_sim$gc,
#'                                         norm_index = which(Gini<=0.12))
#' Yhat.noK.sim=normObj.sim$Yhat
#' beta.hat.noK.sim=normObj.sim$beta.hat
#' fGC.hat.noK.sim=normObj.sim$fGC.hat
#' N.sim = normObj.sim$N
#'
#' # Ploidy initialization
#' ploidy.sim =  PreEst_ploidy(Y = Y_sim, Yhat = Yhat.noK.sim, ref = ref_sim)
#' ploidy.sim
#'
#' choiceofT(Y_qc = Y_sim, gc_qc = ref_sim$gc,
#'         norm_index = which(Gini<=0.12), T = 1:7,
#'         ploidyInt = ploidy.sim,
#'         beta0 = beta.hat.noK.sim,
#'         filename = 'choiceofTdemo.pdf')
#'
#' @author Rujin Wang \email{rujin@email.unc.edu}
#' @import grDevices stats
#' @export
choiceofT = function(Y_qc, gc_qc, norm_index, T, ploidyInt, beta0,
    minCountQC = 20, filename) {
    Y.nonzero <- Y_qc[apply(Y_qc, 1, function(x) {
        !any(x == 0)
    }), ]
    pseudo.sample <- apply(Y.nonzero, 1, function(x) {
        exp(1/length(x) * sum(log(x)))
    })
    Ntotal <- apply(apply(Y.nonzero, 2, function(x) {
        x/pseudo.sample
    }), 2, median)
    N <- round(Ntotal/median(Ntotal) * median(colSums(Y_qc)))
    Nmat <- matrix(nrow = nrow(Y_qc), ncol = ncol(Y_qc), data = N, byrow = TRUE)

    # Get initialization
    gcfit.temp = Y_qc/Nmat/beta0
    alpha0 = matrix(nrow = nrow(Y_qc), ncol = ncol(Y_qc))
    for (j in seq_len(ncol(alpha0))) {
        loe.fit = loess(gcfit.temp[, j] ~ gc_qc)
        gcfit.null = loe.fit$fitted/(ploidyInt[j]/2)
        alpha0[, j] = gcfit.temp[, j]/gcfit.null * 2
    }

    offset = Nmat * matrix(nrow = nrow(Y_qc), ncol = ncol(Y_qc),
        data = beta0, byrow = FALSE)

    pdf(file = filename, width = 8, height = 10)
    for (j in seq_len(ncol(Y_qc))) {
        cat(j, "\t")
        if (j %in% norm_index) {
            fGCj = getfGCj(gcfit.tempj = gcfit.temp[, j], gctemp = gc_qc,
                Yj = Y_qc[, j], offsetj = offset[, j],
                T = 1, draw.plot = TRUE, alphaj = alpha0[, j], minCountQC)
        } else {
            fGCj = getfGCj(gcfit.tempj = gcfit.temp[, j], gctemp = gc_qc,
                Yj = Y_qc[, j], offsetj = offset[, j],
                T = T, draw.plot = TRUE, alphaj = alpha0[, j], minCountQC)
        }
    }
    dev.off()
}
