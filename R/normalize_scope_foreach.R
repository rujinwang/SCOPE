#' @title Normalization of read depth with latent factors using
#' Expectation-Maximization algorithm under the case-control
#' setting in parallel
#'
#' @description Fit a Poisson generalized linear model to normalize
#' the raw read depth data from single-cell DNA sequencing,
#' with latent factors under the case-control setting. Model GC
#' content bias using an expectation-maximization algorithm,
#' which accounts for the different copy number states.
#'
#' @usage
#' normalize_scope_foreach(Y_qc, gc_qc, K, norm_index, T,
#'     ploidyInt, beta0, minCountQC = 20, nCores = NULL)
#' @param Y_qc read depth matrix after quality control
#' @param gc_qc vector of GC content for each bin after quality control
#' @param K Number of latent Poisson factors
#' @param norm_index indices of normal/diploid cells
#' @param T a vector of integers indicating number of CNV groups.
#'  Use BIC to select optimal number of CNV groups. If \code{T = 1},
#'  assume all reads are from normal regions so that EM algorithm is
#'  not implemented. Otherwise, we assume there is always a CNV group
#'  of heterozygous deletion and a group of null region. The rest groups
#'  are representative of different duplication states.
#' @param ploidyInt a vector of initialized ploidy return
#'  from \code{PreEst_Ploidy}
#' @param beta0 a vector of initialized bin-specific biases returned
#'  from CODEX2 without latent factors
#' @param minCountQC the minimum read coverage required for normalization
#'  and EM fitting. Defalut is \code{20}
#' @param nCores number of cores to use. If \code{NULL}, number of cores
#'  is detected. Default is \code{NULL}.
#'
#' @return A list with components
#'     \item{Yhat}{A list of normalized read depth matrix with EM}
#'     \item{alpha.hat}{A list of absolute copy number matrix}
#'     \item{fGC.hat}{A list of EM estimated GC content bias matrix}
#'     \item{beta.hat}{A list of EM estimated bin-specific bias vector}
#'     \item{g.hat}{A list of estimated Poisson latent factor}
#'     \item{h.hat}{A list of estimated Poisson latent factor}
#'     \item{AIC}{AIC for model selection}
#'     \item{BIC}{BIC for model selection}
#'     \item{RSS}{RSS for model selection}
#'     \item{K}{Number of latent Poisson factors}
#'
#' @examples
#' Gini = getGini(Y_sim)
#'
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
#' ploidy.sim = PreEst_ploidy(Y = Y_sim,
#'                             Yhat = Yhat.noK.sim,
#'                             ref = ref_sim)
#' ploidy.sim
#'
#' # Specify nCores = 2 only for checking examples
#' normObj.scope.sim = normalize_scope_foreach(Y_qc = Y_sim, gc_qc = ref_sim$gc,
#'                         K = 1, ploidyInt = ploidy.sim,
#'                         norm_index = which(Gini<=0.12), T = 1:7,
#'                         beta0 = beta.hat.noK.sim, nCores = 2)
#' Yhat.sim = normObj.scope.sim$Yhat[[which.max(normObj.scope.sim$BIC)]]
#' fGC.hat.sim = normObj.scope.sim$fGC.hat[[which.max(normObj.scope.sim$BIC)]]
#'
#' @author Rujin Wang \email{rujin@email.unc.edu}
#' @import stats foreach parallel doParallel
#' @export
normalize_scope_foreach = function(Y_qc, gc_qc, K, norm_index, T,
    ploidyInt, beta0, minCountQC = 20, nCores = NULL) {
    if (max(K) > length(norm_index))
        stop("Number of latent Poisson factors K cannot
            exceed the number of normal samples!")
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
    Nmat <- matrix(nrow = nrow(Y_qc), ncol = ncol(Y_qc),
        data = N, byrow = TRUE)

    # Get initialization
    gcfit.temp = Y_qc/Nmat/beta0
    alpha0 = matrix(nrow = nrow(Y_qc), ncol = ncol(Y_qc))
    for (j in seq_len(ncol(alpha0))) {
        loe.fit = loess(gcfit.temp[, j] ~ gc_qc)
        gcfit.null = loe.fit$fitted/(ploidyInt[j]/2)
        alpha0[, j] = gcfit.temp[, j]/gcfit.null * 2
    }

    Yhat = vector("list", length(K))
    fGC.hat <- vector("list", length(K))
    alpha.hat <- vector("list", length(K))
    beta.hat <- vector("list", length(K))
    g.hat <- vector("list", length(K))
    h.hat <- vector("list", length(K))
    AIC <- rep(NA, length = length(K))
    BIC <- rep(NA, length = length(K))
    RSS <- rep(NA, length = length(K))

    # Initialization
    message("Initialization ...")
    gcfit.temp = Y_qc/Nmat/beta0
    offset = Nmat * matrix(nrow = nrow(Y_qc), ncol = ncol(Y_qc),
        data = beta0, byrow = FALSE)
    fhat.temp = getfGC_foreach(gcfit.temp = gcfit.temp, gctemp = gc_qc,
        Y = Y_qc, norm_index = norm_index, offset = offset,
        T = T, alpha = alpha0, minCountQC = minCountQC, nCores = nCores)
    fhat0 = fhat.temp$fGC.hat
    alpha0 = fhat.temp$alpha

    for (ki in seq_len(length(K))) {
        k = K[ki]
        message("Computing normalization with k = ", k,
            " latent factors ...", sep = "")
        message("k = ", k)
        maxiter = 10
        maxhiter = 50
        BHTHRESH = 1e-04
        HHTHRESH = 1e-05
        iter = 1
        fhat = matrix(nrow = nrow(Y_qc), ncol = ncol(Y_qc), data = 0)
        betahat = beta0
        betahatmat = matrix(nrow = nrow(Y_qc), ncol = ncol(Y_qc),
            data = betahat, byrow = FALSE)
        ghat = matrix(0, nrow = nrow(Y_qc), ncol = k)
        hhat = matrix(0, nrow = ncol(Y_qc), ncol = k)
        bhdiff = rep(Inf, maxiter)
        fhdiff = rep(Inf, maxiter)

        betahatlist = vector("list", maxiter)
        fhatlist = vector("list", maxiter)
        ghatlist = vector("list", maxiter)
        hhatlist = vector("list", maxiter)
        alphahatlist = vector("list", maxiter)

        while (iter <= maxiter) {
            if (iter == 1) {
                fhatnew = fhat0
                alpha = alpha0
            }
            if (iter > 1) {
                gcfit.temp = Y_qc/Nmat/betahat/exp(ghat %*% t(hhat))
                offset = Nmat * matrix(nrow = nrow(Y_qc), ncol = ncol(Y_qc),
                    data = betahat, byrow = FALSE) * exp(ghat %*% t(hhat))
                fhat.temp = getfGC_foreach(gcfit.temp = gcfit.temp,
                    gctemp = gc_qc,
                    Y = Y_qc, norm_index = norm_index,
                    offset = offset, T = T, alpha = alpha0,
                    minCountQC = minCountQC,
                    nCores = nCores)
                fhatnew = fhat.temp$fGC.hat
                alpha = fhat.temp$alpha
            }
            fhatnew[fhatnew < quantile(fhatnew, 0.005)] = quantile(
                fhatnew, 0.005)
            betahatnew = apply((Y_qc/(fhatnew * Nmat *
                exp(ghat %*% t(hhat))))[, norm_index], 1, median)
            betahatnew[betahatnew <= 0] = min(betahatnew[betahatnew > 0])
            bhdiff[iter] = sum((betahatnew - betahat)^2)/length(betahat)
            fhdiff[iter] = sum((fhatnew - fhat)^2)/length(fhat)
            if (fhdiff[iter] > min(fhdiff))
                break
            message("Iteration ", iter, "\t", "beta diff =",
                signif(bhdiff[iter], 3), "\t", "f(GC) diff =",
                signif(fhdiff[iter],
                3))
            fhat = fhatnew
            betahat = betahatnew
            betahatmat = matrix(nrow = nrow(Y_qc), ncol = ncol(Y_qc),
                data = betahat, byrow = FALSE)
            L = log(Nmat * fhat * betahatmat * alpha/2)
            logmat = log(pmax(Y_qc, 1)) - L
            logmat = logmat - matrix(nrow = nrow(Y_qc), ncol = ncol(Y_qc),
                data = apply(logmat, 1, mean), byrow = FALSE)
            hhat = svd(logmat, nu = k, nv = k)$v
            hhatnew = hhat
            hiter = 1
            hhdiff = rep(Inf, maxhiter)
            while (hiter <= maxhiter) {
                for (s in seq_len(nrow(Y_qc))) {
                    temp = try(glm(formula = Y_qc[s, norm_index] ~
                        hhat[norm_index, ] -
                        1, offset = L[s, norm_index],
                        family = poisson)$coefficients, silent = TRUE)
                    if (is.character(temp)) {
                        temp = lm(log(pmax(Y_qc[s, norm_index],
                            1)) ~ hhat[norm_index, ] - 1,
                            offset = log(L[s, norm_index]))$coefficients
                    }
                    ghat[s, ] = temp
                }
                # avoid overflow or underflow of the g latent factors
                ghat[is.na(ghat)] = 0
                if (max(ghat) >= 30) {
                    ghat = apply(ghat, 2, function(z) {
                        z[z > quantile(z, 0.995)] = min(quantile(z, 0.995), 30)
                        z
                    })
                }
                if (min(ghat) <= -30) {
                    ghat = apply(ghat, 2, function(z) {
                        z[z < quantile(z, 0.005)] = max(quantile(z, 0.005), -30)
                        z
                    })
                }
                for (t in seq_len(ncol(Y_qc))) {
                    hhatnew[t, ] = glm(formula = Y_qc[, t] ~ ghat - 1,
                        offset = L[, t], family = poisson)$coefficients
                }
                gh = ghat %*% t(hhatnew)
                gh <- scale(gh, center = TRUE, scale = FALSE)
                hhatnew = svd(gh, nu = k, nv = k)$v
                hhdiff[hiter] = sum((hhatnew - hhat)^2)/length(hhat)
                message("\t\t\t", "hhat diff =", signif(hhdiff[hiter], 3))
                hhat = hhatnew
                if (hhdiff[hiter] < HHTHRESH)
                    break
                if (hiter > 10 & (rank(hhdiff))[hiter] <= 3)
                    break
                hiter = hiter + 1
            }
            alphahatlist[[iter]] = alpha
            fhatlist[[iter]] = fhat
            betahatlist[[iter]] = betahat
            ghatlist[[iter]] = ghat
            hhatlist[[iter]] = hhat
            if (bhdiff[iter] < BHTHRESH)
                break
            if (iter > 5 & bhdiff[iter] > 1)
                break
            iter = iter + 1
        }
        optIter = which.min(fhdiff)
        message(paste("Stop at Iteration ", optIter, ".", sep = ""))
        alpha = alphahatlist[[optIter]]
        fhat = fhatlist[[optIter]]
        betahat = betahatlist[[optIter]]
        ghat = ghatlist[[optIter]]
        hhat = hhatlist[[optIter]]
        betahatmat = matrix(nrow = nrow(Y_qc), ncol = ncol(Y_qc),
            data = betahat, byrow = FALSE)
        Yhat[[ki]] = pmax(round(fhat * Nmat * betahatmat *
            exp(ghat %*% t(hhat)), 0), 1)
        alpha.hat[[ki]] <- alpha
        fGC.hat[[ki]] <- signif(fhat, 3)
        beta.hat[[ki]] <- signif(betahat, 3)
        h.hat[[ki]] <- signif(hhat, 3)
        g.hat[[ki]] <- signif(ghat, 3)
        Yhat.temp = Yhat[[ki]] * alpha/2
        AIC[ki] = 2 * sum(Y_qc * log(pmax(Yhat.temp, 1)) - Yhat.temp) -
            2 * (length(ghat) + length(hhat))
        BIC[ki] = 2 * sum(Y_qc * log(pmax(Yhat.temp, 1)) - Yhat.temp) -
            (length(ghat) + length(hhat)) * log(length(Y_qc))
        RSS[ki] = sum((Y_qc - Yhat.temp)^2/length(Y_qc))
        message("AIC", k, " = ", round(AIC[ki], 3))
        message("BIC", k, " = ", round(BIC[ki], 3))
        message("RSS", k, " = ", round(RSS[ki], 3), "\n")
    }
    list(Yhat = Yhat, alpha.hat = alpha.hat, fGC.hat = fGC.hat,
        beta.hat = beta.hat, g.hat = g.hat, h.hat = h.hat,
        AIC = AIC, BIC = BIC, RSS = RSS, K = K)
}



getfGCj_foreach = function(gcfit.tempj, gctemp, Yj, offsetj, T,
    draw.plot = NULL, alphaj, minCountQC) {
    alphaj = pmax(1, round(alphaj))
    if (is.null(draw.plot)) {
        draw.plot = FALSE
    }
    fGCi = fitGC(gctemp, gcfit.tempj)
    resid = abs(gcfit.tempj - fGCi)

    bin.filter = which(resid > (median(resid) +
        5 * mad(resid)) | Yj < minCountQC)
    if (length(bin.filter) == 0) {
        bin.filter = which.max(gcfit.tempj)
    }
    Ti = T
    if (Ti == 1) {
        Z = matrix(nrow = length(gcfit.tempj), ncol = Ti, data = 1/Ti)
        vec_pi = 1
        loe.fit.temp = loess(gcfit.tempj[-bin.filter] ~ gctemp[-bin.filter])
        fGCi = predict(loe.fit.temp, newdata = gctemp, se = TRUE)$fit
        temp = min(fGCi[!is.na(fGCi) & fGCi > 0])
        fGCi[fGCi <= 0 | is.na(fGCi)] = temp
    }
    if (Ti >= 2) {
        Z = matrix(nrow = length(gcfit.tempj), ncol = Ti, data = 0)
        mintemp = pmin(Ti, alphaj)
        Z[cbind(seq_len(nrow(Z)), mintemp)] = 1
        vec_pi = colSums(Z)/nrow(Z)

        loe.fit.temp = loess((gcfit.tempj/(Z %*%
            as.matrix(seq_len(Ti)/2)))[-bin.filter] ~ gctemp[-bin.filter])
        fGCi = predict(loe.fit.temp, newdata = gctemp, se = TRUE)$fit
        temp = min(fGCi[!is.na(fGCi) & fGCi > 0])
        fGCi[fGCi <= 0 | is.na(fGCi)] = temp

        diff.GC = Inf
        diff.Z = Inf
        iter = 1
        while (iter <= 3 | diff.GC > 5e-06 | diff.Z > 0.005) {
            Mtemp = Mstep(Z, gcfit.tempj, gctemp)
            vec_pi.new = Mtemp$vec_pi
            fGCi.new = Mtemp$fGCi
            Z.new = Estep(fGCi.new, vec_pi.new, Yj, offsetj)
            diff.GC = sum((fGCi - fGCi.new)^2)/length(fGCi)
            diff.Z = sum((Z.new - Z)^2)/length(Z)
            vec_pi = vec_pi.new
            Z = Z.new
            fGCi = fGCi.new
            iter = iter + 1
            if (iter >= 50)
                break
        }
    }

    if (Ti == 1) {
        loe.fit.plot = loess(gcfit.tempj ~ gctemp)
        fGCi.plot = loe.fit.plot$fitted
        temp = min(fGCi.plot[!is.na(fGCi.plot) & fGCi.plot > 0])
        fGCi.plot[fGCi.plot <= 0 | is.na(fGCi.plot)] = temp
        df = predict(loe.fit.plot, newdata = gctemp, se = TRUE)$df

        loglik = sum(dpois(Yj[-bin.filter],
            lambda = (offsetj * fGCi)[-bin.filter], log = TRUE))
        BIC = 2 * loglik - (length(gcfit.tempj) - df) *
            log(length(gcfit.tempj))
    } else {
        loe.fit.plot = loess((gcfit.tempj/(Z %*%
            as.matrix(seq_len(Ti)/2))) ~ gctemp)
        fGCi.plot = loe.fit.plot$fitted
        temp = min(fGCi.plot[!is.na(fGCi.plot) & fGCi.plot > 0])
        fGCi.plot[fGCi.plot <= 0 | is.na(fGCi.plot)] = temp
        df = predict(loe.fit.plot, newdata = gctemp, se = TRUE)$df

        loglik = sum(dpois(Yj[-bin.filter],
            lambda = (offsetj * fGCi * (Z %*%
            as.matrix(seq_len(Ti)/2)))[-bin.filter],
            log = TRUE))
        BIC = 2 * loglik - (length(gcfit.tempj) - df + Ti - 1) *
            log(length(gcfit.tempj))
    }

    if (draw.plot) {
        smoothScatter(gctemp[-bin.filter], gcfit.tempj[-bin.filter],
            xlab = "GC content", ylab = "Y/beta/N/exp(gxh)",
            nrpoints = 0, main = paste("T =", Ti))
        if (Ti == 1) {
            points(gctemp[order(gctemp)], fGCi.plot[order(gctemp)],
                lty = 2, col = 2, type = "l", lwd = 2)
            points(gctemp, gcfit.tempj, cex = 0.4, col = 2, pch = 16)
        } else {
            for (k in seq_len(Ti)) {
                points(gctemp[order(gctemp)], fGCi.plot[order(gctemp)] *
                    k/2, lty = 2, col = k, type = "l", lwd = 2)
                points(gctemp[which((round(Z))[, k] == 1)],
                    (gcfit.tempj)[which((round(Z))[, k] == 1)], cex = 0.4,
                    col = k, pch = 16)
            }
        }
    }

    fGCi.obj = fGCi
    Z.obj = Z
    vec_pi.obj = vec_pi
    return(list(fGCi.obj = fGCi.obj, Z.obj = Z.obj,
        vec_pi.obj = vec_pi.obj, bin.filter = bin.filter,
        loglik = loglik, BIC = BIC))
}



if (getRversion() >= "2.15.1") {
    utils::globalVariables(c("Ti"))
}
getfGC_foreach = function(gcfit.temp, gctemp, Y, norm_index, offset, T,
    alpha, minCountQC, nCores) {
    fGC.hat = matrix(ncol = ncol(Y), nrow = nrow(Y))
    for (j in seq_len(ncol(Y))) {
        cat(j, "\t")
        if (j %in% norm_index) {
            alpha[, j] = 2
            loe.fit = loess(gcfit.temp[, j] ~ gctemp)
            fGC.hat[, j] = loe.fit$fitted
        } else {
            # begin foreach
            if (is.null(nCores)) {
                nCores <- detectCores() - 1
            }
            registerDoParallel(nCores)

            TList <- foreach(Ti = T, .export = c("fitGC", "Estep", "Mstep",
                "getfGCj_foreach")) %dopar% {
                getfGCj_foreach(gcfit.tempj = gcfit.temp[, j], gctemp = gctemp,
                    Yj = Y[, j], offsetj = offset[,j], T = Ti,
                    draw.plot = FALSE,
                    alphaj = alpha[, j], minCountQC = minCountQC)
            }

            # When you're done, clean up the cluster
            stopImplicitCluster()
            # end foreach

            fGCi.obj = lapply(TList, function(x) x[["fGCi.obj"]])
            Z.obj = lapply(TList, function(x) x[["Z.obj"]])
            vec_pi.obj = lapply(TList, function(x) x[["vec_pi.obj"]])
            loglik = vapply(TList, function(x) x[["loglik"]], numeric(1))
            BIC = vapply(TList, function(x) x[["BIC"]], numeric(1))

            fGCj = list(fGCi.obj = fGCi.obj, Z.obj = Z.obj,
                vec_pi.obj = vec_pi.obj, loglik = loglik, BIC = BIC)
            if (which.max(fGCj$BIC) == 1) {
                alpha[, j] = 2
            } else {
                alpha[, j] = apply(fGCj$Z.obj[[which.max(fGCj$BIC)]],
                    1, which.max)
            }
            fGC.hat[, j] = fGCj$fGCi.obj[[which.max(fGCj$BIC)]]
        }
    }
    return(list(fGC.hat = fGC.hat, alpha = alpha))
}
