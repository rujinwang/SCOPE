#' @title Group-wise normalization of read depth with latent factors using
#' Expectation-Maximization algorithm and shared clonal memberships
#'
#' @description Fit a Poisson generalized linear model to normalize
#' the raw read depth data from single-cell DNA sequencing, with
#' latent factors and shared clonal memberships. Model GC content
#' bias using an expectation-maximization algorithm, which accounts for
#' clonal specific copy number states.
#'
#' @usage
#' normalize_scope_group(Y_qc, gc_qc, K, norm_index, groups, T, 
#'                         ploidyInt, beta0, minCountQC = 20)
#' @param Y_qc read depth matrix after quality control
#' @param gc_qc vector of GC content for each bin after quality control
#' @param K Number of latent Poisson factors
#' @param norm_index indices of normal/diploid cells using group/clone
#' labels
#' @param groups clonal membership labels for each cell
#' @param T a vector of integers indicating number of CNV groups.
#'  Use BIC to select optimal number of CNV groups. If \code{T = 1},
#'  assume all reads are from normal regions so that EM algorithm is
#'  not implemented. Otherwise, we assume there is always a CNV group
#'  of heterozygous deletion and a group of null region. The rest
#'  groups are representative of different duplication states.
#' @param ploidyInt a vector of group-wise initialized ploidy return
#'  from \code{initialize_ploidy_group}. Users are also allowed to 
#'  provide prior-knowledge ploidies as the input and to manually 
#'  tune a few cells/clones that have poor fitting
#' @param beta0 a vector of initialized bin-specific biases
#'  returned from CODEX2 without latent factors
#' @param minCountQC the minimum read coverage required for
#'  normalization and EM fitting. Defalut is \code{20}
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
#' normObj.scope.sim.group <- normalize_scope_group(Y_qc = Y_sim, 
#'                                     gc_qc = ref_sim$gc,
#'                                     K = 1, ploidyInt = ploidy.sim.group,
#'                                     norm_index = which(clones=="normal"), 
#'                                     groups = clones, 
#'                                     T = 1:5,
#'                                     beta0 = beta.hat.noK.sim)
#' Yhat.sim.group <- normObj.scope.sim.group$Yhat[[which.max(
#'                                     normObj.scope.sim.group$BIC)]]
#' fGC.hat.sim.group <- normObj.scope.sim.group$fGC.hat[[which.max(
#'                                     normObj.scope.sim.group$BIC)]]
#'
#' @author Rujin Wang \email{rujin@email.unc.edu}
#' @import stats
#' @export
normalize_scope_group <- function(Y_qc, gc_qc, K, norm_index, groups, T, 
                                    ploidyInt, beta0, minCountQC = 20) {
    if (max(K) > length(norm_index))
        stop("Number of latent Poisson factors K cannot exceed the number of
        normal samples. ")
    Y.nonzero <- Y_qc[apply(Y_qc, 1, function(x) {
        !any(x == 0)
    }), , drop = FALSE]
    if(dim(Y.nonzero)[1] <= 10){
        message("Adopt arithmetic mean instead of geometric mean")
        pseudo.sample <- apply(Y_qc, 1, mean)
        Ntotal <- apply(apply(Y_qc, 2, function(x) {
            x/pseudo.sample
            }), 2, median, na.rm = TRUE)
    } else{
        pseudo.sample <- apply(Y.nonzero, 1, function(x) {
            exp(sum(log(x))/length(x))
        })
        Ntotal <- apply(apply(Y.nonzero, 2, function(x) {
            x/pseudo.sample
        }), 2, median)
    }
    N <- round(Ntotal/median(Ntotal) * median(colSums(Y_qc)))
    Nmat <- matrix(nrow = nrow(Y_qc), ncol = ncol(Y_qc), data = N, 
                    byrow = TRUE)

    # Get initialization
    gcfit.temp <- Y_qc/Nmat/beta0
    alpha0 <- matrix(nrow = nrow(Y_qc), ncol = ncol(Y_qc))
    for (j in seq_len(ncol(alpha0))) {
        loe.fit <- loess(gcfit.temp[, j] ~ gc_qc)
        gcfit.null <- loe.fit$fitted/(ploidyInt[j]/2)
        alpha0[, j] <- gcfit.temp[, j]/gcfit.null * 2
    }
    
    Yhat <- vector("list", length(K))
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
    gcfit.temp <- Y_qc/Nmat/beta0
    offset <- Nmat * matrix(nrow = nrow(Y_qc), ncol = ncol(Y_qc),
                            data = beta0, byrow = FALSE)

    fhat.temp <- getfGC2(gcfit.temp = gcfit.temp, gctemp = gc_qc,
                        Y = Y_qc, norm_index = norm_index, groups = groups, 
                        offset = offset,
                        T = T, alpha = alpha0, minCountQC = minCountQC)
    fhat0 <- fhat.temp$fGC.hat
    alpha0 <- fhat.temp$alpha

    for (ki in seq_len(length(K))) {
        k <- K[ki]
        message("Computing normalization with k = ", k,
                " latent factors ...", sep = "")
        message("k = ", k)
        maxiter <- 10
        maxhiter <- 50
        BHTHRESH <- 1e-04
        HHTHRESH <- 1e-05
        iter <- 1
        fhat <- matrix(nrow = nrow(Y_qc), ncol = ncol(Y_qc), data = 0)
        betahat <- beta0
        betahatmat <- matrix(nrow = nrow(Y_qc), ncol = ncol(Y_qc),
                            data = betahat, byrow = FALSE)
        ghat <- matrix(0, nrow = nrow(Y_qc), ncol = k)
        hhat <- matrix(0, nrow = ncol(Y_qc), ncol = k)
        bhdiff <- rep(Inf, maxiter)
        fhdiff <- rep(Inf, maxiter)
    
        betahatlist <- vector("list", maxiter)
        fhatlist <- vector("list", maxiter)
        ghatlist <- vector("list", maxiter)
        hhatlist <- vector("list", maxiter)
        alphahatlist <- vector("list", maxiter)
    
        while (iter <= maxiter) {
            if (iter == 1) {
                fhatnew <- fhat0
                alpha <- alpha0
            }
            if (iter > 1) {
                gcfit.temp <- Y_qc/Nmat/betahat/exp(ghat %*% t(hhat))
                offset <- Nmat * matrix(nrow = nrow(Y_qc), ncol = ncol(Y_qc),
                        data = betahat, byrow = FALSE) * exp(ghat %*% t(hhat))
                fhat.temp <- getfGC2(gcfit.temp = gcfit.temp, gctemp = gc_qc,
                            Y = Y_qc, norm_index = norm_index, groups = groups, 
                            offset = offset, T = T, alpha = alpha0,
                            minCountQC = minCountQC)
                fhatnew <- fhat.temp$fGC.hat
                alpha <- fhat.temp$alpha
            }
            fhatnew[fhatnew < quantile(fhatnew, 0.005)] <- quantile(
                fhatnew, 0.005)
            betahatnew <- apply((Y_qc/(fhatnew * Nmat * exp(ghat %*%
                t(hhat))))[, norm_index], 1, median)
            betahatnew[betahatnew <= 0] <- min(betahatnew[betahatnew > 0])
            bhdiff[iter] <- sum((betahatnew - betahat)^2)/length(betahat)
            fhdiff[iter] <- sum((fhatnew - fhat)^2)/length(fhat)
            if (fhdiff[iter] > min(fhdiff))
                break
            message("Iteration ", iter, "\t", "beta diff =",
                signif(bhdiff[iter], 3), "\t", "f(GC) diff =",
                signif(fhdiff[iter], 3))
            fhat <- fhatnew
            betahat <- betahatnew
            betahatmat <- matrix(nrow = nrow(Y_qc), ncol = ncol(Y_qc),
                        data = betahat, byrow = FALSE)
            L <- log(Nmat * fhat * betahatmat * alpha/2)
            logmat <- log(pmax(Y_qc, 1)) - L
            logmat <- logmat - matrix(nrow = nrow(Y_qc), ncol = ncol(Y_qc),
                        data = apply(logmat, 1, mean), byrow = FALSE)
            hhat <- svd(logmat, nu = k, nv = k)$v
            hhatnew <- hhat
            hiter <- 1
            hhdiff <- rep(Inf, maxhiter)
            while (hiter <= maxhiter) {
                for (s in seq_len(nrow(Y_qc))) {
                    temp <- try(glm(formula = Y_qc[s, norm_index] ~ 
                                hhat[norm_index, ] -
                                1, offset = L[s, norm_index],
                                family = poisson)$coefficients,
                                silent = TRUE)
                    if (is.character(temp)) {
                        temp <- lm(log(pmax(Y_qc[s, norm_index], 1)) ~
                            hhat[norm_index, ] -
                            1, offset = log(L[s, norm_index]))$coefficients
                    }
                    ghat[s, ] <- temp
                }
                # avoid overflow or underflow of the g latent factors
                ghat[is.na(ghat)] <- 0
                if (max(ghat) >= 30) {
                    ghat <- apply(ghat, 2, function(z) {
                        z[z > quantile(z, 0.995)] = min(quantile(z, 0.995), 30)
                        z
                        })
                }
                if (min(ghat) <= -30) {
                    ghat <- apply(ghat, 2, function(z) {
                        z[z < quantile(z, 0.005)] = max(quantile(z,
                        0.005), -30)
                        z
                        })
                }
                for (t in seq_len(ncol(Y_qc))) {
                    hhatnew[t, ] <- glm(formula = Y_qc[, t] ~ ghat - 1,
                        offset = L[, t], family = poisson)$coefficients
                }
                gh <- ghat %*% t(hhatnew)
                gh <- scale(gh, center = TRUE, scale = FALSE)
                hhatnew <- svd(gh, nu = k, nv = k)$v
                hhdiff[hiter] <- sum((hhatnew - hhat)^2)/length(hhat)
                message("\t\t\t", "hhat diff =",
                    signif(hhdiff[hiter], 3))
                hhat <- hhatnew
                if (hhdiff[hiter] < HHTHRESH)
                    break
                if (hiter > 10 & (rank(hhdiff))[hiter] <= 3)
                    break
                hiter <- hiter + 1
            }
        
            alphahatlist[[iter]] <- alpha
            fhatlist[[iter]] <- fhat
            betahatlist[[iter]] <- betahat
            ghatlist[[iter]] <- ghat
            hhatlist[[iter]] <- hhat
            if (bhdiff[iter] < BHTHRESH)
                break
            if (iter > 5 & bhdiff[iter] > 1)
                break
            iter <- iter + 1
        }
        optIter <- which.min(fhdiff)
        message(paste("Stop at Iteration ", optIter, ".", sep = ""))
        alpha <- alphahatlist[[optIter]]
        fhat <- fhatlist[[optIter]]
        betahat <- betahatlist[[optIter]]
        ghat <- ghatlist[[optIter]]
        hhat <- hhatlist[[optIter]]
        betahatmat <- matrix(nrow = nrow(Y_qc), ncol = ncol(Y_qc),
                        data = betahat, byrow = FALSE)
        Yhat[[ki]] <- pmax(round(fhat * Nmat * betahatmat *
                        exp(ghat %*% t(hhat)), 0), 1)
        alpha.hat[[ki]] <- alpha
        fGC.hat[[ki]] <- signif(fhat, 3)
        beta.hat[[ki]] <- signif(betahat, 3)
        h.hat[[ki]] <- signif(hhat, 3)
        g.hat[[ki]] <- signif(ghat, 3)
        Yhat.temp <- Yhat[[ki]] * alpha/2
        AIC[ki] <- 2 * sum(Y_qc * log(pmax(Yhat.temp, 1)) -
                Yhat.temp) - 2 * (length(ghat) + length(hhat))
        BIC[ki] <- 2 * sum(Y_qc * log(pmax(Yhat.temp, 1)) -
                Yhat.temp) - (length(ghat) + length(hhat)) * 
                log(length(Y_qc))
        RSS[ki] <- sum((Y_qc - Yhat.temp)^2/length(Y_qc))
        message("AIC", k, " = ", round(AIC[ki], 3))
        message("BIC", k, " = ", round(BIC[ki], 3))
        message("RSS", k, " = ", round(RSS[ki], 3), "\n")
    }
    list(Yhat = Yhat, alpha.hat = alpha.hat, fGC.hat = fGC.hat,
        beta.hat = beta.hat, g.hat = g.hat, h.hat = h.hat,
        AIC = AIC, BIC = BIC, RSS = RSS, K = K)
}

getfGCG <- function(gcfit.tempg, gctemp, Yg, offsetg, T, draw.plot = NULL,
                    alphag, minCountQC) {
    alphag <- pmax(1, round(alphag))
    if (is.null(draw.plot)) {
        draw.plot <- FALSE
    }
    loglik <- rep(NA, length(T))
    BIC <- rep(NA, length(T))
    fGCi.obj <- vector("list", length(T))
    fGCg.plot <- vector("list", length(T))
    Z.obj <- vector("list", length(T))
    vec_pi.obj <- vector("list", length(T))

    bin.filter.list = vector("list", length = ncol(gcfit.tempg))
    for (j in seq_len(ncol(gcfit.tempg))) {
        fGCi <- fitGC(gctemp, gcfit.tempg[,j])
        resid <- abs(gcfit.tempg[,j] - fGCi)
        bin.filter.list[[j]] <- which(resid > (median(resid) + 5 * 
            mad(resid)) | Yg[,j] < minCountQC)
        if (length(bin.filter.list[[j]]) == 0) {
            bin.filter.list[[j]] <- which.max(gcfit.tempg[,j])
        }
    }
    bin.filter = sort(unique(unlist(bin.filter.list)))

    for (Ti in T) {
        cat(Ti, "\t")
        if (Ti == 1) {
            fGCg <- matrix(NA, ncol = ncol(gcfit.tempg), 
                nrow = nrow(gcfit.tempg))
            Z <- matrix(nrow = nrow(gcfit.tempg), ncol = Ti,
                data = 1/Ti)
            vec_pi <- 1
            for (j in seq_len(ncol(gcfit.tempg))) {
                loe.fit.temp <- loess(gcfit.tempg[-bin.filter,j] ~
                        gctemp[-bin.filter])
                fGCi <- predict(loe.fit.temp, newdata = gctemp,
                        se = TRUE)$fit
                temp <- min(fGCi[!is.na(fGCi) & fGCi > 0])
                fGCi[fGCi <= 0 | is.na(fGCi)] <- temp
                fGCg[,j] = fGCi
            }
        }
        if (Ti >= 2) {
            Z <- matrix(nrow = nrow(gcfit.tempg), ncol = Ti, data = 0)
            mintemp <- pmin(Ti, alphag)
            Z[cbind(seq_len(nrow(Z)), mintemp)] <- 1
            vec_pi <- colSums(Z)/nrow(Z)

            fGCg <- matrix(NA, ncol = ncol(gcfit.tempg), 
                nrow = nrow(gcfit.tempg))
            for (j in seq_len(ncol(gcfit.tempg))) {
                loe.fit.temp <- loess((gcfit.tempg[,j]/(Z %*% 
                                as.matrix(seq_len(Ti)/2)))[-bin.filter] ~
                                gctemp[-bin.filter])
                fGCi <- predict(loe.fit.temp, newdata = gctemp, se = TRUE)$fit
                temp <- min(fGCi[!is.na(fGCi) & fGCi > 0])
                fGCi[fGCi <= 0 | is.na(fGCi)] <- temp
                fGCg[,j] <- fGCi
            }

            diff.GC <- Inf
            diff.Z <- Inf
            iter <- 1
            while (iter <= 3 | diff.GC > 5e-06 | diff.Z > 0.005) {
                Mtemp <- MstepG(Z, gcfit.tempg, gctemp)
                vec_pi.new <- Mtemp$vec_pig
                fGCi.new <- Mtemp$fGCg
                Z.new <- EstepG(fGCi.new, vec_pi.new, Yg, offsetg)
                diff.GC <- sum((fGCg - fGCi.new)^2)/length(fGCg)
                diff.Z <- sum((Z.new - Z)^2)/length(Z)
                vec_pi <- vec_pi.new
                Z <- Z.new
                fGCg <- fGCi.new
                iter <- iter + 1
                if (iter >= 50)
                    break
            }
        }
    
        fGCg.plot[[which(T == Ti)]] <- matrix(NA, ncol = ncol(gcfit.tempg), 
                nrow = nrow(gcfit.tempg))
        if (Ti == 1) {
            loglik.temp <- rep(NA, ncol(gcfit.tempg))
            BIC.temp <- rep(NA, ncol(gcfit.tempg))
            for (j in seq_len(ncol(gcfit.tempg))) {
                loe.fit.plot <- loess(gcfit.tempg[,j] ~ gctemp)
                fGCi.plot <- loe.fit.plot$fitted
                temp <- min(fGCi.plot[!is.na(fGCi.plot) & fGCi.plot > 0])
                fGCi.plot[fGCi.plot <= 0 | is.na(fGCi.plot)] <- temp
                df <- predict(loe.fit.plot, newdata = gctemp, se = TRUE)$df
                fGCg.plot[[which(T == Ti)]][,j] = fGCi.plot
        
                loglik.temp[j] <- sum(dpois(Yg[-bin.filter,j], 
                        lambda = (offsetg[,j] * fGCg[,j])[-bin.filter], 
                        log = TRUE))
                BIC.temp[j] <- 2 * loglik.temp[j] - (nrow(gcfit.tempg) - df) * 
                        log(nrow(gcfit.tempg))
            }
            loglik[which(T == Ti)] <- sum(loglik.temp)
            BIC[which(T == Ti)] <- sum(BIC.temp)
        } else {
            loglik.temp <- rep(NA, ncol(gcfit.tempg)) 
            BIC.temp <- rep(NA, ncol(gcfit.tempg))
                for (j in seq_len(ncol(gcfit.tempg))) {
                    loe.fit.plot <- loess((gcfit.tempg[,j]/(Z %*%
                        as.matrix(seq_len(Ti)/2))) ~ gctemp)
                    fGCi.plot <- loe.fit.plot$fitted
                    temp <- min(fGCi.plot[!is.na(fGCi.plot) & fGCi.plot > 0])
                    fGCi.plot[fGCi.plot <= 0 | is.na(fGCi.plot)] <- temp
                    df <- predict(loe.fit.plot, newdata = gctemp, se = TRUE)$df
                    fGCg.plot[[which(T == Ti)]][,j] <- fGCi.plot
        
                    loglik.temp[j] <- sum(dpois(Yg[-bin.filter,j], lambda = 
                            (offsetg[,j] * fGCg[,j] * (Z %*% 
                            as.matrix(seq_len(Ti)/2)))[-bin.filter], 
                            log = TRUE))
                    BIC.temp[j] <- 2 * loglik.temp[j] - (nrow(gcfit.tempg) - 
                            df + Ti - 1) * log(nrow(gcfit.tempg))
                }
            loglik[which(T == Ti)] <- sum(loglik.temp)
            BIC[which(T == Ti)] <- sum(BIC.temp)
        }
    
        fGCi.obj[[which(T == Ti)]] <- fGCg
        Z.obj[[which(T == Ti)]] <- Z
        vec_pi.obj[[which(T == Ti)]] <- vec_pi
    }


    for (j in seq_len(ncol(gcfit.tempg))) {
        if (draw.plot) {
            par(mfrow = c(5, 2))
            smoothScatter(gctemp[-bin.filter], gcfit.tempg[-bin.filter,j],
                    main = "Original", xlab = "GC content",
                    ylab = "Y/beta/N/exp(gxh)")
        }
    
        for (Ti in T) {
            if (draw.plot) {
                smoothScatter(gctemp[-bin.filter], gcfit.tempg[-bin.filter,j],
                    xlab = "GC content", ylab = "Y/beta/N/exp(gxh)",
                    nrpoints = 0, main = paste("T =", Ti))
                if (Ti == 1) {
                    points(gctemp[order(gctemp)], fGCg.plot[[which(T == Ti)]][
                        order(gctemp),j],
                    lty = 2, col = 2, type = "l", lwd = 2)
                    points(gctemp, gcfit.tempg[,j], cex = 0.4,
                        col = 2, pch = 16)
                } else {
                    for (k in seq_len(Ti)) {
                        points(gctemp[order(gctemp)],
                            fGCg.plot[[which(T == Ti)]][order(gctemp),j] * k/2,
                            lty = 2, col = k, type = "l", lwd = 2)
                        points(gctemp[which((round(Z.obj[[which(T == Ti)
                            ]]))[, k] == 1)],
                            (gcfit.tempg[,j])[which((round(Z.obj[[which(T == Ti)
                            ]]))[, k] == 1)],
                            cex = 0.4, col = k, pch = 16)
                    }
                }
            }
        }
    
        if (draw.plot) {
            plot(T, loglik, type = "b", xlab = "T", ylab = "loglik", 
                main = "Log likelihood")
            abline(v = which.max(loglik), lty = 2)
            plot(T, BIC, type = "b", xlab = "T", ylab = "BIC", main = "BIC")
            abline(v = which.max(BIC), lty = 2)
            par(mfrow = c(1, 1))
        }
    }
    return(list(fGCi.obj = fGCi.obj, Z.obj = Z.obj,
                vec_pi.obj = vec_pi.obj, bin.filter = bin.filter,
                loglik = loglik, BIC = BIC))
}

getfGC2 <- function(gcfit.temp, gctemp, Y, norm_index, groups, offset, T,
                    alpha, minCountQC) {
    fGC.hat <- matrix(ncol = ncol(Y), nrow = nrow(Y))
    for (G in unique(groups)) {
        cat(G, "\t")
        g <- which(groups == G)
        if (!any(is.na(match(norm_index, g))) & 
            !any(is.na(match(g, norm_index)))) {
            alpha[, g] <- 2
            for (j in g) {
                loe.fit <- loess(gcfit.temp[, j] ~ gctemp)
                fGC.hat[, j] <- loe.fit$fitted  
            }
        } else {
            fGCg <- getfGCG(gcfit.tempg = gcfit.temp[, g, drop = FALSE],
                        gctemp = gctemp, Yg = Y[, g, drop = FALSE], 
                        offsetg = offset[, g, drop = FALSE],
                        T = T, draw.plot = FALSE, 
                        alphag = apply(alpha[, g, drop = FALSE], 1, median),
                        minCountQC = minCountQC)
            if (which.max(fGCg$BIC) == 1) {
                alpha[, g] <- 2
            } else {
                alpha[, g] <- apply(fGCg$Z.obj[[which.max(fGCg$BIC)]], 
                                1, which.max)
            }
            fGC.hat[, g] <- fGCg$fGCi.obj[[which.max(fGCg$BIC)]]
        }
    }
    return(list(fGC.hat = fGC.hat, alpha = alpha))
}


EstepG <- function(fGCg, vec_pig, Yg, offsetg) {
    P <- matrix(nrow = nrow(fGCg), ncol = length(vec_pig))
    vec_pig[vec_pig == 0] <- 1e-100
    pPoisson <- vector("list", length = ncol(fGCg))
    for (j in seq_len(ncol(Yg))) {
        lambdaj <- matrix(nrow = nrow(P), ncol = ncol(P),
            data = offsetg[, j, drop = FALSE] * fGCg[, j, drop = FALSE]) * 
            matrix(nrow = nrow(P), ncol = ncol(P),
            data = seq_len(length(vec_pig))/2, byrow = TRUE)
        pPoisson[[j]] <- dpois(matrix(nrow = nrow(P), ncol = ncol(P), 
                    data = Yg[, j, drop = FALSE]),
                    lambda = lambdaj, log = TRUE)
    }
    P <- apply(simplify2array(pPoisson), c(1,2), sum) + matrix(nrow = nrow(P), 
                    ncol = ncol(P), data = log(vec_pig), byrow = TRUE)
    Zg <- matrix(nrow = nrow(fGCg), ncol = length(vec_pig))
    for (k in seq_len(length(vec_pig))) {
        Zg[, k] <- 1/(1 + apply(exp(apply(P, 2, function(x) {
            x - P[, k]
            })[, -k, drop = FALSE]), 1, sum))
    }
    return(Zg)
}



MstepG <- function(Zg, gcfit.tempg, gctemp) {
    fGCg <- matrix(NA, nrow = nrow(gcfit.tempg), ncol = ncol(gcfit.tempg))
    for (j in seq_len(ncol(gcfit.tempg))) {
        gcfit.temp <- gcfit.tempg[,j]/(Zg %*% as.matrix(seq_len(ncol(Zg))/2))
        fGCi <- fitGC(gctemp, gcfit.temp)
        fGCg[,j] <- fGCi
    }
    vec_pig <- colSums(Zg)/nrow(Zg)
    return(list(vec_pig = vec_pig, fGCg = fGCg))
}



fitGC <- function(gctemp, gcfit.temp) {
    loe.fit <- loess(gcfit.temp ~ gctemp)
    fGCi <- loe.fit$fitted
    temp <- min(fGCi[fGCi > 0])
    fGCi[fGCi <= 0] <- temp
    return(fGCi)
}


