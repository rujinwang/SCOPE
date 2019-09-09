#' @title Cross-sample segmentation
#'
#' @description SCOPE offers a cross-sample Poisson likelihood-based
#' recursive segmentation, enabling shared breakpoints across cells
#' from the same genetic background.
#'
#' @param Y raw read depth matrix after quality control procedure
#' @param Yhat normalized read depth matrix
#' @param sampname vector of sample names
#' @param ref GRanges object after quality control procedure
#' @param chr chromosome name. Make sure it is consistent with the
#'  reference genome.
#' @param mode format of returned copy numbers. Only integer mode is
#'  supported for scDNA-seq data.
#' @param max.ns a number specifying how many rounds of nested structure
#'  searching would be performed. Defalut is \code{0}.
#'
#' @return A list with components
#'     \item{poolcall}{Cross-sample CNV callings indicating
#'         shared breakpoints}
#'     \item{finalcall}{Final cross-sample segmented callset of
#'         CNVs with genotyping results}
#'     \item{image.orig}{A matrix giving logarithm of normalized
#'         z-scores}
#'     \item{image.seg}{A matrix of logarithm of estimated
#'         copy number over 2}
#'     \item{iCN}{A matrix of inferred integer copy number profiles}
#'
#' @examples
#' Yhat.sim = normObj.scopeDemo$Yhat[[which.max(normObj.scopeDemo$BIC)]]
#' segment_cs_chr1 = segmentCBScs(Y = Y_sim, Yhat = Yhat.sim,
#'                         sampname = colnames(Y_sim),
#'                         ref = ref_sim, chr = 'chr1', max.ns = 1)
#'
#' @author Rujin Wang \email{rujin@email.unc.edu}
#' @importFrom GenomicRanges GRanges
#' @importFrom IRanges IRanges RangesList Views countOverlaps
#' @importFrom GenomeInfoDb seqnames
#' @export
segmentCBScs = function(Y, Yhat, sampname, ref, chr,
    mode = "integer", max.ns = 0) {
    if (is.na(match(chr, unique(as.character(seqnames(ref)))))) {
        stop("Chromosome not found in the reference genome! Make sure that
            all chromosomes are named consistently.  \n")
    }

    stbin.flag = min(which(as.character(seqnames(ref)) == chr))
    Y = Y[which(as.character(seqnames(ref)) == chr), ]
    Yhat = Yhat[which(as.character(seqnames(ref)) == chr), ]
    ref = ref[which(as.character(seqnames(ref)) == chr)]

    poolcall = NULL
    message("Cross-sample segmenting for ", ncol(Y), " samples.")

    chpts0 = c(1, nrow(Y))
    cs.scan = compute_cs_lratio(Y, Yhat, sampname, chpts0)
    i = cs.scan$i
    j = cs.scan$j
    Z = cs.scan$Z
    init.cs.finalmat = cs.scan$finalmat

    if (!is.null(init.cs.finalmat) && nrow(init.cs.finalmat) > 0) {
        # Further cross-sample nested searching
        chpts = init.cs.finalmat[, c(1, 2)]
        if (is.null(dim(chpts))) {
            chpts = t(as.matrix(chpts))
        }

        keep_going = 1
        # number of nested searching
        num_ns = 1
        while (!all(keep_going == 0) & num_ns <= max.ns) {
            nested.output = search_cs_nested(Y, Yhat, sampname, chpts)
            keep_going = vapply(nested.output,
                function(z) z$is.nested, numeric(1))

            newchpts = NULL
            for (r in seq_len(length(keep_going))) {
                if (keep_going[r]) {
                    newchpts = rbind(newchpts,
                        nested.output[[r]]$finalmat[, c(1, 2)])
                } else {
                    newchpts = rbind(newchpts, chpts[r, ])
                }
            }
            chpts = newchpts
            num_ns = num_ns + 1
        }

        # backward compute Z after nested searching
        temp = NULL
        for (r in seq_len(nrow(chpts))) {
            idx = which(i == chpts[r, 1] & j == chpts[r, 2])
            if (length(idx) != 0) {
                temp = c(temp, Z[idx])
            } else {
                temp = c(temp, rep(NA, 4))
            }
        }
        cs.finalmat = cbind(chpts, temp)
        cs.finalmat = cs.finalmat[!is.na(cs.finalmat[, 3]), ,
            drop = FALSE]
        cs.finalmat = cs.finalmat[order(-cs.finalmat[, 3]), ,
            drop = FALSE]
    } else {
        cs.finalmat = matrix(data = c(chpts0, 0),
            nrow = 1, ncol = 3, byrow = TRUE)
    }

    cs.loglikeij <- rep(NA, nrow(cs.finalmat))
    cs.mBIC <- rep(NA, nrow(cs.finalmat))

    kappa1 = 3/2
    kappa2 = 2.27
    N = ncol(Y)
    T = nrow(Y)

    for (s in seq_len(nrow(cs.finalmat))) {
        tau <- sort(unique(c(as.vector(cs.finalmat[seq_len(s),
            c(1, 2)]), 1, nrow(Y))))
        m <- length(tau) - 2
        if (m > 0) {
            J = matrix(data = NA, nrow = m + 2, ncol = N)

            Y0 = Y
            Yhat0 = Yhat
            Y0[Y0 <= 20] = 20
            Yhat0[Yhat0 <= 20] = 20

            muhat = matrix(data = NA, nrow = m + 1, ncol = N)
            rhat = matrix(data = NA, nrow = m + 1, ncol = N)
            muhat[1, ] = round(apply(Y0[seq_len(tau[2]), ,
                drop = FALSE], 2, sum)/tau[2])
            rhat[1, ] = round(2 * (apply(Y0[seq_len(tau[2]), ,
                drop = FALSE], 2, sum)/apply(Yhat0[seq_len(tau[2]),
                , drop = FALSE], 2, sum)))
            for (r in seq_len(m)) {
                muhat[r + 1, ] = round(apply(
                    Yhat0[(tau[r + 1] + 1):tau[r + 2], ,
                    drop = FALSE], 2, sum)/(tau[r + 2] - tau[r + 1] + 1))
                rhat[r + 1, ] = round((apply(Y0[
                    (tau[r + 1] + 1):tau[r + 2], ,
                    drop = FALSE], 2, sum)/apply(Yhat0[
                    (tau[r + 1] + 1):tau[r + 2], ,
                    drop = FALSE], 2, sum)))
            }
            carriershat = rhat[2:(m + 1), ] - rhat[seq_len(m), ]
            carriershat[carriershat != 0] = 1
            J[2:(m + 1), ] = carriershat
            deltahat = muhat[2:(m + 1), ] - muhat[seq_len(m), ]
            deltahatJ = deltahat * J[2:(m + 1), ]
            if (is.null(dim(deltahatJ))) {
                deltahatJ = matrix(data = deltahat, nrow = 1,
                    ncol = length(deltahat))
            }
            deltahatJ.sq.sum = apply(deltahatJ^2, 1, sum)
            deltahatJ.sq.sum = max(deltahatJ.sq.sum, 1)

            # M and pi
            M = sum(carriershat)
            pihat = M/(N * m)

            cs.temp = create_chptsmat(cs.finalmat[seq_len(s), c(1, 2),
                drop = FALSE], chpts0)
            L = matrix(data = NA, ncol = N, nrow = nrow(cs.temp))
            for (r in seq_len(nrow(L))) {
                yact.temp = apply(Y[cs.temp[r, 1]:cs.temp[r, 2], ,
                    drop = FALSE], 2, sum)
                lambda.temp = apply(Yhat[cs.temp[r, 1]:cs.temp[r, 2], ,
                    drop = FALSE], 2, sum)
                yact.temp[lambda.temp < 20] <- 20
                lambda.temp[lambda.temp < 20] <- 20
                L[r, ] = (1 - round(2 * yact.temp/lambda.temp)/2) *
                    lambda.temp + log((round(2 * (yact.temp/lambda.temp)) +
                    1e-04)/2.0001) * yact.temp
            }
            cs.loglikeij[s] = sum(L)

            term1 = cs.loglikeij[s]
            if (M == 0) {
                term2 = 0
            } else {
                term2 = -M/2 * log(2 * cs.loglikeij[s]/M)
            }
            term3 = -log(choose(T, m))
            # term3 = - m * log(choose(T, m))
            term4 = -M/2
            term5 = -sum(log(deltahatJ.sq.sum))
            term6 = -m * (kappa1 - kappa2)
            if (pihat == 0 || pihat == 1) {
                term7 = 0
            } else {
                term7 = (M * log(pihat) + (N * m - M) * log(1 - pihat))
            }
            mbic = term1 + term2 + term3 + term4 + term5 + term6 + term7
            cs.mBIC[s] <- mbic
        } else {
            cs.mBIC[s] = 0
        }
    }
    cs.mBIC <- round(cs.mBIC, digits = 3)

    cs.finalmat <- (cbind(cs.finalmat, cs.mBIC)[seq_len(which.max(cs.mBIC)), ,
        drop = FALSE])
    poolcall <- create_chptsmat(cs.finalmat[, c(1, 2), drop = FALSE], chpts0)
    colnames(poolcall) = c("st", "end")

    finalcall = NULL
    image.orig = log(pmax(0.001, Y)/Yhat)
    image.orig[image.orig <= -2] = -2
    image.orig[image.orig > 2] = 2

    image.seg = matrix(data = NA, nrow = nrow(image.orig),
        ncol = ncol(image.orig))
    iCN = matrix(data = NA, nrow = nrow(image.orig),
        ncol = ncol(image.orig))
    colnames(image.seg) = colnames(Y)
    colnames(iCN) = colnames(Y)

    for (i in seq_len(nrow(poolcall))) {
        st_bin = poolcall[i, "st"]
        ed_bin = poolcall[i, "end"]

        yact.ind = colSums(Y[st_bin:ed_bin, , drop = FALSE])
        lambda.ind = colSums(Yhat[st_bin:ed_bin, , drop = FALSE])
        if (mode == "integer") {
            chat.ind <- round(2 * (yact.ind/lambda.ind))
        } else if (mode == "fraction") {
            chat.ind <- 2 * (yact.ind/lambda.ind)
        }
        chat.ind[chat.ind > 14] = 14

        image.seg[poolcall[i, 1]:poolcall[i, 2], ] =
            matrix(nrow = (poolcall[i, 2] - poolcall[i, 1] + 1),
            ncol = ncol(image.seg),
            data = log(chat.ind/2), byrow = TRUE)
        iCN[poolcall[i, 1]:poolcall[i, 2], ] =
            matrix(nrow = (poolcall[i, 2] - poolcall[i, 1] + 1),
            ncol = ncol(iCN),
            data = chat.ind, byrow = TRUE)

        temp = cbind(colnames(Y), rep(st_bin, ncol(Y)),
            rep(ed_bin, ncol(Y)), chat.ind)
        temp = temp[chat.ind != 2, ]
        finalcall = rbind(finalcall, temp)
        rownames(finalcall) = NULL
        colnames(finalcall) = c("sample_name", "st_bin",
            "ed_bin", "cnv_no")
    }

    finalcall = as.data.frame(finalcall)

    finalcall$st_bin = as.numeric(paste(finalcall$st_bin))
    finalcall$ed_bin = as.numeric(paste(finalcall$ed_bin))
    finalcall$cnv_no = as.numeric(paste(finalcall$cnv_no))


    poolcall[, "st"] = poolcall[, "st"] + stbin.flag - 1
    poolcall[, "end"] = poolcall[, "end"] + stbin.flag - 1
    finalcall[, "st_bin"] = finalcall[, "st_bin"] + stbin.flag - 1
    finalcall[, "ed_bin"] = finalcall[, "ed_bin"] + stbin.flag - 1
    return(list(poolcall = poolcall, finalcall = finalcall,
        image.orig = image.orig, image.seg = image.seg, iCN = iCN))
}



create_chptsmat = function(mat1, st_end) {
    st = st_end[1]
    end = st_end[2]
    mat1 = mat1[order(mat1[, 1]), , drop = FALSE]
    if (mat1[1, 1] != st) {
        newchptsmat = matrix(data = c(st, mat1[1, 1] - 1, mat1[1, ]),
            ncol = 2, byrow = TRUE)
    } else {
        newchptsmat = t(as.matrix(mat1[1, ]))
    }
    if (nrow(mat1) > 1) {
        for (r in 2:nrow(mat1)) {
            if (mat1[r, 1] != mat1[r - 1, 2] + 1) {
                newchptsmat = rbind(newchptsmat, matrix(data =
                    c(mat1[r - 1, 2] + 1, mat1[r, 1] - 1, mat1[r, ]),
                    ncol = 2, byrow = TRUE))
            } else {
                newchptsmat = rbind(newchptsmat, matrix(data = mat1[r, ],
                    ncol = 2, byrow = TRUE))
            }
        }
    }
    if (mat1[nrow(mat1), 2] != end) {
        newchptsmat = rbind(newchptsmat,
            matrix(data = c(mat1[nrow(mat1), 2] + 1, end), nrow = 1, ncol = 2))
    }
    return(newchptsmat)
}



compute_cs_lratio = function(Y, Yhat, sampname, this.chpts, msgprint = TRUE) {
    Z0 = NULL
    if (this.chpts[2] - this.chpts[1] < 2) {
        return(list(i = NA, j = NA, Z = NA,
            finalmat = t(as.matrix(this.chpts)), is.nested = 0))
    } else {
        lmax = nrow(Y) - 1
        for (sampno in seq_len(ncol(Y))) {
            if (msgprint) {
                message("Calculating scan statistic for sample ",
                    sampno, ": ", sampname[sampno], ".")
            }

            y = Y[, sampno]
            yhat = Yhat[, sampno]

            if (any(yhat < 10)) {
                if (any(y[yhat < 10] < 10)) {
                    y[yhat < 10] <- 10
                }
                yhat[yhat < 10] <- 10
            }

            ysum = sum(y)
            yhatsum = sum(yhat)
            num = length(y)
            y = c(y, rep(0, lmax))
            yhat = c(yhat, rep(0, lmax))
            i = rep(seq_len(num), rep(lmax, num))
            j = rep(seq_len(lmax), num) + i
            yact = rep(0, length(i))
            lambda = rep(0, length(i))
            for (k in seq_len(num)) {
                yact[(lmax * k - (lmax - 1)):(lmax * k)] <- cumsum(
                    y[k:(k + lmax)])[-1]
                lambda[(lmax * k - (lmax - 1)):(lmax * k)] <- cumsum(
                    yhat[k:(k + lmax)])[-1]
            }
            i <- i[j <= num]
            yact <- yact[j <= num]
            lambda <- lambda[j <= num]
            j <- j[j <= num]

            # integer mode
            chat <- round(2 * (yact/lambda))
            idx.noLeft = which(i == this.chpts[1])
            idx.noRight = which(j == this.chpts[2])
            chat.L = rep(NA, length(chat))
            chat.R = rep(NA, length(chat))
            yact.L = rep(NA, length(chat))
            yact.R = rep(NA, length(chat))
            lambda.L = rep(NA, length(chat))
            lambda.R = rep(NA, length(chat))
            for (x in seq_len(length(chat))) {
                if (x %in% idx.noLeft) {
                    yact.L[x] = 0
                    lambda.L[x] = 0
                } else {
                    yact.L[x] = sum(y[seq_len((i[x] - 1))])
                    lambda.L[x] = sum(yhat[seq_len((i[x] - 1))])
                }

                if (x %in% idx.noRight) {
                    yact.R[x] = 0
                    lambda.R[x] = 0
                } else {
                    yact.R[x] = sum(y[(j[x] + 1):length(y)])
                    lambda.R[x] = sum(yhat[(j[x] + 1):length(y)])
                }
            }
            chat.L <- round(2 * (yact.L/lambda.L))
            chat.R <- round(2 * (yact.R/lambda.R))
            chat.L[which(is.na(chat.L))] <- 2
            chat.R[which(is.na(chat.R))] <- 2

            lratio.C = (1 - chat/2) * lambda +
                log((chat + 1e-04)/2.0001) * yact
            lratio.L = (1 - chat.L/2) * lambda.L +
                log((chat.L + 1e-04)/2.0001) * yact.L
            lratio.R = (1 - chat.R/2) * lambda.R +
                log((chat.R + 1e-04)/2.0001) * yact.R
            lratio = lratio.C + lratio.L + lratio.R
            chat[chat > 14] <- 14

            Z0 = cbind(Z0, lratio)
        }
        Z = rowSums(Z0)

        winlag = this.chpts[1] - 1
        i = i + winlag
        j = j + winlag

        if (sum(Z > 0) > 0) {
            if (sum(Z > 0) >= 2) {
                finalmat <- (cbind(i, j, Z))[Z > 0, ]
                finalmat <- finalmat[order(-finalmat[, 3]), ]
                s <- 1
                while (s <= (nrow(finalmat))) {
                    rowstart <- finalmat[s, 1]
                    rowend <- finalmat[s, 2]
                    rowsel <- (finalmat[, 1] <= rowend &
                        finalmat[, 2] >= rowstart)
                    rowsel[s] <- FALSE
                    finalmat <- finalmat[!rowsel, ]
                    if (is.vector(finalmat)) {
                        finalmat <- t(as.matrix(finalmat))
                    }
                    s <- s + 1
                }
            }
            if (sum(Z > 0) == 1) {
                finalmat <- (cbind(i, j, Z))[Z > 0, ]
                finalmat <- t(as.matrix(finalmat))
            }
            finalmat <- finalmat[finalmat[, 2] - finalmat[, 1] > 1, ,
                drop = FALSE]
            finalmat <- finalmat[finalmat[, 3] > 10, , drop = FALSE]
            finalmat <- round(finalmat, digits = 3)

            if (nrow(finalmat) == 0) {
                is.nested = 0
            } else if (nrow(finalmat) == 1 & finalmat[1, 1] ==
                this.chpts[1] & finalmat[1, 2] == this.chpts[2]) {
                is.nested = 0
            } else if (nrow(finalmat) > 2 &
                length(unique(finalmat[, 3])) == 1) {
                is.nested = 0
            } else if (max(finalmat[, 2] - finalmat[, 1]) <= 5) {
                is.nested = 0
            } else if (nrow(finalmat) == 1 & (this.chpts[2] -
                this.chpts[1]) -
                (finalmat[1, 2] - finalmat[1, 1]) == 1) {
                is.nested = 0
            } else {
                is.nested = 1
            }
        } else {
            finalmat = NULL
            is.nested = 0
        }

        # Avoid losing the boundaries
        if (!is.null(finalmat) && nrow(finalmat) > 0) {
            min.st = min(finalmat[, c(1, 2)])
            if (this.chpts[1] < min.st) {
                idx = which(i == this.chpts[1] & j == (min.st - 1))
                idx2 = which(finalmat[, 1] == min.st)
                if (length(idx) != 0) {
                    Z.temp = Z[idx]
                    if (Z.temp > 0) {
                        finalmat[idx2, 1] = this.chpts[1]
                    } else {
                        finalmat = rbind(c(this.chpts[1], min.st - 1, 0),
                        finalmat)
                    }
                } else {
                    finalmat[idx2, 1] = this.chpts[1]
                }
            }
            max.ed = max(finalmat[, c(1, 2)])
            if (this.chpts[2] > max.ed) {
                idx = which(i == (max.ed + 1) & j == this.chpts[2])
                idx2 = which(finalmat[, 2] == max.ed)
                if (length(idx) != 0) {
                    Z.temp = Z[idx]
                    if (Z.temp > 0) {
                        finalmat[idx2, 2] = this.chpts[2]
                    } else {
                        finalmat = rbind(finalmat, c(max.ed + 1,
                            this.chpts[2], 0))
                    }
                } else {
                    finalmat[idx2, 2] = this.chpts[2]
                }
            }
        }

        if (!is.null(finalmat) && nrow(finalmat) > 1) {
            for (r in seq_len(nrow(finalmat))) {
                singlepoint = finalmat[r, 1] - 1
                sg.idx = which(finalmat[, 2] == singlepoint - 1)
                if (length(sg.idx) != 0 && r < sg.idx) {
                    finalmat[sg.idx, 2] = singlepoint
                }
            }
            for (r in seq_len(nrow(finalmat))) {
                singlepoint = finalmat[r, 2] + 1
                sg.idx = which(finalmat[, 1] == singlepoint + 1)
                if (length(sg.idx) != 0 && r < sg.idx) {
                    finalmat[sg.idx, 1] = singlepoint
                }
            }
        }

        if (!is.null(finalmat) && nrow(finalmat) > 1) {
            finalmat = create_chptsmat(finalmat[, c(1, 2),
                drop = FALSE], this.chpts)
        }
        return(list(i = i, j = j, Z = Z, finalmat = finalmat,
            is.nested = is.nested))
    }
}



search_cs_nested = function(Y, Yhat, sampname, chptsmat) {
    message("Performing cross-sample nested search... \n")
    nest.list = vector("list", nrow(chptsmat))
    for (r in seq_len(nrow(chptsmat))) {
        nest.list[[r]] = compute_cs_lratio(Y[chptsmat[r, 1]:chptsmat[r, 2], ,
            drop = FALSE], Yhat[chptsmat[r, 1]:chptsmat[r,
            2], , drop = FALSE], sampname, chptsmat[r, ], msgprint = FALSE)
    }
    return(nest.list)
}




