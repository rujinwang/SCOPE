#' @title Plot post-segmentation copy number profiles of integer values
#'
#' @description Show heatmap of inferred integer copy-number profiles
#' by SCOPE with cells clustered by hierarchical clustering
#'
#' @usage
#' iCNplot(iCNmat, ref, Gini, annotation = NULL,
#'         plot.dendrogram = TRUE, filename)
#' @param iCNmat inferred integer copy-number matrix by SCOPE,
#' with each column being a cell and each row being a genomic bin
#' @param ref GRanges object after quality control procedure
#' @param Gini vector of Gini coefficients for each cell,
#' with the same order as that of cells in columns of \code{iCNmat}
#' @param annotation vector of annotation for each cell,
#' with the same order as that of cells in columns of \code{iCNmat}.
#' Default is \code{NULL}.
#' @param plot.dendrogram logical, whether to plot the dendrogram.
#' Default is \code{TRUE}.
#' @param filename name of the output png file
#'
#' @return png file with integer copy-number profiles across single cells
#' with specified annotations
#'
#' @examples
#' Gini = getGini(Y_sim)
#' iCNplot(iCNmat = iCN_sim,
#'         ref = ref_sim,
#'         Gini = Gini,
#'         filename = 'iCNplotdemo')
#'
#' @author Rujin Wang \email{rujin@email.unc.edu}
#' @import graphics stats
#' @importFrom GenomeInfoDb seqnames
#' @importFrom gplots colorpanel
#' @importFrom RColorBrewer brewer.pal
#' @importFrom grDevices png dev.off
#' @export
iCNplot = function(iCNmat, ref, Gini, annotation = NULL,
    plot.dendrogram = TRUE, filename) {
    smart_image = function(mat, ...) {
        image(t(mat[rev(seq(nrow(mat))), ]), ...)
    }
    hm_col = c("#2166AC", "#92C5DE", "#FDFDFD", "#FDDBC7",
        "#F4A582", "#D6604D", "#B2182B", "#67001F")

    if (!is.matrix(iCNmat)) {
        stop("Invalid plot object: must be an integer matrix! \n")
    }
    if (length(ref) != nrow(iCNmat)) {
        stop("Invalid GRanges object: length of ref and # of
            rows in iCNmat must be the same")
    }
    if (!is.null(annotation)) {
        if (!is.null(dim(annotation))) {
            stop("Invalid annotation object: has to be a vector or
                factor with the same # of cells as that of iCNmat")
        }
        if (length(annotation) != ncol(iCNmat)) {
            stop("Invalid annotation object: length of annotation and
                # of cells in iCNmat must be the same")
        }
    }
    if (length(Gini) != ncol(iCNmat)) {
        stop("Invalid Gini object: length of Gini coefficient and
            # of cells in iCNmat must be the same")
    }

    # page setup
    if (is.null(annotation)) {
        if (plot.dendrogram) {
            mm = matrix(c(0, 0, 4, 0,
                2, 3, 1, 5,
                2, 3, 1, 0,
                2, 3, 1, 6), nrow = 4, byrow = TRUE)
            mh = c(2, 20, 20, 20)
            mh = mh/sum(mh)
            mw = c(0.25, 0.1, 5, 0.5)
            mw = mw/sum(mw)
        } else {
            mm = matrix(c(0, 0, 3, 0,
                0, 2, 1, 4,
                0, 2, 1, 0,
                0, 2, 1, 5), nrow = 4, byrow = TRUE)
            mh = c(2, 20, 20, 20)
            mh = mh/sum(mh)
            mw = c(0.25, 0.1, 5, 0.5)
            mw = mw/sum(mw)
        }
    } else {
        if (plot.dendrogram) {
            mm = matrix(c(0, 0, 0, 5, 0,
                2, 3, 4, 1, 6,
                2, 3, 4, 1, 7,
                2, 3, 4, 1, 8), nrow = 4, byrow = TRUE)
            mh = c(2, 20, 20, 20)
            mh = mh/sum(mh)
            mw = c(0.25, 0.1, 0.1, 5, 0.5)
            mw = mw/sum(mw)
        } else {
            mm = matrix(c(0, 0, 0, 4, 0,
                0, 2, 3, 1, 5,
                0, 2, 3, 1, 6,
                0, 2, 3, 1, 7), nrow = 4, byrow = TRUE)
            mh = c(2, 20, 20, 20)
            mh = mh/sum(mh)
            mw = c(0.25, 0.1, 0.1, 5, 0.5)
            mw = mw/sum(mw)
        }
    }
    png(paste0(filename, ".png"), width = 2500,
        height = 1600, pointsize = 25)
    layout(mm, widths = mw, heights = mh)
    par(mar = rep(0, 4))

    iCNmat = round(iCNmat)
    if (!is.null(annotation)) {
        annotation = as.factor(annotation)
    }

    chr.pos = rep(NA, length(unique(seqnames(ref))))
    for (chri in seq_len(22)) {
        chr.pos[chri] = length(ref[which(as.character(
            seqnames(ref)) == paste0("chr", chri))])
    }
    chr.pos = cumsum(chr.pos)
    xpos = round(c(0, chr.pos[seq_len(21)]) +
        (chr.pos - c(0, chr.pos[seq_len(21)]))/2)

    # 1) iCN heatmap
    dat = t(iCNmat)
    dat[dat >= 7] = 7
    dat[dat <= 0] = 0
    iCNtab = as.numeric(names(table(dat)))
    rclust = hclust(dist(dat))
    dat = dat[rclust$order, ]

    smart_image(dat, col = hm_col[iCNtab + 1], xaxs = "i",
        yaxs = "i", axes = FALSE)
    abline(v = 0, lwd = 2)
    for (i in seq_len(length(chr.pos))) {
        abline(v = chr.pos[i]/length(ref), lwd = 2)
    }

    # 2) hclust
    if (plot.dendrogram) {
        plot(rev(as.dendrogram(rclust)), leaflab = "none",
            horiz = TRUE, axes = FALSE, yaxs = "i")
    }

    # 3) Gini annotation
    anno.Gini = matrix(Gini[rclust$order], nrow = nrow(dat),
        ncol = 1)
    col.Gini = gplots::colorpanel(50, "#F7FBFF", "#084594")
    smart_image(anno.Gini, col = col.Gini, xaxs = "i",
        yaxs = "i", axes = FALSE)

    if (!is.null(annotation)) {
        # 4) Customized annotation
        anno.level = levels(annotation)
        anno.mat = matrix(match(annotation[rclust$order], anno.level),
            nrow = nrow(dat), ncol = 1)
        col.anno = brewer.pal(n = 12, name = "Set3")[
            sort(unique(match(annotation[rclust$order], anno.level)))]
        smart_image(anno.mat, col = col.anno, xaxs = "i",
            yaxs = "i", axes = FALSE)
    }

    # 5) chromosome
    anno.chrom = NULL
    for (i in seq_len(22)) {
        if (i%%2 == 1) {
            temp = matrix(rep(1, length(which(as.character(
                seqnames(ref)) == paste0("chr", i)))), nrow = 1)
            anno.chrom = cbind(anno.chrom, temp)
        } else {
            temp = matrix(rep(2, length(which(as.character(
                seqnames(ref)) == paste0("chr", i)))), nrow = 1)
            anno.chrom = cbind(anno.chrom, temp)
        }
    }
    image(t(anno.chrom), col = c("gray", "black"),
        xaxs = "i", yaxs = "i", axes = FALSE)
    pos.text = xpos/length(ref)
    text(pos.text, 0.2, seq(22), col = c("black", "grey"), cex = 1.5)

    # 6) Gini legend
    par(mar = c(2, 2, 2, 4))
    image(1, seq_len(length(brewer.pal(n = 8, name = "Blues"))),
        t(as.matrix(seq_len(length(brewer.pal(n = 8, name = "Blues"))))),
        col = brewer.pal(n = 8, name = "Blues"), xlab = "", ylab = "",
        xaxt = "n", yaxt = "n", bty = "n")
    axis(4, at = c(1, length(brewer.pal(n = 8, name = "Blues"))),
        labels = round(c(min(Gini), max(Gini)), 2), col.ticks = "white",
        col = NA, lwd.ticks = 0, cex.axis = 1.5, las = 2, font = 2)
    title("Gini", cex.main = 1.5)

    if (!is.null(annotation)) {
        plot(0, 0, type = "n", axes = FALSE)
        legend("center", legend = sort(unique(annotation)),
            col = col.anno, pch = 15, bty = "n", cex = 1.5)
    }

    # 8) iCN legend
    par(mar = c(2, 2, 2, 4))
    image(1, seq_len(length(hm_col[iCNtab + 1])),
        t(as.matrix(seq_len(length(hm_col[iCNtab + 1])))),
        col = hm_col[iCNtab + 1], xlab = "", ylab = "",
        xaxt = "n", yaxt = "n", bty = "n")
    axis(4, at = seq_len(length(hm_col[iCNtab + 1])),
        labels = c(0:7)[iCNtab + 1], col.ticks = "white", col = NA,
        lwd.ticks = 0, cex.axis = 1.5, las = 2, font = 2)
    title("integer CN", cex.main = 1.5)
    dev.off()
}







