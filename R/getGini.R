#' @title Compute Gini coefficients for single cells
#'
#' @description Gini index is defined as two times the area
#' between the Lorenz curve and the diagonal.
#'
#' @param Y raw read depth matrix after quality control procedure
#'
#' @return
#'   \item{Gini}{Vector of Gini coefficients for single cells
#'     from scDNA-seq}
#'
#' @examples
#' Gini = getGini(Y_sim)
#'
#' @author Rujin Wang \email{rujin@email.unc.edu}
#' @import DescTools
#' @export
getGini = function(Y) {
    Gini = rep(NA, ncol(Y))
    for (i in seq_len(ncol(Y))) {
        y = sort(Y[, i])
        x = c(0, seq_len(length(y))/length(y))
        y = c(0, cumsum(y)/sum(y))
        Gini[i] = 2 * round(0.5 - AUC(x, y), 4)
    }
    return(Gini)
}

