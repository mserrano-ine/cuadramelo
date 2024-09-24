#' Modify matrix by blocks
#'
#' Applies a function to a matrix by horizontal or vertical blocks.
#' @param Y Matrix
#' @param MARGIN Blocks are distributed: 1 horizontally, 2 vertically.
#' @param L Number of lines of the block.
#' @param FUN Funtion to apply to the block.
#' @param ... Arguments to be passed to FUN.
#' @returns A matrix.
#' @export
apply_by_block <- function(Y, MARGIN, L, FUN, ...){
  Z <- Y
  m <- ifelse(MARGIN == 1, ncol(Y), nrow(Y))
  if (m %% L != 0) {
    stop("L does not divide the length of the chosen MARGIN")
  }
  k <- floor(m/L)
  for (i in 1:k) {
    if (MARGIN == 1) {
      X <- Y[, (L*(i-1)+1):(L*i)]
      Z[, (L*(i-1)+1):(L*i)] <- FUN(X, ...)
    } else {
      X <- Y[(L*(i-1)+1):(L*i), ]
      Z[(L*(i-1)+1):(L*i), ] <- FUN(X, ...)
    }
  }
  return(Z)
}

#' Round matrix by blocks
#'
#' Applies \code{round_matrix()} to the blocks using \code{apply_by_block()}.
#' @param Y Matrix
#' @param MARGIN 1 horizontal, 2 vertical.
#' @param L Number of lines of the block.
#' @param digits Number of decimal places to be rounded to.
#' @export
round_by_blocks <- function(Y, MARGIN, L, digits = 0) {
  Z <- apply_by_block(Y, MARGIN, L, round_matrix, digits = digits)
  return(Z)
}
