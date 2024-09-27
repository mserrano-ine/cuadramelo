#' Modify matrix by blocks
#'
#' Applies a function to a matrix by horizontal or vertical blocks.
#' @param Y Matrix
#' @param layout Blocks are distributed: 1 horizontally, 2 vertically.
#' @param L Number of lines of the block.
#' @param FUN Funtion to apply to the block.
#' @param ... Arguments to be passed to FUN.
#' @returns A matrix.
#' @export
apply_by_block <- function(Y, layout, L, FUN, ...){
  Z <- Y
  m <- ifelse(layout == 1, ncol(Y), nrow(Y))
  if (m %% L != 0) {
    stop("L does not divide the length of the chosen MARGIN")
  }
  k <- floor(m/L)
  for (i in 1:k) {
    if (layout == 1) {
      X <- Y[, (L*(i-1)+1):(L*i)]
      Z[, (L*(i-1)+1):(L*i)] <- FUN(X, ...)
    } else if (layout == 2) {
      X <- Y[(L*(i-1)+1):(L*i), ]
      Z[(L*(i-1)+1):(L*i), ] <- FUN(X, ...)
    } else {
      stop("layout must be 1 or 2!")
    }
  }
  return(Z)
}

#' Round matrix by blocks
#'
#' Applies \code{round_matrix()} to equally-sized blocks that partition the
#' matrix either vertically or horizontally.
#' @param Y Matrix.
#' @param layout The blocks are distributed: 1 horizontally, 2 vertically.
#' @param L Number of lines that a block encompasses.
#' @param digits Number of decimal places to be rounded to.
#' @param MARGIN_BLOCK For each block
#' \itemize{
#'  \item{0} Preserves the rounded colSums and rowSums.
#'  \item{1} Preserves the rounded rowSums independently of each other.
#'  \item{2} Preserves the rounded colSums independently of each other.
#' }
#' @examples
#' set.seed(10)
#' Y <- (rnorm(32)*10) |> matrix(ncol = 2) |> round(3)
#' X <- round_by_blocks(Y, 2, 4)
#' U <- Y[5:8,] |> round_matrix()
#' X[5:8,] - U
#' @export
round_by_blocks <- function(Y, layout, L, digits = 0, MARGIN_BLOCK = 0) {
  Z <- apply_by_block(Y, layout, L, round_matrix, digits = digits, MARGIN = MARGIN_BLOCK)
  return(Z)
}

#' Balance matrix by blocks
#'
#' Applies \code{balance_matrix()} to equally-sized blocks that partition the
#' matrix either vertically or horizontally.
#' @param Y Matrix to be balanced.
#' @param col_totals Desired colSums for each block. See details.
#' @param row_totals Desired rowSums for each block. See details.
#' @param layout The blocks are distributed: 1 horizontally, 2 vertically.
#' @param L Number of lines that a block encompasses.
#' @details
#' When Y is composed of **vertically** stacked blocks, col_totals must be
#' a matrix whose rows are the colSums for each block, and row_totals
#' just a (vertical) vector.
#'
#' When Y is composed of blocks arraged **horizontally**, col_totals is a
#' (horizontal) vector, and row_totals is a matrix whose columns are the rowSums
#' for each block.
#'
#' @examples
#' set.seed(10)
#' Y <- (rnorm(32)*10) |> matrix(ncol = 2) |> round(3)
#' v <- aggregate(Y, by = list(rep(1:4, times = rep(4,4))), FUN = sum)[, -1] |>
#'   round() |> as.matrix()
#' X <- balance_by_blocks(Y, v, layout = 2, L = 4)
#' U <- Y[5:8,] |> balance_matrix(v[2,])
#' X[5:8,] - U
#' @export
balance_by_blocks <- function(Y, col_totals = NULL, row_totals = NULL,
                              layout, L){
  Z <- Y
  m <- ifelse(layout == 1, ncol(Y), nrow(Y))
  if (m %% L != 0) {
    stop("L does not divide the length of the chosen MARGIN")
  }
  k <- floor(m/L)
  for (i in 1:k) {
    if (layout == 1) {
      Z <- balance_by_blocks(t(Y), row_totals, col_totals,
                             layout, L) |> t()
    } else if (layout == 2) {
      row_totals <- as.vector(row_totals)
      X <- Y[(L*(i-1)+1):(L*i), ]
      v <- col_totals[i, ] |> as.vector()
      h <- row_totals[(L*(i-1)+1):(L*i)]
      Z[(L*(i-1)+1):(L*i), ] <- balance_matrix(X, col_totals = v, row_totals = h)
    } else {
      stop("layout must be 1 or 2!")
    }
  }
  return(Z)
}
