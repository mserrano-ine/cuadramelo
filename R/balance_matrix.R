#' Stack matrix columns
#'
#' Stack the columns of a matrix into a vector.
#' @param X A numerical matrix.
#' @returns A vector whose i*j-th element is the i-th
#' element of the j-th column.
stack.matrix <- function(X) {
  Y <- X |> as.data.frame() |> stack()
  return(Y$values)
}

#' Balance matrix
#'
#' Balances a matrix so that the columns and/or rows add up
#' to a certain vector.
#' @param Y Matrix to be balanced.
#' @param v (optional) Desired sum of columns.
#' @param h (optional) Desired sum of rows.
#' @details
#' Balancing is done according to the criteria of minimum sum
#' of squares.
#'
#' If neither \code{v} nor \code{h} is given, the same matrix will be
#' returned. If only one of them is given, only that axis will be
#' balanced.
#' @returns A list containing \code{X}, the balanced matrix, and a
#' measure of the change the balancing did to the original matrix.
#' @examples
#' set.seed(2)
#' Y <- rnorm(3*5) |> matrix(3,5) |> round(3)
#' v <- c( 0.876, -1.078, 3.452, 0.261, 1.349)
#' h <- c(-1.851, 0.243, 6.468)
#' Y
#' X1 <- balance_matrix(Y, v, h)$X
#' h
#' rowSums(X1)
#' v
#' colSums(X1)
#' print(paste("The change has norm2 =", balance_matrix(Y, v, h)$change))
#' X2 <- balance_matrix(Y, v = v)
#' v
#' colSums(X2)
#' X3 <- balance_matrix(Y, h = h)
#' h
#' rowSums(X3)
#' @importFrom dplyr near
#' @importFrom MASS ginv
#' @export
balance_matrix <- function(Y, v = NULL, h = NULL) {
  y <- stack.matrix(Y)
  n <- ncol(Y)
  m <- nrow(Y)
  a <- rep(1,m) # vector de agregación vertical
  b <- t(rep(1,n)) # vector de agregación horizontal
  C <- diag(n) %x% t(a) # verticales
  R <- b %x% diag(m) # horizontales
  z <- c(v,h)
  if (!is.null(h) & !is.null(v)) {
    ok <- dplyr::near(sum(h), sum(v), tol = 1e-8)
    if (!ok) {
      stop("sum(v) != sum(h) so balancing is infeasible!")
    }
    H <- rbind(C,
               R)
  } else if (is.null(h)){
    H <- C
  } else {
    H <- R
  }
  G <- rbind(cbind(diag(n*m), t(H)),
             cbind(H,diag(0,nrow(H))))
  invG <- MASS::ginv(G)
  w <- c(y,z)
  x <- (invG %*% w)[1:(m*n),]
  X <- matrix(x, m, n)
  result <- list(X = X,
                 change = norm(X-Y,"2"))
  return(result)
}

#' Make non-negative
#'
#' Modifies as little as possible the entries of a matrix
#' in order to make them non-negative, keeping row and column totals unchanged.
#' @param X Matrix to be positivized.
#' @returns A non-negative matrix, except if it is imposible to balance the
#' matrix.
#' @examples
#' Y <- c(1,2,-1,1,
#'        2,2,3,1,
#'        1,1,-2,3) |>
#'        matrix(nrow = 3)
#' X <- make_non_negative(Y)
#' Y
#' X |> round(2)
#' rowSums(Y)
#' rowSums(X)
#' colSums(Y)
#' colSums(X)
#' set.seed(2)
#' Y <- rnorm(3*5) |> matrix(3,5) |> round(3)
#' make_non_negative(Y)
#' make_non_negative(Y, allowSlack = T) |> round()
#' @import CVXR
#' @export
make_non_negative <- function(Y, allowSack = FALSE) {
  library(CVXR)
  n <- ncol(Y)
  m <- nrow(Y)
  X <- Variable(m,n)
  v <- colSums(Y)
  h <- rowSums(Y)
  slackH <- Variable(m)
  slackV <- Variable(n)
  cons <- list(sum_entries(X, 1) == h + as.numeric(allowSlack)*slackH,
               sum_entries(X, 2) == v + as.numeric(allowSlack)*slackV,
               X >= 0)
  obj <- Minimize(sum_squares(X-Y)+sum_squares(slackH)+sum_squares(slackV))
  p <- Problem(obj, cons)
  sol <- CVXR::solve(p)
  if (sol$status != "optimal") {
    stop(paste("Optimal solution not found. Solution status:", sol$status))
  } else {
    return(sol$getValue(X))
  }
}
