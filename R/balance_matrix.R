#' Balance matrix
#'
#' Balances a matrix so that the columns and/or rows add up
#' to a certain vector.
#' @param Y Matrix to be balanced.
#' @param v (optional) Desired sum of columns.
#' @param h (optional) Desired sum of rows.
#' @param allow_negative Are negative entries in the balanced matrix allowed?
#' @details
#' Balancing is done according to the criteria of minimum sum
#' of squares.
#'
#' If neither \code{v} nor \code{h} is given, the same matrix will be
#' returned. If only one of them is given, only that axis will be
#' balanced.
#' @returns A list containing \code{X}, the balanced matrix, and the 2-norm
#' of the change the balancing did to the original matrix.
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
#' @import CVXR
#' @export
balance_matrix <- function(Y, v = NULL, h = NULL, allow_negative = TRUE) {
  if (!methods::is(Y, "matrix")) {
    Y <- as.matrix(Y)
  }
  require(CVXR)
  n <- ncol(Y)
  m <- nrow(Y)
  X <- Variable(m,n)
  if (is.null(h) & is.null(v)) {
    result <- list(X = Y,
                   norm_of_change = 0)
    return(result)
  }
  if (!is.null(h) & !is.null(v)) {
    ok <- dplyr::near(sum(h), sum(v), tol = 1e-8)
    if (!ok) {
      stop("sum(v) != sum(h) so balancing is infeasible!")
    }
    cons <- list(sum_entries(X, 1) == h,
                 sum_entries(X, 2) == v,
                 X >= 0)
  } else if (is.null(h)){
    cons <- list(sum_entries(X, 2) == v,
                 X >= 0)
  } else if (is.null(v)){
    cons <- list(sum_entries(X, 1) == h,
                 X >= 0)
  }
  obj <- Minimize(sum_squares(X-Y))
  p <- Problem(obj, cons)
  sol <- CVXR::solve(p)
  if (sol$status != "optimal") {
    stop(paste("Optimal solution not found. Solution status:", sol$status))
  } else {
    X <- sol$getValue(X) |> round(8)
  }
  result <- list(X = X,
                 norm_of_change = norm(X-Y,"2"))
  return(result)
}

#' Make non-negative
#'
#' Modifies as little as possible the entries of a matrix
#' in order to make them non-negative, keeping row and column totals unchanged.
#' @param Y Matrix to be positivized.
#' @param allowSlack Can colSums and rowSums be modified?
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
make_non_negative <- function(Y, allowSlack = FALSE) {
  if (!methods::is(Y, "matrix")) {
    Y <- as.matrix(Y)
  }
  require(CVXR)
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
    X <- sol$getValue(X) |> round(8)
    return(X)
  }
}
