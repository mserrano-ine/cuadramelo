#' Round a matrix
#'
#' Returns an integer matrix that preserves the rounded colSums and rowSums.
#' @param Y A matrix.
#' @returns The rounded matrix.
#' @examples
#' set.seed(2)
#' Y <- rnorm(3*5) |> matrix(3,5) |> round(3)
#' X <- round_matrix(Y)$X
#' Y
#' X
#' colSums(Y) |> round()
#' colSums(X)
#' rowSums(Y) |> round()
#' rowSums(X)
#' @import CVXR
#' @export
round_matrix <- function(Y) {
  require(CVXR)
  if (!methods::is(Y, "matrix")) {
    Y <- as.matrix(Y)
  }
  require(CVXR)
  n <- ncol(Y)
  m <- nrow(Y)
  X <- Variable(m,n, integer = TRUE)
  S <- Variable(m,n, integer = TRUE)
  Z <- X + S
  v <- colSums(Y) |> round()
  h <- rowSums(Y) |> round()
  ok <- dplyr::near(sum(h), sum(v), tol = 1e-8)
  if (!ok) {
    stop("sum(v) != sum(h) so balancing is infeasible!")
  }
  cons <- list(sum_entries(Z, 1) == h,
               sum_entries(Z, 2) == v)
  obj <- Minimize(sum_squares(X-Y)+sum_squares(S))
  p <- Problem(obj, cons)
  sol <- CVXR::solve(p)
  if (sol$status != "optimal") {
    stop(paste("Optimal solution not found. Solution status:", sol$status))
  } else {
    X <- sol$getValue(X) |> round(8)
    S <- sol$getValue(S) |> round(8)
  }
  if (norm(S,"2")>0) {
    warning("Rounding could not be done exactly. A slack was needed.")
  }
  result <- list(X = X,
                 S = S,
                 norm_of_change = norm(X-Y,"2"))
  return(result)
}
