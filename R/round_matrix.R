round_matrix_bivariate <- function(Y, digits=0) {
  mat_up <- Y * (10**digits)
  # Suma redondeadas de filas y columnas originales
  original_row_sums <- round(rowSums(mat_up))
  original_col_sums <- round(colSums(mat_up))

  rounded_mat_up <- floor(mat_up)
  # Parte decimal de la matriz original
  diff_mat <- mat_up - rounded_mat_up

  # Residuos para preservar suma de filas y columnas en la matriz redondeada
  row_diffs <- original_row_sums - rowSums(rounded_mat_up)
  col_diffs <- original_col_sums - colSums(rounded_mat_up)

  # Si el total de elementos a redondear en filas y columnas es distinto,
  # no hay solución exacta
  posible.sol <- sum(row_diffs) == sum(col_diffs)
  if(!posible.sol)
    warning("There is no exact solution.")

  # Elements of the matrix in decreasing order by decimal part
  idx <- order(diff_mat, decreasing = T)

  if(posible.sol) # if there exist a solution, search for the best
    idx <- bg(idx, row_diffs, col_diffs) # select elements to round up
  else if(idx[1] < 0){ # if there isn't an exact solution or can't find the best
    idx <- order(diff_mat, decreasing = T) # prepare for next loop
  }
  # if idx is a solution, accept it
  # if it is the vector with every element in the matrix, a greedy algorithm
  # is run to choose the best first solution
  for(i in idx){
    # Row and column of the element
    row <- (i-1)%%nrow(Y)+1
    col <- floor((i-1)/nrow(Y))+1
    if(row_diffs[row]>0 && col_diffs[col]>0){
      # round the element up
      rounded_mat_up[row,col] <- rounded_mat_up[row,col] + 1
      row_diffs[row] <- row_diffs[row]-1
      col_diffs[col] <- col_diffs[col]-1
    }
  }

  rounded_mat <- rounded_mat_up / (10**digits)

  return(rounded_mat)
}

bg <- function(idx, row_diffs, col_diffs, elementos_marcados = list()){
  # If everything selected: return solution
  if(sum(c(row_diffs, col_diffs)) == 0)
    return(elementos_marcados)
  # If there are no elements left: out
  if(length(idx) == 0)
    return( -1 )

  # Choose the most critical element to be rounded up.
  e <- idx[1]

  # Row and column of the chosen element
  row <- (e-1)%%length(row_diffs)+1
  col <- floor((e-1)/length(row_diffs))+1

  # If the element can be rounded up, it is marked
  if(row_diffs[row]>0 && col_diffs[col]>0){
    row_diffs[row] <- row_diffs[row]-1
    col_diffs[col] <- col_diffs[col]-1
    elementos_marcados <- append(elementos_marcados, e)
  }

  s <- bg(idx[-1], row_diffs, col_diffs, elementos_marcados)

  # If s=-1, the chosen element must be eliminated because
  # it doesn't lead to a valid solution
  if(s[1] < 0 && e == elementos_marcados[[length(elementos_marcados)]]){
    row_diffs[row] <- row_diffs[row]+1
    col_diffs[col] <- col_diffs[col]+1
    elementos_marcados <- head(elementos_marcados, -1)
    s <- bg(idx[-1], row_diffs, col_diffs, elementos_marcados)
  }

  return(s)
}

#' Round univariate
#'
#' Rounds a vector preserving the rounded sum.
#' @param x A vector.
#' @param digits Number of decimal places to be rounded to.
#' @export
round_vector <- function(x, digits = 0){
  x <- as.vector(x)
  up <- 10**digits
  x <- x*up
  y <- floor(x)
  indices <- tail(order(x-y), round(sum(x)) - sum(y))
  y[indices] <- y[indices] + 1
  return(y/up)
}

#' Round a matrix
#'
#' Returns an integer matrix that preserves the rounded colSums and rowSums.
#' @param Y A matrix.
#' @param digits Decimal places to round to.
#' @param MARGIN One of
#' \itemize{
#'  \item{0} Preserves the rounded colSums and rowSums.
#'  \item{1} Preserves the rounded rowSums independently of each other.
#'  \item{2} Preserves the rounded colSums independently of each other.
#' }
#' @returns The rounded matrix.
#' @details
#' The function will throw a *warning* if the problem is infeasable. To be able
#' to round the matrix in this fashion, the following things must be equal:
#' \itemize{
#'  \item {the sum of the differences between the row totals and
#'  the rounded row totals}
#'  \item {the sum of the differences between the column totals and
#'  the rounded row totals}
#' }
#' @examples
#' set.seed(6)
#' Y <- rnorm(3*5)*10 |> matrix(3,5) |> round(3)
#' X <- round_matrix(Y)
#' Y
#' X
#' colSums(Y) |> round()
#' colSums(X)
#' rowSums(Y) |> round()
#' rowSums(X)
#' @export
round_matrix <- function(Y, digits = 0, MARGIN = 0) {
  if (MARGIN == 0) {
    X <- round_matrix_bivariate(Y,digits)
  } else  if (MARGIN == 1) {
    X <- apply(Y, MARGIN = 1, round_vector, digits = digits) |> t()
  } else  if (MARGIN == 2) {
    X <- apply(Y, MARGIN = 2, round_vector, digits = digits)
  } else {
    stop("MARGIN must be 0, 1 or 2.")
  }
  return(X)
}
