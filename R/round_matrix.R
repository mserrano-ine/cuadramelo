#' Round a matrix
#'
#' Returns an integer matrix that preserves the rounded colSums and rowSums.
#' @param Y A matrix.
#' @param digits Decimal places to round to.
#' @returns The rounded matrix.
#' @details
#' The function will throw a *warning* if the problem is infeasable. To be able
#' to round the matrix in this fashion, the following things must be equal:
#' \itemize{
#'  \item{a}{the sum of the differences between the row totals and
#'  the rounded row totals}
#'  \item{b}{the sum of the differences between the column totals and
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
#' set.seed(2)
#' Y <- rnorm(3*5)*10 |> matrix(3,5) |> round(3)
#' X <- round_matrix(Y)
#' Y
#' X
#' colSums(Y) |> round()
#' colSums(X)
#' rowSums(Y) |> round()
#' rowSums(X)
#' @export
round_matrix <- function(Y, digits=0) {
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
    warning("No hay solución exacta")

  # Elementos de la matriz ordenados decrecientemente por parte decimal
  idx <- order(diff_mat, decreasing = T)

  if(posible.sol) # Si hay posible solución, se busca la mejor
    idx <- bg(idx, row_diffs, col_diffs) # Elementos seleccionados para redondear al alza
  else if(idx[1] < 0){ # En caso de no existir solución exacta o no encontrar la mejor
    idx <- order(diff_mat, decreasing = T) # Preparo los elementos para el siguiente bucle
  }
  # Si idx es una solución, simplemente se materializa
  # Si es el vector de todos los elementos de la matriz, se ejecuta un algoritmo
  # voraz para seleccionar la mejor primera solución
  for(i in idx){
    # Fila y columna del elemento de la matriz
    row <- (i-1)%%nrow(Y)+1
    col <- floor((i-1)/nrow(Y))+1
    if(row_diffs[row]>0 && col_diffs[col]>0){
      # Redondeo al alza del elemento
      rounded_mat_up[row,col] <- rounded_mat_up[row,col] + 1
      row_diffs[row] <- row_diffs[row]-1
      col_diffs[col] <- col_diffs[col]-1
    }
  }

  rounded_mat <- rounded_mat_up / (10**digits)

  return(rounded_mat)
}

bg <- function(idx, row_diffs, col_diffs, elementos_marcados = list()){
  # Si se ha seleccionado todo: devolver la solución
  if(sum(c(row_diffs, col_diffs)) == 0)
    return(elementos_marcados)
  # Si no  quedan elementos por seleccionar: salir
  if(length(idx) == 0)
    return( -1 )

  # Selecciona el elemento más prioritario a redondear al alza
  e <- idx[1]

  # Fila y columna del elemento seleccionado
  row <- (e-1)%%length(row_diffs)+1
  col <- floor((e-1)/length(row_diffs))+1

  # Si el elemento se puede redondear al alza, se marca
  if(row_diffs[row]>0 && col_diffs[col]>0){
    row_diffs[row] <- row_diffs[row]-1
    col_diffs[col] <- col_diffs[col]-1
    elementos_marcados <- append(elementos_marcados, e)
  }

  s <- bg(idx[-1], row_diffs, col_diffs, elementos_marcados)

  # Si s=-1 entonces hay que eliminar el elemento seleccionado porque no
  # conduce a una solución válida
  if(s[1] < 0 && e == elementos_marcados[[length(elementos_marcados)]]){
    row_diffs[row] <- row_diffs[row]+1
    col_diffs[col] <- col_diffs[col]+1
    elementos_marcados <- head(elementos_marcados, -1)
    s <- bg(idx[-1], row_diffs, col_diffs, elementos_marcados)
  }

  return(s)
}
