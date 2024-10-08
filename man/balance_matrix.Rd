% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/balance_matrix.R
\name{balance_matrix}
\alias{balance_matrix}
\title{Balance matrix}
\usage{
balance_matrix(Y, col_totals = NULL, row_totals = NULL, allow_negative = TRUE)
}
\arguments{
\item{Y}{Matrix to be balanced.}

\item{col_totals}{(optional) Desired sum of columns.}

\item{row_totals}{(optional) Desired sum of rows.}

\item{allow_negative}{Are negative entries in the balanced matrix allowed?}
}
\value{
The balanced matrix.
}
\description{
Balances a matrix so that the columns and/or rows add up
to a certain vector.
}
\details{
Balancing is done according to the criteria of minimum sum
of squares.

If neither \code{col_totals} nor \code{row_totals} is given, the same matrix will be
returned. If only one of them is given, only that axis will be
balanced.
}
\examples{
set.seed(2)
Y <- rnorm(3*5) |> matrix(3,5) |> round(3)
v <- c( 0.876, -1.078, 3.452, 0.261, 1.349)
h <- c(-1.851, 0.243, 6.468)
X1 <- balance_matrix(Y, v, h)
Y
X1
h
rowSums(X1)
v
colSums(X1)
X3 <- balance_matrix(Y, col_totals = v)
v
colSums(X3)
X4 <- balance_matrix(Y, row_totals = h)
h
rowSums(X4)
}
