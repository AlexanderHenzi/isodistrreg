#' Componentwise partial order relation
#'
#' @description
#' Compares all rows of a numeric matrix or data frame \code{X}
#' with respect to the componentwise order. A row \code{X[i, ]} is smaller or
#' equal to a row \code{X[j, ]} in the componentwise order, if \code{X[i, k] <=
#' X[j, k]} for \code{k = 1, ..., ncol(X)}.
#'
#' @usage compOrd(X)
#'
#' @param X a numeric matrix or a data frame containing numeric or ordered
#'   factor variables with at least two columns.
#'   
#' @details
#' The columns of \code{X} are sorted sequentially: First all constraints
#' based on only the first column \code{X[, 1]} are activated (set \code{TRUE}),
#' then constraints are dropped based on the orders of the remaining columns.
#' This avoids \code{nrow(X)^2 / 2} pairwise comparisons.
#'
#' @return A list containing
#'
#' \item{\code{paths}}{a two-column matrix giving all pairs of indices
#'   \code{(i,j)} which satisfy \code{all(X[i, ] <= X[j, ])}.}
#'
#' \item{\code{colOrder}}{a matrix of the columnwise orders of \code{X}. Used to
#'   compute paths and required for other function calls.}
#'
#' @keywords internal
compOrd <- function(X) {
 d <- ncol(X)
 m <- nrow(X)
 smaller <- integer(m * (m - 1) / 2)
 nSmaller <- integer(m)
 lower <- 1L

 colOrder <- matrix(nrow = m, ncol = d) 
 ranks <- matrix(nrow = m, ncol = d) 
 for (j in seq_len(d)) {
  colOrder[, j] <- order(X[, j])
  ranks[, j] <- rank(X[, j], ties.method = "max")
 }
 
 for (k in seq_len(m)) {
  nonZeros <- logical(m)
  # Find all indices k such that X[i, 1] >= X[k, 1]
  nonZeros[colOrder[seq_len(ranks[k, 1]), 1]] <- TRUE
  for (l in 2:d) {
   # Remove TRUE for all indices k with X[i, l] < X[k, l]
   if (ranks[k, l] < m) nonZeros[colOrder[(ranks[k, l] + 1):m, l]] <- FALSE
  }
  nonZeros <- which(nonZeros)
  nNonZero <- length(nonZeros)
  upper <- lower + nNonZero - 1L
  ind <- lower:upper
  nSmaller[k] <- nNonZero
  smaller[ind] <- nonZeros
  lower <- upper + 1L
 }
 paths <- c(smaller[seq_len(upper)], rep.int(seq_len(m), times = nSmaller))
 dim(paths) <- c(length(paths) / 2, 2)
 dimnames(paths) <- list(NULL, c("smaller", "greater"))
 list(paths = paths, colOrder = colOrder)
}

#' Transitive Reduction of path matrix
#'
#' @description
#' Computes the transitive reduction of the path matrix of a directed acyclic
#' graph.
#'
#' @usage
#' trReduc(paths, n)
#'
#' @param paths a two column integer matrix containing all pairs of nodes
#'   which are linked by a path of edges.
#' @param n the number of nodes.
#'
#' @details
#' 
#' For each \code{k}, all indirect paths going through node \code{k} are
#' removed (if there are any).
#'
#' @return
#' A two column matrix giving all pairs of edges \code{(i,j)} such that there
#' exists a directed edge from node \code{i} to node \code{j}.
#' 
#' @keywords internal
trReduc <- function(paths, n) {
  # Construct path matrix from node pairs
  edges <- matrix(nrow = n, ncol = n, FALSE)
  edges[paths] <- TRUE
  diag(edges) <- FALSE
  for (k in 1:n) {
    # Remove all indirect paths which go through node k
    edges[edges[, k], edges[k, ]] <- FALSE
  }
  edges <- which(edges, arr.ind = TRUE)
  colnames(edges) <- c("from", "to")
  edges
}

#' Neighbor points with respect to componentwise order
#'
#' @description Find the neighbor points of the rows of a matrix \code{x} within
#' the rows of \code{X} in the componentwise partial order. That is, for each
#' row \code{x[i, ]}, find all indices \code{k} such that either \code{x[i, ] >=
#' X[k, ]} in all components and \code{x[i, ] >= X[j, ] >= X[k, ]} holds for no
#' \code{j} different from \code{k}, or \code{x[i, ] <= X[k, ]} in all
#' components and \code{x[i ,] <= X[j, ] <= X[k, ]} holds for no \code{j}
#' different from \code{k}.
#'
#' @usage
#' neighborPoints(x, X, orderX)
#'
#' @param x numeric matrix with at least two columns.
#' @param X numeric matrix with same number of columns as \code{x}.
#' @param orderX output of \code{compOrd(X)}.
#'
#' @return
#' Lists of length \code{nrow(x)} giving for each \code{x[i, ]} the indices
#' of the smaller and the greater neighbor points within the rows of \code{X}.
#' 
#' @keywords internal
neighborPoints <- function(x, X, orderX) {
  colOrder <- orderX$colOrder
  nx <- nrow(x)
  k <- ncol(x)
  n <- nrow(X)
  ranksGreater <- ranksSmaller <- matrix(nrow = nx, ncol = k)
  
  for (j in 1:k) {
    # Find positions of the x[i, j] in the sorted vectors X[, j] for all j
    ranksGreater[, j] <- findInterval(x = x[, j], vec = X[colOrder[, j], j])
    ranksSmaller[, j] <- findInterval(x = x[, j], vec = X[colOrder[, j], j],
                                      left.open = TRUE)
  }
  
  xGeqX <- matrix(nrow = n, ncol = nx, FALSE)
  xLeqX <- matrix(nrow = n, ncol = nx, TRUE)
  
  for (i in 1:nx) {
    # Find all indices k such that X[i, 1] >= X[k, 1]
    # (or X[i, 1] <= X[k, 1])
    if (ranksGreater[i, 1] > 0) 
      xGeqX[colOrder[seq_len(ranksGreater[i, 1]), 1], i] <- TRUE
    if (ranksSmaller[i, 1] > 0) 
      xLeqX[colOrder[seq_len(ranksSmaller[i, 1]), 1], i] <- FALSE
    for (j in 2:k) {
      # Remove TRUE for all indices k with X[i, j] < X[k, j]
      # (or X[i, 1] > X[k, 1])
      if (ranksGreater[i, j] < n) 
        xGeqX[colOrder[(ranksGreater[i, j] + 1):n, j], i] <- FALSE
      if (ranksSmaller[i, j] > 0) 
        xLeqX[colOrder[seq_len(ranksSmaller[i, j]), j], i] <- FALSE
    }
  }
  
  # Transitive reduction to filter out adjacent points
  paths <- matrix(nrow = n, ncol = n, FALSE)
  paths[orderX$paths] <- TRUE
  diag(paths) <- FALSE
  for (k in which(apply(xLeqX | xGeqX, 1, any))) {
    xLeqX[paths[k, ], xLeqX[k, ]] <- FALSE
    xGeqX[paths[, k], xGeqX[k, ]] <- FALSE
  }
  
  smaller <- lapply(asplit(xGeqX, 2), which) # lapply(split(t(xGeqX), seq_len(nx)), which)
  greater <- lapply(asplit(xLeqX, 2), which) # lapply(split(t(xLeqX), seq_len(nx)), which)
  list(smaller = unname(smaller), greater = unname(greater))
}