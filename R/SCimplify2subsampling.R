require(SuperCell)
require(igraph)

#' rom Scimplify result to subsampling
#'
#' @export

SCimple2Subsampling <- function(X, SC, gamma, seed = 12345){

  N.c        <- ncol(X)
  N.SC       <- round(N.c/gamma)
  n.pc       <- SC$n.pc
  k.knn      <- SC$k.knn
  N.keep.genes <- length(SC$genes.use)

  set.seed(seed)
  samples <- unique(SC$sc.cell.samples)
  keep.cells.idx <- c()
  for(samp in samples){
      drawable <- which(SC$sc.cell.samples == samp)
      N.sample <- round(length(drawable) / gamma)
      idx <- sample(drawable, N.sample)
      keep.cells.idx <- c(keep.cells.idx, idx)
  }
  keep.cells.idx <- sort(keep.cells.idx)
  keep.cells.ids <- colnames(X)[keep.cells.idx]

  X <- X[, keep.cells.ids]

  keep.genes.var <- sort(apply(X, 1, var), decreasing = T)[1:N.keep.genes]
  keep.genes     <- names(keep.genes.var[keep.genes.var>0])

  X          <- X[keep.genes, ]



  X.for.pca               <- scale(Matrix::t(X))
  if(N.SC > 100){
   PCA          <- irlba::irlba(X.for.pca, nv = max(n.pc))
   PCA$x        <- PCA$u %*% diag(PCA$d)
  } else {
   PCA          <- prcomp(X.for.pca, rank. = max(n.pc), scale. = F, center = F)
  }

  sc.nw <- build_knn_graph(X = PCA$x[,n.pc], k = k.knn, from = "coordinates", use.nn2 = TRUE, dist_method = "euclidean")

  res <- list(graph.singlecell = sc.nw,
              gamma = gamma,
              n.pc = n.pc,
              k.knn = k.knn,
              genes.use = keep.genes,
              cells.use.ids = keep.cells.ids,
              cells.use.idx = keep.cells.idx,
              supercell_size = rep(1, length(keep.cells.idx)),
              simplification.algo  = "Subsampling")
  return(res)

}
