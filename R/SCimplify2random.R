require(SuperCell)
require(igraph)

# As every supercell should contain at least one cell, this corrects for zero-
# cell supercell (computationally inefficient)
sample.correction <- function(m, samp, gamma){
    samples <- names(samp)
    for(s in samples){
        idx <- which(m %in% samp[[s]])
        toRedraw <- setdiff(samp[[s]], m[idx])
        while(length(toRedraw) != 0){
            for(el in toRedraw){
                idx.sampled <- sample(idx, round(gamma/2))
                m[idx.sampled] <- el
            }
            toRedraw <- setdiff(samp[[s]], m[idx])
        }
    }
    return(m)
}


#' from Scimplify result to random simplification, modified in order to randomize
#' memberships per sample
#'
#' @export

SCimple2Random <- function(SC, gamma, seed = 12345){
  if(is.null(SC$graph.singlecell)){
    stop("SCimply2random is only available for SCiplify result whith return.singlecell.NW == TRUE")
  }

  N.c          <- length(SC$membership)
  N.SC         <- round(N.c/gamma)

  set.seed(seed)
  if(!is.null(SC$sc.cell.samples)){
      samples = unique(SC$sc.cell.samples)
      SC.per.samples <- lapply(table(SC$sc.cell.samples), function(x) round(x / gamma))
      membership.per.samples <- list()
      curCount = 0
      for(sample in samples){
        membership.per.samples[[sample]] <- (1 + curCount):(SC.per.samples[[sample]] + curCount)
        curCount <- curCount + SC.per.samples[[sample]]
      }
      r.membership <- sapply(SC$sc.cell.samples, function(x) sample(membership.per.samples[[x]], 1))
      r.membership <- sample.correction(r.membership, membership.per.samples, gamma)
  }else{
    r.membership <- c(1:N.SC, sample(N.SC, N.c-N.SC, replace = T))
  }
  SC.random                      <- SC
  SC.random$membership           <- r.membership
  SC.random$supercell_size       <- as.vector(table(r.membership))
  SC.random$simplification.algo  <- "Random_grouping"

  SC.NW                          <- igraph::contract(SC.random$graph.singlecell, r.membership)
  SC.NW                          <- igraph::simplify(SC.NW, remove.loops = T, edge.attr.comb="sum")

  igraph::E(SC.NW)$width         <- sqrt(igraph::E(SC.NW)$weight/10)
  igraph::V(SC.NW)$size          <- SC.random$supercell_size
  igraph::V(SC.NW)$sizesqrt      <- sqrt(igraph::V(SC.NW)$size)

  SC.random$graph.supercells     <- SC.NW
  SC.random$seed                 <- seed

  return(SC.random)
}
