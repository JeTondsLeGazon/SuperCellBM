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


# Randomize supercell memberships per groups (label or sample) in order to keep
# grouped cell together and not randomize totally different cells, which could
# yield really bad results but would not be meaningful
randomizeGroups <- function(membership, groups){
    
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
  
  groups <- unique(SC$SC.cell.annotation.)
  r.membership <- c()
  sc.names <- names(SC$sc.cell.annotation.)
  for(group in groups){  # groups may be label / samples / etc.
      super.id <- which(SC$SC.cell.annotation. == group)
      sc.id <- which(SC$sc.cell.annotation. == group)
      for(supercell.group in super.id){
          N <- SC$supercell_size[supercell.group]
          chosen.cells <- sample(sc.id, N)
          sc.id <- setdiff(sc.id, chosen.cells)
          new_member <- rep(supercell.group, N)
          names(new_member) <- sc.names[chosen.cells]
          r.membership <- c(r.membership, new_member)
      }
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
