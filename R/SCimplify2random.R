require(SuperCell)
require(igraph)


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
      
      # Attribute one cell to each available supercell first to ensure each
      # supercell has at least one element
      chosen.cells <- sample(sc.id, length(super.id), replace = FALSE)
      sc.id <- setdiff(sc.id, chosen.cells)
      new_members <- super.id
      names(new_members) <- sc.names[chosen.cells]
      r.membership <- c(r.membership, new_members)
      
      # Pick among all supercells with replacement to assign single-cell
      new_members <- sample(super.id, length(sc.id), replace = TRUE)
      names(new_members) <- sc.names[sc.id]
      r.membership <- c(r.membership, new_members)
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
  
  SC.random$SC.cell.annotation. <- supercell_assign(SC$sc.cell.annotation., r.membership)
  return(SC.random)
}
