#' DR statistic (R version)
#'
#' Computes the Jetz DR rate. Code modified from fisse
#'
#' @param x an ape::phylo object
#'
#' @return vector of rates for each tip
#' @export
#' @examples
#' DR_statistic(ape::rcoal(10))

DR_statistic <- function(x) {
  rootnode <- length(x$tip.label) + 1
  sprates <- numeric(length(x$tip.label))
  names(sprates) <- x$tip.label
  # Reference these vectors directly; saves an expensive
  # hash table lookup in the hot loop below.
  elv <- x$edge.length
  ev1 <- x$edge[, 1]
  ev2 <- x$edge[, 2]
  # For each tip...
  for (i in 1:length(sprates)) {
      # Our current node of interest is the tip
      node <- i
      index <- 1
      qx <- 0
      # Traverse up to the root
      while (node != rootnode) {
        # Get the length of the edge subtending the focal node
        el <- elv[ev2 == node]
        # Set the new focal node to the current node's parent (traversal step)
        node <- ev1[ev2 == node]
        # Calculate statistic at this node
        qx <- qx + el* (1 / 2^(index - 1))
        # Advance the current index to permit decay of statistic
        index <- index + 1
      }
      sprates[i] <- 1/qx
  }
  return(sprates)
}
