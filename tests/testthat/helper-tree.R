# Flatten the tab tree into a list of leaves with their tab path
tree_leaves <- function(nodes, path = character()) {
  out <- list()
  for (nd in nodes) {
    p <- c(path, nd$title)
    if (.xp_is_leaf(nd)) {
      out[[length(out) + 1L]] <- list(path = p, fn = nd$content)
    } else {
      out <- c(out, tree_leaves(nd$children, p))
    }
  }
  out
}

tree_titles <- function(nodes) vapply(nodes, function(n) n$title, character(1))
