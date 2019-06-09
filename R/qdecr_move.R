qdecr_move <- function(results, stack, so) {
  for (i in seq_along(results)){
    for (j in seq_along(stack$names)){
      bsfbm2mgh(results[[i]], stack[[so[i]]][[j]], filter = j)
    }
  }
}