qdecr_make_stack <- function(vw, stats, mcz_thr, clusterfiles = c("cluster.mgh", "cluster.summary", "ocn.annot", "ocn.mgh", "masked.mgh", "voxel.mgh")){

    stack <- list()
    stack$names <- attr(vw$model$mm[[1]], "dimnames")[[2]]
    nn1 <- c(stats, clusterfiles)
    stats2 <- paste0(stats, ".mgh")
    clusterfiles2 <- paste0("cache.th", mcz_thr, ".abs.sig.", clusterfiles)
    nn2 <- c(stats2, clusterfiles2)
    for (n in nn1) stack[[n]] <- list()
    for(i in seq_along(stack$names)){
      for (j in seq_along(nn1)){
        stack[[nn1[j]]][[stack$names[i]]] <- paste0(vw$paths$final_path, ".stack", i, ".", nn2[j])
      }
    }
    stack
}