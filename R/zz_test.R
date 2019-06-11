qdecr_read_ocn <- function(vw, stack) load.mgh(vw$stack$ocn.mgh[[stack]])
qdecr_read_ocn_mask <- function(vw, stack) qdecr_read_ocn(vw, stack)$x > 0

qdecr_read_coef <- function(vw, stack) load.mgh(vw$stack$coef[[stack]])
qdecr_read_p <- function(vw, stack) load.mgh(vw$stack$p[[stack]])
qdecr_read_t <- function(vw, stack) load.mgh(vw$stack$t[[stack]])
qdecr_read_se <- function(vw, stack) load.mgh(vw$stack$se[[stack]])

`%AND%` <- function(map1, map2) map1 & map2
`%OR%` <- function(map1, map2) map1 | map2
`%XOR%` <- function(map1, map2) (map1 & !map2) | (!map1 & map2)
`%MASK%` <- function(map1, map2) map1 * map2

AND <- function(map1, map2) map1 %AND% map2
OR <- function(map1, map2) map1 %OR% map2
XOR <- function(map1, map2) map1 %XOR% map2
MASK <- function(map1, map2) map1 %MASK% map2



# demo_plot <- function(vw, stack = NULL, software = c("tksurfer", "freeview"), verbose = TRUE){
#   if (is.null(stack)) stop("`stack` not defined. Please choose: ", paste(stacks(vw), collapse = ", "))
#   software <- match.arg(software)
#   if (software == "freeview"){
# 		#freeview --surface $SUBJECTS_DIR/fsaverage/surf/${hemi}.inflated:overlay=./${OBN2}.cache.th13.abs.sig.cluster.mgh:overlay_threshold=1.3,4 ${cmds}
#     cmdStr <- paste0("freeview --surface ", vw$paths$dir_subj, "/fsaverage/surf/", vw$input$hemi, ".inflated:overlay=", vw$stack$cluster.mgh[[stack]])
#     system(cmdStr, ignore.stdout = !verbose)
#   } else if (software == "tksurfer") {
#     #cmdStr <- paste("tksurfer", vw$input$target, vw$input$hemi, "inflated", "-overlay", vw$stack$coef[[stack]], "-overlay", vw$stack$cluster.mgh[[stack]])
#     cmdStr <- paste("tksurfer", vw$input$target, vw$input$hemi, "inflated", "-overlay", vw$stack$cluster.mgh[[stack]])
#     message2(verbose = verbose, cmdStr)
#     system(cmdStr, ignore.stdout = !verbose)
#   }
# }

