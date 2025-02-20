#' assign_control
#' to be documented
#' @usage assign_control(control=list(),return=FALSE, env = parent.frame())
#' @param control to be documented
#' @param env env = parent.frame()
#' @noRd
#' @return to be documented
reassign_control <- function(control=list(),return=FALSE, env = parent.frame()) {
    con=list()
    nmsC <- names(control)
    con[(namc <- names(control)[names(control)!=""])] <- control[names(control)!=""]
    if (length(noNms <- namc[!namc %in% nmsC]))  warning("unknown names in control: ", paste(noNms, collapse = ", "))
    if(!return) for(i in names(con))
    {
      assign(i,con[[i]],envir =parent.frame())
    } else {
      for(i in names(con))
    {
      control[[i]]<-get(i,envir =parent.frame())
      }
    control
    }
  }
