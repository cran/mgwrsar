cell.quadtree <-
function(q, xylim, ...) {
	  i <- q$index
	  j <- 3 - q$index
	  clip <- function(xylim.clip, i, upper) {
	    if (upper) xylim.clip[1, i] <- max(q$threshold, xylim.clip[1,i]) else 
	      xylim.clip[2,i] <- min(q$threshold, xylim.clip[2,i])
	    xylim.clip
	  } 
	  d <- data.frame(id=NULL, x=NULL, y=NULL)
	  if(q$threshold > xylim[1,i]) d <- cell(q$lower, clip(xylim, i, FALSE), ...)
	  if(q$threshold < xylim[2,i]) d <- rbind(d, cell(q$upper, clip(xylim, i, TRUE), ...))
	  d
	}
