cell.quadtree.leaf <-
function(q, xylim) {
	  data.frame(id = q$id, 
	             x = c(xylim[1,1], xylim[2,1], xylim[2,1], xylim[1,1], xylim[1,1]),
	             y = c(xylim[1,2], xylim[1,2], xylim[2,2], xylim[2,2], xylim[1,2]))
	}
