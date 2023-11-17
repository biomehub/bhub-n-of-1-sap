compute_iauc <- function(t, y) {
  stopifnot("Must have at least 3 time points" = length(t) > 2)
  baseline <- y[1]
  yc <- y - baseline
  seg.type <- auc <- numeric(length = length(t) - 1)
  auc[1] <- yc[2] * (t[2]-t[1])/2 * as.numeric(yc[2] > 0)
  for (i in 3:length(t)) {
    delta_t <- t[i] - t[i-1]
    current_y <- y[i]
    current_yc <- yc[i]
    previous_y <- y[i-1]
    previous_yc <- yc[i-1]
    if (current_yc >= 0 & previous_yc >= 0) { 
      auc[i-1] <- delta_t * (current_yc + previous_yc)/2 
      seg.type[i-1] <- 1
    } else if (current_yc >= 0 & previous_yc < 0) { 
      auc[i-1] <- delta_t * (current_yc^2 / (current_y-previous_y))/2 
      seg.type[i-1] <- 2
    } else if (current_yc < 0 & previous_yc >= 0) {  
      auc[i-1] <- delta_t * (previous_yc^2/(previous_y-current_y))/2
      seg.type[i-1] <- 3
    } else if (current_yc < 0 & previous_yc < 0) { 
      auc[i-1] <- 0             
      seg.type[i-1] <- 4
    } else {
      # The above cases are exhaustive, so this should never happen
      stop(paste0("i:", i, "Error: No condition met\n"))
    }
  }
  return(list(auc=sum(auc), segments=auc, seg.type=seg.type))
}
