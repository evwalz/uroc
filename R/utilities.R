# compute area under curve using trapezoidal rule


Trapezoidal <- function(farate, hitrate) {

  diffs.far <- farate[-1] - farate[-length(farate)]
  means.vert <- (hitrate[-1] + hitrate[-length(hitrate)]) * 0.5

  return(sum(means.vert * diffs.far))
}



# create names for png file based on iteraton index

rename_figure <- function(indx) {
    return(name <- paste('animation_', indx, '.png', sep = ''))
}
