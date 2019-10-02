# compute area under curve using trapezoidal rule


Trap <- function(farate, hitrate) {

  diffs.far <- farate[-1] - farate[-length(farate)]
  means.vert <- (hitrate[-1] + hitrate[-length(hitrate)]) * 0.5

  return(sum(means.vert * diffs.far))
}



# create names for png file based on iteraton index

rename <- function(i) {

  if (i <= 10) {
    i <- i - 1
    return(name <- paste('ani', '000', i, '.png', sep = ''))
  }

  if (i <= 100 && i > 10) {
    i <- (i - 1)
    return(name <- paste('ani', '00', i, '.png', sep = ''))
  }

  if (i <= 1000 && i > 100) {
    i <- (i - 1)
    return(name <- paste('ani', '0', i, '.png', sep = ''))
  }

  if (i > 1000) {
    return(name <- paste('ani', i, 'png', sep = ''))
  }
}
