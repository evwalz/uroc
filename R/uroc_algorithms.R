# This function computes the uROC curve as an weighted average over the ROC curves from the associated rocm

Exact <- function(response, predictor) {

  n <- length(response)
  thresholds <- unique(response)[-1]
  N <- length(thresholds)+1

  ncontrols <- (which(duplicated(response) == FALSE) - 1)[-1]
  ncases <- (n - ncontrols)
  weights_all <- (ncases * ncontrols)
  weights = sum(weights_all)

  pre.order <- order(predictor, decreasing=TRUE)
  predictor <- predictor[pre.order]
  response <- response[pre.order]
  dups <- rev(duplicated(rev(predictor)))
  
  response.binary <- lapply(thresholds, function(x, response) {as.numeric(response >= x)}, response)
  Hitrate.all <- lapply(response.binary, function(x, dups) {c(0, (cumsum(x == 1) * sum(x == 0))[!dups])}, dups)
  Farate.all <- lapply(response.binary, function(x, dups) {c(0, (cumsum(x == 0) / sum(x == 0))[!dups])}, dups)
  farate.unique <- sort(unique(c(unlist(Farate.all), 0)))

  #Initialization
  final.hit.min.new <- rep(0, length(farate.unique))
  final.hit.max.new <- rep(0, length(farate.unique))
  hit.min.new <- rep(0, length(farate.unique))
  hit.max.new <- rep(0, length(farate.unique))
  
  for (i in 1:(N-1)) {
    hitrate <- Hitrate.all[[i]]
    farate <- Farate.all[[i]]

    far.group.min <- which(as.numeric(duplicated(farate)) == 0)

    #only select far-values with steps in graph and corresponding min and max value of hitrate
    far.min <- farate[far.group.min]
    hit.min <- hitrate[far.group.min]
    
    roc_curve_max <- approxfun(x = farate, y = hitrate, method = "linear", ties = "ordered")
    hit.max.new <- roc_curve_max(farate.unique)

    hit.min.new <- hit.max.new
    indx.min <- match(far.min, farate.unique)
    hit.min.new[indx.min] <- hit.min
    
    final.hit.min.new <- final.hit.min.new + hit.min.new
    final.hit.max.new <- final.hit.max.new + hit.max.new

  }

  return(list(Farate = sort(c(farate.unique, farate.unique)), Hitrate = sort(c(final.hit.min.new, final.hit.max.new) / weights)))
}


# This function computes an approximation of the uROC curve and takes into account ties in the predictor

Approx1 = function(response, predictor, space.size) {

  n <- length(response)
  N <- length(unique(response))

  ncontrols <- (which(duplicated(response) == FALSE) - 1)[-1]
  ncases <- (n-ncontrols)
  weights_all <- (ncases*ncontrols)
  weights <- sum(weights_all)

  pre.order <- order(predictor, decreasing = TRUE)
  predictor <- predictor[pre.order]
  dups <- rev(duplicated(rev(predictor)))

  if (is.null(space.size)) {
    space.size <- 100
  }

  InterPoint <- seq(0, 1, (1 / space.size))
  Hit.weighted <- rep(0, length(InterPoint))

  for (i in 1:(N-1)) {
    controls <- ncontrols[i]
    cases <- ncases[i]
    response.binary <- c(rep(0, controls), rep(1, cases))[pre.order]
    hitrate <- c(0, cumsum(response.binary == 1)[!dups] * controls)
    farate <- c(0, cumsum(response.binary == 0)[!dups] / controls)

    roc_curve <- approxfun(x = farate, y = hitrate, method = "linear", ties = "ordered")
    Hit.weighted <- roc_curve(InterPoint) + Hit.weighted
  }
  return(list(Farate = c(0, InterPoint), Hitrate = sort(c(0, Hit.weighted / weights))))
}


# This function computes an approximation of the uROC curve by ignoring ties in the predictor

Approx2 <- function(response, predictor, space.size) {

  n = length(response)
  N = length(unique(response))

  ncontrols <- (which(duplicated(response) == FALSE) - 1)[-1]
  ncases <- (n - ncontrols)
  weights_all <- (ncases * ncontrols)
  weights = sum(weights_all)

  pre.order <- order(predictor, decreasing=TRUE)
  predictor <- predictor[pre.order]
  dups <- rev(duplicated(rev(predictor)))

  if(is.null(space.size)){
    space.size <- 500
  }

  InterPoint <- seq(0, 1, (1 / space.size))
  Hit.weighted <- rep(0, length(InterPoint))

  for (i in 1:(N-1)) {
    controls <- ncontrols[i]
    cases <- ncases[i]
    resBin <- c(rep(0, controls), rep(1, cases))[pre.order]
    hitrate <- c(0, cumsum(resBin == 1)[!dups] * controls)
    farate <- c(0, cumsum(resBin == 0)[!dups] / controls)

    roc_curve <- stepfun(x = farate[-1], y = hitrate)  
    Hit.weighted <- Hit.weighted + roc_curve(InterPoint)
  }

  return(list(Farate = c(0, InterPoint), Hitrate = sort(c(0, Hit.weighted / weights))))
}

