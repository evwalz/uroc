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

Approx <- function(response, predictor, n, N, split, space_size) {
  
  # compute controls and weights for the N/split transformed problems 
  ncontrols <- (which(duplicated(response) == FALSE) - 1)[-1]
  ncontrols_split <- ncontrols[seq(1, length(ncontrols), split)]
  group_length <- c(ncontrols_split) - c(0, ncontrols_split[-length(ncontrols_split)])
  groups <- rep(seq(1, length(group_length),1), group_length)
  
  weights <- sum(ncontrols_split * (n - ncontrols_split))
  
  # indices to define which transformed problems are used in the computation of uroc
  indx_transformation <- seq(1, (N-1), split)
  
  Ranking_predictor_1 <- rank(predictor, ties.method = "first")
  Classes_predictor <- cumsum(duplicated(sort(predictor))==FALSE)[Ranking_predictor_1]
  
  Split_classes_predictor <- split(Classes_predictor[1:length(groups)], groups)
  Split_classes_predictor_ordered <- lapply(Split_classes_predictor, function(x){c(0, sort(x))})
  rm(Split_classes_predictor)
  Split_classes_predictor_ordered_diff <- lapply(Split_classes_predictor_ordered, function(x){x[-1] - x[-length(x)]})
  rm(Split_classes_predictor_ordered)
  
  
  # compute firs roc curve  
  order_predictor <- order(predictor, decreasing = TRUE)
  first_threshold <- min(response)
  response_binary <- response[order_predictor] > first_threshold

  dups <- rev(duplicated(rev(predictor[order_predictor])))
  tp <- c(0, cumsum(response_binary == 1)[!dups])
  fp <- c(0, cumsum(response_binary == 0)[!dups])
  truepositive <- rev(tp)
  falsepositive <- rev(fp)
  first_roc_curve <- approxfun(x = fp/ncontrols[1], y = tp, method = "linear", ties = "ordered")
  
  sum_tp_fp <- falsepositive + truepositive   
  InterPoint <- seq(0, 1, (1 / space_size))
  Hit_weighted <- first_roc_curve(InterPoint) * ncontrols_split[1]

  for (i in 1:length(indx_transformation)) {
    m <- length(Split_classes_predictor_ordered_diff[[i]])
    
    sum_indicator <- rep(seq(m,1,-1), Split_classes_predictor_ordered_diff[[i]])
    
    sequence_to_change <- length(sum_indicator)
    
    truepositive[1:sequence_to_change] <- truepositive[1:sequence_to_change] - sum_indicator
    
    farate <- (sum_tp_fp - truepositive) / ncontrols_split[i]
    
    roc_curve <- approxfun(x = rev(farate), y = rev(truepositive) * ncontrols_split[i], method = "linear", ties = "ordered")
    
    Hit_weighted <- sort(roc_curve(InterPoint)) + Hit.weighted

  }
  return(list(Farate = c(0, InterPoint), Hitrate = sort(c(0, Hit_weighted / weights))))
}


