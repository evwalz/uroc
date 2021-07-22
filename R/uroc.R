#' @title Computes a UROC curve.
#' @description This function builds a UROC curve and returns a "uroc" object, a list of class "uroc".
#' @details The default option to compute uroc curve generates an approximation to the UROC curve by using linear interpolation of each ROC curve. For small datasets or binary response the exact uroc curve is computed. Setting option \code{approx = TRUE} uses a faster approximation algorithm where computation time can be further reduced by setting the paramter \code{split} which selects only a subset of ROC curves in the computation.
#' @param response a numeric vector of real valued responses.
#' @param predictor a numeric vector of the same length than \code{response}, containing real valued predictions for each observation.
#' @param approx Boolean. If TRUE approximates true roc curve with faster algorithm.
#' @param split integer value with a default of \code{split = 1}. Computes uroc curve by considering only a subset of all N-1 available ROC curves to reduce computation time. The split parameter defines the distance between a set of equidistant indices which are then used to select particular ROC curves among the N-1.
#'
#' @importFrom graphics text lines
#' @importFrom stats approx
#' @importFrom stats aggregate
#'
#' @return If \code{object = TRUE} this function returns a list of class "uroc".
#' @export
#'
#' @examples
#' data(longley)
#' response = longley$Employed
#' predictor = longley$GNP
#' uroc(response, predictor)


uroc <- function(response, predictor, approx = FALSE, split = 1) {

  if (!is.vector(predictor) || !is.vector(response)) {
    stop("Input must be a vector")
  }

  if (!(is.numeric(predictor) || is.logical(predictor))) {
    stop("predictor must be numeric")
  }

  if (!(is.numeric(response) || is.logical(response))) {
    stop("response must be numeric")
  }

  if (anyNA(response) || anyNA(predictor)){
    stop("missing values in the data")
  }

  response_order <- order(response, decreasing=FALSE)
  response <- response[response_order]
  predictor <- predictor[response_order]
  
  thresholds <- unique(response)
  N <- length(unique(response))
  n <- length(response)

  split <- floor(split)
  if (split > N || split < 0) {
    stop("invalid value for split")
  }

  if (N < 2) {
    stop("response must have more than one level")
  }

  if (n != length(predictor)) {
    stop("response and predictor should have the same length")
  }

  if (N == 2) {
    uroc_object <- uroc_exact(response, predictor, thresholds, n, N)
    return(structure(list(farate = uroc_object$farate, hitrate = uroc_object$hitrate), class = "uroc"))
  }

  if (approx) {
    uroc_object <- uroc_approx(response, predictor, n, N, split)
    return(structure(list(farate = uroc_object$farate, hitrate = uroc_object$hitrate), class = "uroc"))
  }

  thresh_boolean <- !duplicated(response)
  thresholds_index <- which(thresh_boolean) - 1
  ncontrols <- thresholds_index[-1]
  weights <- ncontrols*(n-ncontrols)

  weight_s <- sum(weights)

  class_predictor <- cumsum(duplicated(sort(predictor))==FALSE)[rank(predictor, ties.method = "first")]

  splitted_classes = aggregate(data.frame(classes = class_predictor), by = data.frame(thresholds = cumsum(thresh_boolean)), FUN = identity)$classes

  # compute first roc
  cat("Estimating uroc...\n")
  pb <- utils::txtProgressBar(style = 1)

  response_binary <- rep(1, n)
  response_binary[1:ncontrols[1]] <- 0

  predictor_order <- order(predictor, decreasing=TRUE)
  response_binary <- response_binary[predictor_order]
  predictor_sorted <- rev(predictor[predictor_order])

  predictor_unique <- unique(predictor_sorted)
  predictor_unique_indx <- which(!duplicated(predictor_sorted))

  dups <- n +1 - rev(predictor_unique_indx)

  tpr <- c(0, cumsum(response_binary)[dups])
  fpr <- c(0, cumsum(response_binary == 0)[dups])
  tpr_weight <- rev(tpr)
  fpr_weight <- rev(fpr)
  interpoint <- seq(0, 1, (1 / 1000))


  tpr_interpolated <- approx(x = fpr/(ncontrols[1]), y = tpr, xout = interpoint, method = "linear", ties = "ordered")$y * ncontrols[1]
  sum_tpr_fpr <- tpr_weight + fpr_weight

  for (i in 2:(N-1)) {
    utils::setTxtProgressBar(pb, i/(N-2))
    diff_split_element <- diff(sort(c(0, splitted_classes[[i]])))
    m = length(diff_split_element)
    sum_indicator <- rep(c(m:1), diff_split_element)
    seq_change <- length(sum_indicator)
    tpr_weight[1:seq_change] <- tpr_weight[1:seq_change] - sum_indicator
    fpr <- (sum_tpr_fpr - tpr_weight) / ncontrols[i]
    tpr_interpolated <- approx(x = rev(fpr), y = rev(tpr_weight), xout = interpoint, method = "linear", ties = "ordered")$y * ncontrols[i] + tpr_interpolated
  }
  close(pb)

  tpr_interpolated_weight <- tpr_interpolated / weight_s

  return(structure(list(farate = c(0, interpoint), hitrate = c(0, tpr_interpolated_weight)), class = "uroc"))
}


#' @title Plot of UROC curve.
#' @description This function plots a UROC curve.
#' @method plot uroc
#' @param x object of class "uroc"
#' @param ... further arguments to \code{\link{plot}}
#' @importFrom graphics plot
#' @return plot
#' @export


plot.uroc <- function(x, ...) {

  farate <- x$farate
  hitrate <- x$hitrate
  cpa_approx <- trap(farate, hitrate)
  cpa_value <- round(cpa_approx, 2)
  plot(x = farate, y = hitrate, type = "l", xlab = "1-Specificity", ylab = "Sensitivity",...)
  lines(x = c(0, 1), y = c(0, 1), lty = 2)
  text(x = 0.6, y = 0.4, labels = paste("CPA:",cpa_value))

}

uroc_exact <- function(response, predictor, thresholds, n, N) {
  thresholds <- thresholds[-1]
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
  cat("Estimating uroc...\n")
  pb <- utils::txtProgressBar(style = 1)
  for (i in 1:(N-1)) {
    utils::setTxtProgressBar(pb, i/(N-1))

    hitrate <- Hitrate.all[[i]]
    farate <- Farate.all[[i]]

    far.group.min <- which(as.numeric(duplicated(farate)) == 0)

    #only select far-values with steps in graph and corresponding min and max value of hitrate
    far.min <- farate[far.group.min]
    hit.min <- hitrate[far.group.min]

    roc_curve_max <- approx(x = farate, y = hitrate, xout = farate.unique, method = "linear", ties = "ordered")
    hit.max.new <- roc_curve_max$y

    hit.min.new <- hit.max.new
    indx.min <- match(far.min, farate.unique)
    hit.min.new[indx.min] <- hit.min

    final.hit.min.new <- final.hit.min.new + hit.min.new
    final.hit.max.new <- final.hit.max.new + hit.max.new

  }
  close(pb)
  return(list(farate = sort(c(farate.unique, farate.unique)), hitrate = sort(c(final.hit.min.new, final.hit.max.new) / weights)))

}


uroc_approx <- function(response, predictor, n, N, split) {
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


  # compute first roc curve
  cat("Estimating uroc...\n")
  pb <- utils::txtProgressBar(style = 1)

  order_predictor <- order(predictor, decreasing = TRUE)
  first_threshold <- min(response)
  response_binary <- response[order_predictor] > first_threshold

  dups <- rev(duplicated(rev(predictor[order_predictor])))
  tp <- c(0, cumsum(response_binary == 1)[!dups])
  fp <- c(0, cumsum(response_binary == 0)[!dups])
  truepositive <- rev(tp)
  falsepositive <- rev(fp)

  space_size <- 1000
  sum_tp_fp <- falsepositive + truepositive
  InterPoint <- seq(0, 1, (1 / space_size))

  first_roc_curve <- approx(x = fp/ncontrols[1], y = tp, xout = InterPoint, method = "linear", ties = "ordered")
  Hit_weighted <- first_roc_curve$y * ncontrols_split[1]
  lit <- length(indx_transformation)
  for (i in 2:lit) {
    utils::setTxtProgressBar(pb, i/(lit-1))
    m <- length(Split_classes_predictor_ordered_diff[[i]])
    sum_indicator <- rep(seq(m,1,-1), Split_classes_predictor_ordered_diff[[i]])
    sequence_to_change <- length(sum_indicator)
    truepositive[1:sequence_to_change] <- truepositive[1:sequence_to_change] - sum_indicator
    farate <- (sum_tp_fp - truepositive) / ncontrols_split[i]
    roc_approx <- approx(x = rev(farate), y = rev(truepositive) * ncontrols_split[i], xout = InterPoint, method = "linear", ties = "ordered")
    Hit_weighted <- roc_approx$y + Hit_weighted

  }
  close(pb)
  return(list(farate = c(0, InterPoint), hitrate = c(0, Hit_weighted / weights)))
}
