  #' @title Builds the ROC movie (ROCM) an animated sequence of ROC curves.
  #' @description This function computes the sequence of ROC curves which form the ROC Movie and produces a GIF animated ROCM
  #' @details The ROC movie can be used to visualize the performance of a real valued foreacsting problem. Therefore, a sequence of ROC curves is generated which can than be combined into a GIF animation. Each entry of the list consist of two vectors of length 1000 containing the values of farate (1-Specificity) and hitrate (sensitivity) and three values, namely the associated auc value, the weight and the threshold
  #' @param response a numeric vector of real valued responses
  #' @param predictor a numeric vector of the same length than \code{response}, containing real valued predictions for each observation
  #' @param a selects a subset of all ROC curves for the ROC movie with at least \code{a} and at most \code{a+b} ROC curves (default \code{a=400})
  #' @param b selects a subset of all ROC curves for the ROC movie with at least \code{a} and at most \code{a+b} ROC curves (default \code{b=100})
  #' @param object if TRUE a list of ROC curves is returned (default \code{object = TRUE})
  #' @param gif if TRUE a gif animation is created
  #' @param movie.name name of the movie (with extension)
  #' @param ... parameters to control the behavior of the GIF animation using the external function ani.option from \link{animation}.
  #'
  #' @importFrom stats approx
  #' @importFrom animation ani.options saveGIF
  #' @importFrom graphics plot
  #'
  #' @return if \code{object = TRUE}, this function returns a list of ROC curves
  #' @export
  #'
  #' @examples
  #' \dontrun{
  #' data(longley)
  #' response = longley$Employed
  #' predictor = longley$GNP
  #' rocm(response, predictor)}

  rocm <-  function(response,
                    predictor,
                    a = NULL,
                    b = NULL,
                    object = TRUE,
                    gif = FALSE,
                    movie.name = "animation.gif",
                    ...) {

    if (!is.vector(predictor) || !is.vector(response)) {
      stop("Input must be a vector")
    }

    if (!(is.numeric(predictor) || is.logical(predictor))) {
      stop("predictor must be numeric")
    }

    if (!(is.numeric(response) || is.logical(response))) {
      stop("response must be numeric")
    }

    if (anyNA(response) || anyNA(predictor)) {
      stop("missing values in the data")
    }

    if (!is.logical(object)) {
      stop("invalid input for object")
    }

    if (!is.logical(gif)) {
      stop("invalid input for object")
    }

    n <- length(response)

    if (n != length(predictor)) {
      stop("response and predictor should have the same length")
    }

    response_order <- order(response, decreasing=FALSE)
    response <- response[response_order]
    predictor <- predictor[response_order]


    Encoding = rle(response)
    thresholds <- Encoding$values[-1]
    N <- length(thresholds)

    if (N == 0) {
      stop("response must have more than one level")
    }

    if (is.null(a) && N <= 400) {
      a = N
    } else if (is.null(a) && N > 400) {
      a = 400
    }

    if (is.null(b)) {
      b = 1
    }

    if (a < 2 || a > N || round(a)!=a) {
      stop("invalid value for a")
    }

    if (b <= 0 || b > n || round(b)!=b) {
      stop("invalid value for b")
    }

    # find set C of ROC curves

    # falls NULL then b = 1 and a = N
    class_length <- Encoding$lengths[-1]
    s <- ceiling((N - 1) / (a - 1))
    indx_setCa <- seq(1, (1 + (a - 1) * s), s)
    indx_setCb <- which(class_length>n/b)
    indxsetC <- sort(unique(c(indx_setCa, indx_setCb)))

    ncontrols <- (which(duplicated(response) == FALSE) - 1)[-1]
    ncontrols_split <- ncontrols[indxsetC]
    thresholds_split <- thresholds[indxsetC]
    group_length <- c(ncontrols_split) - c(0, ncontrols_split[-length(ncontrols_split)])
    groups <- rep(seq(1, length(group_length),1), group_length)

    weights <- (ncontrols_split * (n - ncontrols_split))
    weights_scaled <- weights/ max(ncontrols*(n-ncontrols))
    # indices to define which transformed problems are used in the computation of uroc

    Ranking_predictor_1 <- rank(predictor, ties.method = "first")
    Classes_predictor <- cumsum(duplicated(sort(predictor))==FALSE)[Ranking_predictor_1]

    Split_classes_predictor <- split(Classes_predictor[1:length(groups)], groups)
    Split_classes_predictor_ordered <- lapply(Split_classes_predictor, function(x){c(0, sort(x))})
    rm(Split_classes_predictor)
    Split_classes_predictor_ordered_diff <- lapply(Split_classes_predictor_ordered, function(x){x[-1] - x[-length(x)]})
    rm(Split_classes_predictor_ordered)


    # compute first roc curve
    order_predictor <- order(predictor, decreasing = TRUE)
    first_threshold <- thresholds_split[1]
    response_binary <- response[order_predictor] >= first_threshold

    dups <- rev(duplicated(rev(predictor[order_predictor])))
    tp <- c(0, cumsum(response_binary == 1)[!dups])
    fp <- c(0, cumsum(response_binary == 0)[!dups])
    truepositive <- rev(tp)
    falsepositive <- rev(fp)
    space_size <- 1000
    InterPoint <- seq(0, 1, (1 / space_size))
    first_roc_curve <- approx(x = fp/ncontrols_split[1], y = tp/(n-ncontrols[1]), xout = InterPoint, method = "linear", ties = "ordered")
    hitrate <- first_roc_curve$y
    # include first plot
    sum_tp_fp <- falsepositive + truepositive
    name <- paste("roc_curve_",1,sep="")
    auc <- round(Trapezoidal(InterPoint, hitrate),2)
    w <- round(weights_scaled[1],2)
    z <- round(thresholds_split[1],2)
    rocm_list = list()
    roc_single <- list(farate = c(0,InterPoint), hitrate = c(0,hitrate), auc = auc, weight = w, threshold = z)
    rocm_list[[name]] <- roc_single

    for (i in 2:length(indxsetC)) {

      m <- length(Split_classes_predictor_ordered_diff[[i]])
      sum_indicator <- rep(seq(m,1,-1), Split_classes_predictor_ordered_diff[[i]])
      sequence_to_change <- length(sum_indicator)
      truepositive[1:sequence_to_change] <- truepositive[1:sequence_to_change] - sum_indicator
      farate <- (sum_tp_fp - truepositive) / ncontrols_split[i]
      roc_approx <- approx(x = rev(farate), y = rev(truepositive) / (n-ncontrols_split[i]), xout = InterPoint, method = "linear", ties = "ordered")
      hitrate <- roc_approx$y

      auc <- round(Trapezoidal(InterPoint, hitrate),2)
      w <- round(weights_scaled[i],2)
      z <- round(thresholds_split[i],2)
      name <- paste("roc_curve_",i,sep="")
      roc_single <- list(farate = c(0,InterPoint), hitrate = c(0,hitrate), auc = auc, weight = w, threshold = z)
      rocm_list[[name]] <- roc_single
    }

    if(object == TRUE) {
      return(rocm_list)
    }

    if(gif == TRUE) {
      uroc_object <- uroc(response, predictor, object = TRUE, plot = FALSE)
      auc <- round(Trapezoidal(uroc_object$Farate, uroc_object$Hitrate), 2)
      rocm_list[["uroc"]] <- list(farate = uroc_object$Farate, hitrate = uroc_object$Hitrate, auc = auc)

      ani.options(loop = 1, interval = 0.1, ...)

      saveGIF({
        for(i in 1:length(rocm_list)) {
          plot(rocm_list[[i]]$farate, rocm_list[[i]]$hitrate, type="l", xlab = "1-Specificity", ylab = "Sensitivity")
          lines(x = c(0,1), y = c(0,1))
          if(i < length(rocm_list)) {text(x=0.75, y=0.2, labels = paste("AUC:", rocm_list[[i]]$auc), adj = 0)
            text(x = 0.0, y = 0.95, labels = paste("z =", rocm_list[[i]]$threshold), adj = 0)
            text(x = 0.2, y = 0.95, labels = paste("w =", rocm_list[[i]]$weight), adj = 0)}
          if(i == length(rocm_list)) {text(x = 0.75, y = 0.2, labels = paste("CPA:", rocm_list[[i]]$auc), adj = 0)
            text(x = 0, y = 0.95, labels = "UROC curve", adj = 0)}
       }}, movie.name = movie.name)
    }
  }
