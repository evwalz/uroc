#' @title Builds a ROC movie (ROCM)
#' @description This function computes an animated ROC movie.
#' @details The ROC movie can be used to visualize the performance of a real valued foreacsting problem. Therefore, a sequence of png-files is generated and combined into a GIF animation using the external software ImageMagick (\url{https://imagemagick.org/}) and the R package \code{animation}.
#' @param response a numeric vector of real valued responses
#' @param predictor a numeric vector of the same length than \code{response}, containing real valued predictions for each observation
#' @param path folder path
#' @param output name of GIF animation
#' @param b default is \code{b=100}
#' @param a default is \code{a=400}
#' @param clean if FALSE png files are not deleted
#' @param convert convert command for the function  \link[animation]{im.convert}
#' @param cmd.fun a function to invoke OS command in  \link[animation]{im.convert}
#' @param interval a postve number to set the time interval of the animation (unit in second) in \link[animation]{ani.options}
#' @param ... plotting arguments
#'
#' @importFrom animation ani.options im.convert
#' @importFrom grDevices png dev.off
#' @importFrom graphics par plot lines text
#'
#' @return GIF animation
#' @export
#'
#' @examples
#' \dontrun{
#' data(longley)
#' response = longley$Employed
#' predictor = longley$GNP
#' rocm2(response, predictor, path="/home")}

rocm2 <-  function(response,
                  predictor,
                  path,
                  output = "animation.gif",
                  b = 100,
                  a = 400,
                  clean = TRUE,
                  convert = "convert",
                  cmd.fun = if (.Platform$OS.type == "windows") shell else system,
                  interval = 0.1,
                  ...) {

  if (dir.exists(path) == FALSE) {
    stop(paste("The directory", path, "does not exist"))
  }

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

  # find set C of ROC curves
  class_length <- Encoding$lengths[-1]
  s <- ceiling(N/a)
  indx_setCa <- seq(1,N,s)
  indx_setCb <- which(class_length>n/b)
  indxsetC <- sort(unique(c(indx_setCa, indx_setCb)))

  ncontrols <- (which(duplicated(response) == FALSE) - 1)[-1]
  ncontrols_split <- ncontrols[indxsetC]
  group_length <- c(ncontrols_split) - c(0, ncontrols_split[-length(ncontrols_split)])
  groups <- rep(seq(1, length(group_length),1), group_length)

  weights <- sum(ncontrols_split * (n - ncontrols_split))
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
  first_threshold <- min(response)
  response_binary <- response[order_predictor] > first_threshold

  dups <- rev(duplicated(rev(predictor[order_predictor])))
  tp <- c(0, cumsum(response_binary == 1)[!dups])
  fp <- c(0, cumsum(response_binary == 0)[!dups])
  truepositive <- rev(tp)
  falsepositive <- rev(fp)
  space_size <- 1000
  InterPoint <- seq(0, 1, (1 / space_size))
  first_roc_curve <- approx(x = fp/ncontrols[1], y = tp/(n-ncontrols[1]), xout = InterPoint, method = "linear", ties = "ordered")
  hitrate <- first_roc_curve$y
  # include first plot
  sum_tp_fp <- falsepositive + truepositive
  animation::ani.options(loop = 1, interval = interval)

  for (i in 2:length(indxsetC)) {

    m <- length(Split_classes_predictor_ordered_diff[[i]])
    sum_indicator <- rep(seq(m,1,-1), Split_classes_predictor_ordered_diff[[i]])
    sequence_to_change <- length(sum_indicator)
    truepositive[1:sequence_to_change] <- truepositive[1:sequence_to_change] - sum_indicator
    farate <- (sum_tp_fp - truepositive) / ncontrols_split[i]
    roc_approx <- approx(x = rev(farate), y = rev(truepositive) * ncontrols_split[i], xout = InterPoint, method = "linear", ties = "ordered")
    hitrate <- roc_approx$y

    auc <- Trapezoidal(InterPoint, hitrate)
    w <- round(weights_scaled[i],2)
    z <- round(thresholds[i],2)

    name <- rename_figure(i)

    png(filename = paste(path, "/", name, sep=""), res=72, width = 480, height = 480)
    plot(x = InterPoint, y = hitrate, las = 1, type = "l", lwd = 2, xlim = c(0, 1), ylim = c(0, 1) ,
         xlab = "False alarm rate", ylab = "Hit rate", cex.lab = 3, cex.axis = 2, cex.main = 2.5)#, ...)
    lines(x = c(0, 1), y = c(0, 1), lty = 2)
    text(x = 0, y = 1, labels = paste("Threshold:", format(round(thresholds[i], 2), nsmall = 2),
                                      "\nWeight:", format(w, nsmall = 2)),
         adj = c(0, 1), col = "black" ,cex = 2.5)
    text(x = 0.6, y = 0.4, labels = paste("AUC:", format(auc, nsmall = 2)), cex = 2.5)
    dev.off()

  }

  name <- rename_figure(length(indxsetC)+1)
  uroc_curve <- uroc(response, predictor, plot = FALSE, object = TRUE)
  uFar <- uroc_curve$Farate
  uHitrate <- uroc_curve$Hitrate
  cpa_value <- Trapezoidal(uFar, uHitrate)
  png(filename = paste(path, "/", name, sep=""))
  plot(x = uFar, y = uHitrate, las = 1, type = "l", lwd = 2, xlim = c(0,1), ylim = c(0,1) ,
       xlab = "False alarm rate", ylab = "Hit rate", cex.lab = 3, cex.axis = 2, cex.main = 2.5)
  lines(x = c(0, 1), y = c(0, 1), lty = 2)
  text(x = 0.6, y = 0.4, paste("CPA:", format(round(cpa_value, 2), nsmall = 2)), cex = 2.5)
  dev.off()

  # create GIF animation
  animation::im.convert(paste(path, "/", "*.png", sep=""), output = paste(path, "/", output,sep = ""), convert = convert,
                        cmd.fun = cmd.fun, extra.opts = "",clean = clean)

}
