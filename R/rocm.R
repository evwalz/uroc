#' @title Builds a ROC movie (ROCM)
#' @description This function computes an animated ROC movie.
#' @details The ROC movie can be used to visualize the performance of a real valued foreacsting problem. Therefore, a sequence of png-files is generated and combined into a GIF animation using the external software ImageMagick (\url{https://imagemagick.org/}) and the R package \code{animation}.
#' @param response a numeric vector of real valued responses
#' @param predictor a numeric vector of the same length than \code{response}, containing real valued predictions for each observation
#' @param path folder path
#' @param output name of GIF animation
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
#' rocm(response, predictor, path="/home")}

rocm <-  function(response,
                  predictor,
                  path,
                  output = "animation.gif",
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

  thresholds <- unique(sort(response))[-1]
  N <- length(thresholds)
  n <- length(response)

  if (N == 0) {
    stop("response must have more than one level")
  }

  if (n != length(predictor)) {
    stop("response and predictor should have the same length")
  }

  response_order <- order(response, decreasing=FALSE)
  response <- response[response_order]
  predictor <- predictor[response_order]

  # compute uROC curve
  uroc_curve <- Exact(response, predictor)

  ncontrols <- (which(duplicated(response) == FALSE) - 1)[-1]
  ncases <- (n - ncontrols)
  weights_all <- (ncases * ncontrols)
  weights_scaled <- weights_all / max(weights_all)
  weights_sum <- sum(weights_all)

  pre_order <- order(predictor, decreasing = TRUE)
  predictor <- predictor[pre_order]
  dups <- rev(duplicated(rev(predictor)))

  # set animation options
  animation::ani.options(loop = 1, interval = interval)

  # create png files
  for (i in 1:N) {
    controls <- ncontrols[i]
    cases <- (n - controls)
    response_binary <- c(rep(0, controls), rep(1, cases))[pre_order]
    hitrate <- c(0, cumsum(response_binary == 1)[!dups] / cases)
    farate <- c(0, cumsum(response_binary == 0)[!dups] / controls)

    auc <- Trap(farate, hitrate)

    name <- rename(i)

    png(filename = paste(path, "/", name, sep=""), width = 1920, height = 1080, units = "px")
    par(mgp = c(5, 1, 0), mar = c(9, 9, 9, 9))
    plot(x = farate, y = hitrate, las = 1, type = "l", lwd = 2, xlim = c(0, 1), ylim = c(0, 1) ,
         xlab = "False alarm rate", ylab = "Hit rate", cex.lab = 3, cex.axis = 2, cex.main = 2.5, ...)
    lines(x = c(0, 1), y = c(0, 1), lty = 2)
    text(x = 0, y = 1, labels = paste("Threshold:", format(round(thresholds[i], 2), nsmall = 2),
                             "\nWeight:", format(weights_scaled[i], nsmall = 2)),
         adj = c(0, 1), col = "black" ,cex = 2.5)
    text(x = 0.6, y = 0.4, labels = paste("AUC:", format(round(auc, 2), nsmall = 2)), cex = 2.5)
    dev.off()
  }

  # create png of uROC curve
  name <- rename(N+1)
  uFar <- uroc_curve$Far
  uHitrate <- uroc_curve$Hitrate
  cpa_value <- Trap(uFar, uHitrate)
  png(filename = paste(path, "/", name, sep = ""), width = 1920, height = 1080, units = "px")
  par(mgp = c(5, 1, 0), mar = c(9, 9, 9, 9))
  plot(x = uFar, y = uHitrate, las = 1, type = "l", lwd = 2, xlim = c(0,1), ylim = c(0,1) ,
       xlab = "False alarm rate", ylab = "Hit rate", cex.lab = 3, cex.axis = 2, cex.main = 2.5)
  lines(x = c(0, 1), y = c(0, 1), lty = 2)
  text(x = 0.6, y = 0.4, paste("CPA:", format(round(cpa_value, 2), nsmall = 2)), cex = 2.5)
  dev.off()

  # create GIF animation
  animation::im.convert(paste(path, "/", "*.png", sep=""), output = paste(path, "/", output,sep = ""), convert = convert,
             cmd.fun = cmd.fun, extra.opts = "",clean = clean)

}





