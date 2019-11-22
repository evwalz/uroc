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
#' rocm(response, predictor, path="/home")}

rocm <-  function(response,
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

  # sort response and predictor
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

  # compute controls
  ncontrols <- (which(duplicated(response) == FALSE) - 1)[-1]

  # select controls by set C and compute weights
  ncontrols_split <- ncontrols[indxsetC]
  ncases_split <- (n - ncontrols_split)
  weights_split <- (ncases_split * ncontrols_split)
  weights_scaled <- weights_split / max(weights_split)
  thresholds_split <- thresholds[indxsetC]

  # order by predictor to compute farate and hitrate
  predictor_order <- order(predictor, decreasing = TRUE)
  predictor <- predictor[predictor_order]
  dups <- rev(duplicated(rev(predictor)))
  response <- response[predictor_order]
  #
  InterPoint <- seq(0,1,1/1000)


  List_of_hitrate <- lapply(ncontrols_split, function(x, n, dups, predictor_order) {
    cases = n-x
    binary <- c(rep(0, x), rep(1, cases))[predictor_order]
    hitrate = c(0, (cumsum(binary == 1) / sum(x == binary))[!dups])
    farate = c(0, (cumsum(binary == 0) / sum(binary == 0))[!dups])
    hit = approx(farate, hitrate, InterPoint, method = "linear", ties="ordered" )
    return(hit$y)},n, dups, predictor_order)

  List_of_auc <- lapply(List_of_hitrate, function(x, InterPoint) {round(Trapezoidal(InterPoint, x),2)}, InterPoint)

  animation::ani.options(loop = 1, interval = interval)

  # create png files
  for (i in 1:length(indxsetC)) {

    auc <- List_of_auc[[i]]

    w = round(weights_scaled[i],2)
    z = round(thresholds[i],2)

    name <- rename_figure(i)

    png(filename = paste(path, "/", name, sep=""), res=72, width = 480, height = 480)
    plot(x = InterPoint, y = List_of_hitrate[[i]], las = 1, type = "l", lwd = 2, xlim = c(0, 1), ylim = c(0, 1) ,
         xlab = "False alarm rate", ylab = "Hit rate", cex.lab = 3, cex.axis = 2, cex.main = 2.5)#, ...)
    lines(x = c(0, 1), y = c(0, 1), lty = 2)
    text(x = 0, y = 1, labels = paste("Threshold:", format(round(thresholds[i], 2), nsmall = 2),
                             "\nWeight:", format(w, nsmall = 2)),
         adj = c(0, 1), col = "black" ,cex = 2.5)
    text(x = 0.6, y = 0.4, labels = paste("AUC:", format(auc, nsmall = 2)), cex = 2.5)
    dev.off()
  }

  # create png of uROC curve
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





