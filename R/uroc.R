#' @title Computes a uROC curve
#' @description This function builds a uROC curve and returns a "uroc" object, a list of class "uroc".
#' @details There are 3 different algorithms available to create a uROC curve. The input argument \code{algo="exact"} computes the exact uROC curve. Using \code{algo="approx1"} or \code{algo="approx2"} generates an approximation to the uROC curve by computing the y-values of the curve only on specific x-values. The x-values are equidistant over the interval [0,1] and the number of x-values used in the computation can be set by \code{space.size}. Calling \code{algo="approx1"} generates an approximation with a correction for ties in the predictor vaiable whereas \code{algo="approx2"} ignores ties in the predictor variable but results in a faster computation. Therefore, it is recommended to either use \code{algo="exact"} or \code{algo="approx2"} if the input vector for \code{predictor} contains a lot of tied values. If the type of algorithm is not specified, the \code{\link{uroc}} function choses one of the three versions based on the input arguments in \code{response} and \code{predictor}.
#' @param response a numeric vector of real valued responses
#' @param predictor a numeric vector of the same length than \code{response}, containing real valued predictions for each observation
#' @param object if TRUE an object of type uroc is returned containg the false alarm rate and the hitrate of the uROC curve
#' @param plot plot the uROC curve? if \code{FALSE} the curve is not displayed
#' @param algo optional argument to select an algorithm for the computation of the uROC curve. See Details.
#' @param space.size optional argument to set the number of x-values for which the corresponding value in the approximation algorithm for the uROC curve is computed. It is the inverse value of the distance between equidistant points within the interval [0,1]
#'
#' @importFrom graphics text lines
#'
#' @return If \code{object = TRUE} this function returns a list of class "uroc".
#' @export
#'
#' @examples
#' data(longley)
#' response = longley$Employed
#' predictor = longley$GNP
#' uroc(response, predictor)

uroc <- function(response,
                 predictor,
                 object = FALSE,
                 plot = TRUE,
                 algo = NULL,
                 space.size = NULL) {

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

  N <- length(unique(response))
  n <- length(response)

  if (N == 0) {
    stop("response must have more than one level")
  }

  if (n != length(predictor)) {
    stop("response and predictor should have the same length")
  }

  if (!is.null(algo) && algo != "exact" && algo != "approx1" && algo != "approx2") {
      stop("invalid argument for algo")
  }


  cpa_exact <- cpa(response, predictor)
  response_order <- order(response, decreasing=FALSE)
  response <- response[response_order]
  predictor <- predictor[response_order]


  if (is.null(algo)) {
      if (N <= 400) {
      algo <- "exact"
    } else if (n<=1000) {
      algo <- "approx1"
    } else if (n <= 10000 & N <= 50) {
      algo <- "approx1"
    } else if (max(rle(duplicated(sort(predictor)))$length)<200) {
      algo <- "approx1"
    } else {
      algo <- "approx2"
    }
  }

  if (algo == "exact") {
    uROC <- Exact(response, predictor)
  }  else if (algo == "approx1") {
    uROC <- Approx1(response, predictor, space.size)
  } else if (algo == "approx2") {
    uROC <- Approx2(response, predictor, space.size)
  }

  Farate <- uROC$Farate
  Hitrate <- uROC$Hitrate

  cpa_approx <- Trap(Farate, Hitrate)

  if (abs(cpa_approx-cpa_exact) > 0.05) {
    warning("poor approximation of uROC curve")
  }

  #plot
  if (plot) {
    cpa_value <- round(cpa_approx, 2)
    plot(x = Farate, y = Hitrate, type = "l", xlab = "False alarm rate", ylab = "Hitrate")
    lines(x = c(0, 1), y = c(0, 1), lty = 2)
    text(x = 0.6, y = 0.4, labels = paste("CPA:",cpa_value))
  }

  #object
  if (object) {
    uroc <- list(Farate = Farate, Hitrate = Hitrate)
    class(uroc) <- "uroc"
    return(uroc)
  }

}
