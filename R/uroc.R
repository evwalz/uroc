#' @title Computes a UROC curve
#' @description This function builds a UROC curve and returns a "uroc" object, a list of class "uroc".
#' @details There are 2 different algorithms available to create a UROC curve. The default option is  \code{algo="approx"} which generates an approximation to the UROC curve by using linear interpolation of each ROC curve. To reduce computation time the paramter \code{split} can be specified to select a subset of ROC curves in the computation. The input argument \code{algo="exact"} computes the exact UROC curve and should only be used on small data.
#' @param response a numeric vector of real valued responses
#' @param predictor a numeric vector of the same length than \code{response}, containing real valued predictions for each observation
#' @param object if TRUE an object of type uroc is returned containg the false alarm rate and the hitrate of the UROC curve
#' @param plot plot the UROC curve? if \code{FALSE} the curve is not displayed
#' @param algo optional argument to select an algorithm for the computation of the UROC curve. See Details.
#' @param split a integer value with a default of \code{split = 1}. Computes uroc curve by considering only a subset of all N-1 available ROC curves to reduce computation time. The split parameter defines the distance between a set of equidistant indices which are then used to select particular ROC curves among the N-1.
#'
#' @importFrom graphics text lines
#' @importFrom stats approx
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
                 algo = "approx",
                 split = 1) {

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

  if (!is.null(algo) && algo != "exact" && algo != "approx") {
      stop("invalid argument for algo")
  }

  split <- floor(split)
  if (split > N || split < 0) {
    stop("invalid value for split")
  }

  if (!is.logical(object)) {
    stop("invalid input for object")
  }

  if (!is.logical(plot)) {
    stop("invalid input for object")
  }

  response_order <- order(response, decreasing=FALSE)
  response <- response[response_order]
  predictor <- predictor[response_order]


  if(algo=='exact'){
    uroc_object <- compute_uroc_exact(response, predictor,n,N)
  } else {
    uroc_object <- compute_uroc_approx(response, predictor,n,N, split)
  }


  Farate <- uroc_object$Farate
  Hitrate <- uroc_object$Hitrate

  cpa_approx <- Trapezoidal(Farate, Hitrate)

  #plot
  if (plot) {
    cpa_value <- round(cpa_approx, 2)
    plot(x = Farate, y = Hitrate, type = "l", xlab = "1-Specificity", ylab = "Sensitivity")
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
