#' @title Computes a uROC curve
#' @description This function builds a uROC curve and returns a "uroc" object, a list of class "uroc".
#' @details There are 2 different algorithms available to create a uROC curve. The input argument \code{algo="exact"} computes the exact uROC curve and \code{algo="approx"} generates an approximation to the uROC curve by computing the y-values of the curve only on specific x-values. The x-values are equidistant points over the interval [0,1] and the number of x-values can be set by \code{space.size}. If the type of algorithm is not specified, the \code{\link{uroc}} function choses one of the two versions based on the input arguments in \code{response} and \code{predictor}. 
#' @param response a numeric vector of real valued responses
#' @param predictor a numeric vector of the same length than \code{response}, containing real valued predictions for each observation
#' @param object if TRUE an object of type uroc is returned containg the false alarm rate and the hitrate of the uROC curve
#' @param plot plot the uROC curve? if \code{FALSE} the curve is not displayed
#' @param algo optional argument to select an algorithm for the computation of the uROC curve. See Details.
#' @param space.size optional argument to set the number of x-values for which the corresponding value in the approximation algorithm for the uROC curve is computed. It is the inverse value of the distance between equidistant points within the interval [0,1]
#'
#' @importFrom graphics text lines
#' @importFrom stats approxfun 
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
                 split = NULL,
                 space.size = NULL,
                 algo = NULL) {

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


  response_order <- order(response, decreasing=FALSE)
  response <- response[response_order]
  predictor <- predictor[response_order]
  
  
  if(algo=='exact'){
    uroc_object <- compute_uroc_exact(response, predictor,n,N) 
  } else {
    uroc_object <- compute_uroc_approx(response, predictor,n,N)  
  }
  

  Farate <- uroc_object$Farate
  Hitrate <- uroc_objectHitrate

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
