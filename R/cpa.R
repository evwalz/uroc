#' @title Computes coefficient of predictive ability (CPA).
#' @description This function computes the coefficient of predictive ability which is equalivalent to the area under the UROC curve. Two syntaxes are possible: one object of class "uroc" or two vectors, the response and the predictor.
#' @aliases cpa.default cpa.uroc
#' @param response a numeric vector of real valued responses.
#' @param predictor a numeric vector of the same length as \code{response}, containing real valued predictions for each observation.
#' @param uroc an object of class "uroc" contaning the values of the false alarm rate (1-specificity) and the hitrate (sensitivity) of the UROC curve.
#' @param ... ignored
#' @details The CPA is an asymmetric measure that is linearly related to the correlation between the classes of the response variable and the ranks of the predictor.
#'
#' @importFrom stats cov
#'
#' @return The numeric CPA value.
#' @rdname cpa
#' @export
#'
#' @examples
#'
#' data(longley)
#' response = longley$Employed
#' predictor = longley$GNP
#' cpa(response, predictor)


cpa <- function(...) {
  UseMethod("cpa")
}

#' @rdname cpa
#' @export
cpa.default <- function(response, predictor, ...) {

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

  if (length(unique(response)) < 2) {
    stop("response must have more than one level")
  }

  if (length(response) != length(predictor)) {
    stop("response and predictor should have the same length")
  }

  # order decreasing by response
  response_order <- order(response, decreasing = FALSE)
  response_sort <- response[response_order]
  predictor_sort <- predictor[response_order]

  # compute ranks and classes
  predictor_Rank <- rank(x = predictor_sort, ties.method = "average")
  response_Rank <- rank(x = response_sort, ties.method = "average")
  response_Class <- cumsum(!duplicated(response_Rank))

  # compute cpa using definition of class correlation coefficient
  cpa <- (cov(response_Class, predictor_Rank)/cov(response_Class, response_Rank) + 1) * 0.5

  return(cpa)
}

#' @rdname cpa
#' @export

cpa.uroc <- function(uroc, ...) {

  Farate <- uroc$Farate
  Hitrate <- uroc$Hitrate
  cpa <- Trapezoidal(Farate, Hitrate)

  return(cpa)
}
