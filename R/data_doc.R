#' KidneyMimic data set
#'
#' @name KidneyMimic
#' @docType data
#' @keywords datasets
#' @description A simulated data frame with 200 individuals, distributed in 40 clinics (clusters) with 5 individuals per clinic. Simulated dataset that mimic a kidney result study. In this simulated dataset it is hypothetically suggested the following causes of outcome: event of interest: death due to kidney failure; dependent censoring: if the patient received a transplant; administrative censoring: all other outcomes.
#' @format A data frame with 200 rows and 9 variables:
#' \itemize{
#'   \item ident: Variable that indicates the clinic (cluster) of the patient.
#'   \item time: time observed until the occurrence of the outcome.
#'   \item event: variable that indicates the occurrence of the event of interest, event=1 if the event of interest occurs.
#'   \item x1: covariate 1, generated from a Uniform Distribution. Can denote, for example, a exam result.
#'   \item x2: covariate 2, generated from a Binomial Distribution. Can denote, for example, an treatment
#'   \item x3: covariate 3, generated from a Normal Distribution. can denote, for example, an standardized age.
#'   \item cens: variable that indicates the outcome, cens=1 if the event of interest occurred (death due to kidney failure); cens=2 if the  dependent censoring occurred (patient received a transplant); cens=3 if administrative censoring.
#'   \item delta_t: indicator function of the event of interest.
#'   \item delta_c: indicator function of the dependent censoring.
#' }
NULL
