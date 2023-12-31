% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/BARModel.R
\name{fastlognormdens}
\alias{fastlognormdens}
\title{#' @export
fastnormdens <- function(x, mean = 0, sd = 0){
  (1/(sd*sqrt(2*pi)))*exp(-0.5*((x-mean)/sd)^2)
}}
\usage{
fastlognormdens(x, mean = 0, sd = 0)
}
\description{
#' @export
fastnormdens <- function(x, mean = 0, sd = 0){
  (1/(sd*sqrt(2*pi)))*exp(-0.5*((x-mean)/sd)^2)
}
}
