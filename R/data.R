#' A simulated digital elevation model
#' 
#' A dataset mimicking a DEM of ground surface data to serve as an example for calculations in \code{lithosphere}. Resolution is 10.
#' 
#' @format 
#' \describe{
#' \item{x}{x coordinate}
#' \item{y}{y coordinate}
#' \item{z}{z coordinate}
#' }
#' 
#' @source made with \code{\link[NLMR]{nlm_fbm}}, see \href{https://ropensci.github.io/NLMR/}{NLMR website}
"synthetic_dem_1"


#' A set of simulated well markers
#' 
#' A dataset of 50 fictional wells and associated stratigraphy to serve as an example for calculations in \code{lithosphere}
#' #' @format 
#' \describe{
#' \item{well_name}{a fictional name for the well to serve as unique identifier}
#' \item{x}{x coordinate}
#' \item{y}{y coordinate}
#' \item{z}{z coordinate of the bottom of a stratigraphic unit}
#' \item{stratigraphy}{a unified geological description of the lithology at xyz}
#' }
#' 
#' @source made with \code{\link[NLMR]{NLMR-package}}, see \href{https://ropensci.github.io/NLMR/}{NLMR website}
"synthetic_welldata_1"