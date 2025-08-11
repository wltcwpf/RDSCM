
#' Calculate Spatial Correlation Between Two Locations in California
#'
#' Computes the estimated spatial correlation and its confidence interval between two geographic locations in California, based on a specified oscillator period.
#'
#' @param location_1 A numeric vector of length two. The first element is longitude and the second is latitude.
#' Coordinates must be in decimal degrees using the WGS84 coordinate reference system (EPSG:4326).
#' Note: This function currently only supports locations within California, USA.
#' @param location_2 A numeric vector of length two. The first element is longitude and the second is latitude, using the same CRS as \code{location_1}.
#' @param period A numeric value indicating the oscillator period. Accepted values include \code{-1}, \code{0}, or any value in the interval \code{[0, 10]}.
#' If \code{period = -1}, it refers to peak ground velocity (PGV);
#' if \code{period = 0}, it refers to peak ground acceleration (PGA).
#' @param alpha A numeric value specifying the significance level. Default is \code{0.05}.
#'
#' @return A data frame with two columns:
#' \describe{
#'   \item{\code{estimate}}{The estimated spatial correlation between the two locations.}
#'   \item{\code{interval}}{A character string representing the confidence interval for the estimate, formatted as \code{"[lower, upper]"}.}
#' }
#'
#' @importFrom sf st_point st_sfc st_linestring st_transform st_crs
#' @importFrom stats approx qnorm
#' @importFrom units set_units
#'
#' @export
#'
#' @examples
#' location_1 <- c(-122.4194, 37.7749)  # San Francisco
#' location_2 <- c(-121.8863, 37.3382)  # San Jose
#' spa_cor <- CA_spa_cor(location_1, location_2, period = -1)
CA_spa_cor <- function(location_1, location_2, period, alpha = 0.05){

  if(period == -1){
    polys <- RDSCM::pgv
    polys_2 <- NULL
  }else if(period == 0){
    polys <- RDSCM::pga
    polys_2 <- NULL
  }else if(period > 10 | period < 0){
    return(warning("Exceed applicable period range."))
  }else if(period <=0.1){
    polys <- RDSCM::psa0p1
    polys_2 <- NULL
  }else if(period > 0.1 & period <= 0.3){
    polys <- RDSCM::psa0p1
    polys_2 <- RDSCM::psa0p3
    x <- c(0.1, 0.3)
  }else if(period > 0.3 & period <= 1){
    polys <- RDSCM::psa0p3
    polys_2 <- RDSCM::psa1
    x <- c(0.3, 1)
  }else if(period > 1 & period <= 3){
    polys <- RDSCM::psa1
    polys_2 <- RDSCM::psa3
    x <- c(1, 3)
  }else{
    polys <- RDSCM::psa3
    polys_2 <- NULL
  }

  p1 <- st_point(location_1)
  p2 <- st_point(location_2)
  line <- st_sfc(st_linestring(rbind(p1, p2)), crs = 4326)
  # make the line and the map have the same projected CRS (here: EPSG=3857)
  line_proj  <- st_transform(line,  st_crs(polys))


  segments <- find_intersections(polys = polys, line_proj = line_proj)
  mu <- sum(segments$lamb*as.numeric(segments$seg_length))
  sigma2 <- sum(segments$lamb_sd^2*as.numeric(segments$seg_length)^2)


  if(!is.null(polys_2)){
    segments_2 <- find_intersections(polys = polys_2, line_proj = line_proj)
    mu_2 <- sum(segments_2$lamb*as.numeric(segments_2$seg_length))
    sigma2_2 <- sum(segments_2$lamb_sd^2*as.numeric(segments_2$seg_length)^2)

    # Linear interpolation
    y_mu <- c(mu, mu_2)
    y_sigma2 <- c(sigma2, sigma2_2)
    interp_mu <- approx(log(x), y_mu, xout = log(period))$y
    interp_sigma2 <- approx(log(x), y_sigma2, xout = log(period))$y

    rho_hat <- 0.9*exp(-interp_mu)
    interval <- 0.9*exp(-interp_mu + c(-1,1) * qnorm(1-alpha/2) *sqrt(interp_sigma2))
    return(data.frame(estimate = rho_hat,
                      interval = paste0("[",round(interval[1],3),", ", round(interval[2],3), "]")))
  }

  rho_hat <- 0.9*exp(-mu)
  interval <- 0.9*exp(-mu + c(-1,1) * qnorm(1-alpha/2) *sqrt(sigma2))

  return(data.frame(estimate = rho_hat,
                    interval = paste0("[",round(interval[1],3),", ", round(interval[2],3), "]")))

}
