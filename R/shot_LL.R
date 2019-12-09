#' shot_ll
#'
#' Compute the shot from one point to another
#' @param lat1d numeric Launch latitude
#' @param lng1d numeric Launch longitude
#' @param lat2d numeric Landing latitude
#' @param lng2d numeric Landing longitude
#' @export
#' @examples
#' shot_ll(-30,-40,40,30)
#'
## Update the details for the return value
#' @return This function returns a \code{data.frame} including columns:
#' \itemize{
#'  \item Launch Latitude
#'  \item Launch Longitude
#'  \item Landing Latitude
#'  \item Landing Longitude
#'  \item Distance in KM
#'  \item Launch Bearing
#'  \item Landing Bearing
#'  \item Midpoint Latitude
#'  \item Midpoint Longitude
#' }
#'
#' @author Ken Harmon <harmkenn@gmail.com>

shot_ll <- function(lat1d,lng1d,lat2d,lng2d){
  #in Radians
  lat1r <- lat1d*pi/180
  lng1r <- lng1d*pi/180
  lat2r <- lat2d*pi/180
  lng2r <- lng2d*pi/180

  beta <- lng2r - lng1r

  #Central angle between points
  theta <- acos(cos(lat1r)*cos(lat2r)*cos(beta)+sin(lat1r)*sin(lat2r))
  er <-6371 #mean radius of the earth in km

  #Distance between points
  dist <- theta*er

  #Bearing of the shot in radians
  shotr <- atan2(sin(lng2r-lng1r)*cos(lat2r),
                 cos(lat1r)*sin(lat2r)-sin(lat1r)*cos(lat2r)*cos(lng2r-lng1r))
  #Bearing of the shot in degrees
  shotd <- shotr*180/pi

  #Bearing of the impact in radians
  impr <- pi + atan2(sin(lng1r-lng2r)*cos(lat1r),
                     cos(lat2r)*sin(lat1r)-sin(lat2r)*cos(lat1r)*cos(lng1r-lng2r))
  #Bearing of the impact in degrees
  impd <- impr*180/pi

  # Compute midpoint
  Bx <- cos(lat2r)*cos(lng2r-lng1r)
  By <- cos(lat2r)*sin(lng2r-lng1r)
  latmr <- atan2(sin(lat1r) + sin(lat2r), sqrt((cos(lat1r)+Bx)^2+By^2))
  lngmr <- lng1r + atan2(By, cos(lat1r) + Bx)

  latmd <- latmr*180/pi
  lngmd <- lngmr*180/pi


  data.frame("Launch Latitude" = lat1d, "Launch Longitude" = lng1d,
             "Landing Latitude" = lat2d, "Landing Longitude" = lng2d,
             "Distance km" = dist,
             "Launch Bearing" = paste(round(shotd,2)),
             "Landing Bearing" = paste(round(impd,2)),
             "Midpoint Lat" = round(latmd,2),
             "Midpoint Lng" = round(lngmd,2))
}

