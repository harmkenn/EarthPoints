#' polar_LL
#'
#' Compute a polar shot from a point and a vector
#' @param lat1d numeric launch latitude in degrees
#' @param lng1d numeric launch longitude in degrees
#' @param dist numeric distance in km
#' @param shotd numeric bearing in degrees
#' @export
#' @examples
#' polar_LL(-30,-40,10,30)
#'
## Update the details for the return value
#' @return This function returns a \code{data.frame} including columns:
#' \itemize{
#'  \item Launch Latitude
#'  \item Launch Longitude
#'  \item Distance in KM
#'  \item Launch Bearing
#'  \item Landing Latitude
#'  \item Landing Longitude
#'  \item Landing Bearing
#'  \item Midpoint Latitude
#'  \item Midpoint Longitude
#' }
#'
#' @author Ken Harmon <harmkenn@gmail.com>

polar_ll <- function(lat1d,lng1d,dist,shotd){
  #in Radians
  lat1r <- lat1d*pi/180
  lng1r <- lng1d*pi/180
  shotr <- shotd*pi/180
  er <- 6371 #Earth Radius in km
  delta <- dist/er

  lat2r <- asin(sin(lat1r)*cos(delta)+cos(lat1r)*sin(delta)*cos(shotr))
  lng2r <- lng1r + atan2(sin(shotr)*sin(delta)*cos(lat1r),
                         cos(delta)-sin(lat1r)*sin(lat2r))

  lat2d <- lat2r*180/pi
  lng2d <- lng2r*180/pi

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

  data.frame("Launch Latitude" = lat1d,
             "Launch Longitude" = lng1d,
             "Distance km" = dist,
             "Launch Bearing" = paste(round(shotd,2)),
             "Landing Latitude" = lat2d,
             "Landing Longitude" = lng2d,
             "Landing Bearing" = round(impd,2),
             "Midpoint Lat" = round(latmd,2),
             "Midpoint Lng" = round(lngmd,2))
}

