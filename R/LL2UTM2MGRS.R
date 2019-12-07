#' LL2UTM2MGRS
#'
#' Convert Lat Long into MGRS through UTM
#' @param lat numeric latitude
#' @param lng numeric longitude
#' @export
#' @examples
#' LL2UTM2MGRS(-30.4,-40.8)
#'
## Update the details for the return value
#' @return This function returns a \code{data.frame} including columns:
#' \itemize{
#'  \item Latitude 
#'  \item Longitude
#'  \item Easting
#'  \item Northing
#'  \item Zone
#'  \item Hemisphere
#'  \item MGRS
#' }
#'
#' @author Ken Harmon <harmkenn@gmail.com>

LL2UTM2mgrs<-function(lat,lng){
# Outside Functions I used
  mod <- pracma::mod
  match_df <- plyr::match_df
  f_pad_zero <- numform::f_pad_zero
# My Stuff
  hem <- "N"
  if (lat < 0){hem <- "S"}


  # 100km grid Square easting letter
  esqls <- as.matrix(rbind(c("A","B","C","D","E","F","G","H"),
                           c("J","K","L","M","N","P","Q","R"),
                           c("S","T","U","V","W","X","Y","Z")))
  colnames(esqls) <- 1:8
  rownames(esqls) <- c(1,2,0)


  lngr <- lng*pi/180 # Longitude in radians
  latr <- lat*pi/180 # Latitude in radians

  # Lets go find Easting

  sin1 <- pi/(180*3600) #One Second
  a <- 6378137 # Equitorial Radius
  b <- 6356752.31424518 #Polar Radius
  k0 <- .9996 #Scalar Factor Constant
  gzn <- floor(1/6*lng)+31 #Longitude Zone
  Czone <- 6*gzn - 183 #Longitude of the center of the zone
  dlng <- lng - Czone # Longitude from the center of the zone
  p <- dlng*3600/10000 #Hecta seconds?
  e <- sqrt(1-(b/a)^2) #eccentricity
  e1 <- sqrt(a^2-b^2)/b
  e1sq <- e1^2
  c <- a^2/b
  nu <- a/sqrt(1-(e*sin(latr))^2) #r curv 2
  Kiv <- nu*cos(latr)*sin1*k0*10000 #Coef for UTM 4
  Kv <- (sin1*cos(latr))^3*(nu/6)*
    (1-tan(latr)^2+e1sq*cos(latr)^2)*k0*10^12 #Coef for UTM 5
  Easting <- 500000+Kiv*p+Kv*p^3

  # Now let's go find Northing
  n <- (a-b)/(a+b)
  A0 <- a*(1-n+(5*n^2/4)*(1-n)+(81*n^4/64)*(1-n)) # Meridional Arc Length
  B0 <- (3*a*n/2)*(1-n-(7*n^2/8)*(1-n)+55*n^4/64) # Meridional Arc Length
  C0 <- (15*a*n^2/16)*(1-n+(3*n^2/4)*(1-n)) # Meridional Arc Length
  D0 <- (35*a*n^3/48)*(1-n+11*n^2/16) # Meridional Arc Length
  E0 <- (315*a*n^4/51)*(1-n) # Meridional Arc Length
  S <- A0*latr - B0*sin(2*latr) + C0*sin(4*latr) - D0*sin(6*latr) + E0*sin(8*latr) # Meridional Arc
  Ki <- S*k0 #Coef for UTM 1
  Kii <- nu*sin(latr)*cos(latr)*sin1^2*k0*100000000/2 #Coef for UTM 2
  Kiii <- ((sin1^4*nu*sin(latr)*cos(latr)^3)/24)*(5-tan(latr)^2+9*e1sq*cos(latr)^2*cos(latr)^4)*k0*10^16 #Coef for UTM 2
  Northing <- Ki + Kii * p^2 + Kiii * p^4
  if (lat < 0) {Northing <- 10000000 + Northing} # In the Southern Hemisphere is Northing is measured from the south pole instead of from the equator

  Easting <- round(Easting,0)
  Northing <- round(Northing,0)


  ## Now let's turn UTM into MGRS
  gze <- "Odd"
  if((gzn %% 2) == 0) {gze <- "Even"}

  #Latitude gridzone letters
  gznls <- c("C","D","E","F","G","H","J","K","L","M",
             "N","P","Q","R","S","T","U","V","W","X")
  lzc <- (lat+80-(lat+80)%%8)/8+1 # Latitude letter zone count
  gzl <- gznls[lzc]

  gsen <-  plyr::round_any(Easting,100000,floor)/100000 #Grab off the first digit off the Easting
  gsnn <-  plyr::round_any(Northing,100000,floor)/100000 #Grab off the first two digits off the Northing

  emod <- mod(gzn,3) # gzn mod 3 for 100km lookup
  if (emod == 0) {emod <- 3} #turn mod zero into row 3
  gsel <- esqls[emod,gsen] #Here is the 100km easting grid letter

  alln <- data.frame(hem = hem, gze = gze, gzl = gzl, gsnl = NA, gsnn = gsnn)
  fulln<-match_df(northdes,alln,on=c("hem","gze","gzl","gsnn"))

  gsnl <- as.character(fulln[1,4])

  east <- plyr::round_any(Easting,1) - plyr::round_any(Easting,100000,floor) #Snatch the right 5 digits off the easting
  easting5 <- f_pad_zero(round(east,0), width = 5, pad.char = "0") #fill the front with zeros if needed
  north <- plyr::round_any(Northing,1) - plyr::round_any(Northing,100000,floor) #Snatch the right 5 digits off the northing
  northing5 <- f_pad_zero(north, width = 5, pad.char = "0") #fill the front with zeros if needed

  utmz <- f_pad_zero(gzn, width = 2, pad.char = "0") #pad the zone with a zero, if needed
  MGRS <- paste(utmz,gzl,gsel,gsnl,easting5,northing5,sep = "") #Paste all six pieces together

  final <- as.data.frame(rbind(c(lat,lng, Easting,Northing,utmz,hem,MGRS)))
  colnames(final) <- c("Latitude","Longitude","Easting","Northing","zone","Hemisphere","MGRS")


  return(final)
}
