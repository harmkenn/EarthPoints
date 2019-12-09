#' UTM2MGRS2LL
#'
#' Convert UTN into MGRS and Lat Long
#' @param Easting numeric 6 digit Easting
#' @param Northing numeric 7 digit Northing
#' @param gzn numeric 2 digit zone
#' @param hem character Hemisphere "N" or "S"
#' @import plyr
#' @import pracma
#' @importFrom pracma mod
#' @export
#' @examples
#' UTM2MGRS2LL(530000,3910000,11,"N")
#'
## Update the details for the return value
#' @return This function returns a \code{data.frame} including columns:
#' \itemize{
#'  \item Easting
#'  \item Northing
#'  \item Zone
#'  \item Hemisphere
#'  \item Latitude
#'  \item Longitude
#'  \item MGRS
#' }
#'
#' @author Ken Harmon <harmkenn@gmail.com>

UTM2MGRS2LL <- function(Easting,Northing,gzn,hem){
  # Outside Functions I used
  mod <- pracma::mod
  match_df <- plyr::match_df
  f_pad_zero <- numform::f_pad_zero
  # My Stuff
  Easting <- round(Easting,0)
  Northing <- round(Northing,0)

  sin1 <- pi/(180*3600) #One Second
  a <- 6378137 # Equitorial Radius
  b <- 6356752.31424518 #Polar Radius
  k0 <- .9996 #Scalar Factor Constant
  e1 <- sqrt(a^2-b^2)/b
  e1sq <- e1^2
  c <- a^2/b

  NfEQ <- Northing
  if (hem == "S"){NfEQ <- Northing - 10000000}
  Fi <- (NfEQ)/(6366197.724*k0)
  Ni <- (c/(1+e1sq*(cos(Fi))^2)^(1/2))*k0
  Czone <- 6*gzn-183
  dln <- (Easting-500000)/Ni
  A1 <- sin(2*Fi)
  A2 <- A1*(cos(Fi))^2
  J2 <- Fi+(A1/2)
  J4 <- (3*J2+A2)/4
  J6 <- (5*J4+A2*(cos(Fi))^2)/3
  alfa <- 3/4*e1sq
  beta <- 5/3*alfa^2
  gamma <- 35/27*alfa^3
  Bfi <- k0*c*(Fi-(alfa*J2)+(beta*J4)-(gamma*J6))
  BB <- (NfEQ-Bfi)/Ni
  zeta <- ((e1sq*dln^2)/2)*(cos(Fi))^2
  Xi <- dln*(1-(zeta/3))
  Eta <- BB*(1-zeta)+Fi
  ShXi <- (exp(Xi)-exp(-Xi))/2
  dLam <- atan(ShXi/cos(Eta))
  Tau <- atan(cos(dLam)*tan(Eta))
  FiR <- Fi+(1+e1sq*(cos(Fi))^2-(3/2)*e1sq*sin(Fi)*cos(Fi)*(Tau-Fi))*(Tau-Fi)
  lat <- FiR/pi*180
  lng <- dLam/pi*180+Czone

  ## Now let's turn UTM into MGRS
  gze <- "Odd"
  if((gzn %% 2) == 0) {gze <- "Even"}

  #Latitude gridzone letters
  gznls <- c("C","D","E","F","G","H","J","K","L","M",
             "N","P","Q","R","S","T","U","V","W","X")
  lzc <- (lat+80-(lat+80)%%8)/8+1 # Latitude letter zone count
  gzl <- gznls[lzc]

  # 100km grid Square easting letter
  esqls <- as.matrix(rbind(c("A","B","C","D","E","F","G","H"),
                           c("J","K","L","M","N","P","Q","R"),
                           c("S","T","U","V","W","X","Y","Z")))
  colnames(esqls) <- 1:8
  rownames(esqls) <- c(1,2,0)

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

  as.data.frame(rbind(c("Easting" = Easting,"Northing" = Northing, "Zone" = gzn, "Hemisphere" = hem, "Latitude" = round(lat,2), "Longitude"=round(lng,2), "MGRS" = MGRS)))
}
