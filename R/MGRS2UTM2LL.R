#' MGRS2UTM2LL
#'
#' Convert MGRS into Lat Long through UTM
#' @param MGRS String of 15 characters of Military Grid Reference
#' @import stringr
#' @import plyr
#' @import pracma
#' @importFrom pracma mod
#' @export
#' @examples
#' MGRS2UTM2LL("11SNV3000010000")
#'
## Update the details for the return value
#' @return This function returns a \code{data.frame} including columns:
#' \itemize{
#'  \item MGRS
#'  \item Easting
#'  \item Northing
#'  \item Zone
#'  \item Hemisphere
#'  \item Latitude
#'  \item Longitude
#' }
#'
#' @author Ken Harmon <harmkenn@gmail.com>

MGRS2UTM2LL <-function(MGRS){

  # Outside Functions I used
  mod <- pracma::mod
  match_df <- plyr::match_df
  f_pad_zero <- numform::f_pad_zero

  # 100km grid Square easting letter
  esqls <- as.matrix(rbind(c("A","B","C","D","E","F","G","H"),
                           c("J","K","L","M","N","P","Q","R"),
                           c("S","T","U","V","W","X","Y","Z")))
  colnames(esqls) <- 1:8
  rownames(esqls) <- c(1,2,0)

  sin1 <- pi/(180*3600) #One Second
  a <- 6378137 # Equitorial Radius
  b <- 6356752.31424518 #Polar Radius
  k0 <- .9996 #Scalar Factor Constant
  e1 <- sqrt(a^2-b^2)/b
  e1sq <- e1^2
  c <- a^2/b

  gl <- gsub("[^a-zA-Z]", "", MGRS)
  gll <- str_locate(MGRS, gl)
  gzn <- as.numeric(substr(MGRS, 1, gll[1]-1))
  gze <- "Odd"
  if((gzn %% 2) == 0) {gze <- "Even"}
  gzl <- substr(MGRS, gll[1],gll[1])
  gsel <- substr(MGRS, gll[1]+1,gll[1]+1)
  gsnl <- substr(MGRS, gll[1]+2,gll[1]+2)
  hem <- "S"
  if (grepl(gzl, "NPQRSTUVWX")){hem <- "N"}

  east5 <- as.numeric(substr(MGRS, gll[1]+3,gll[1]+7))
  Easting <- east5 + 100000*which(esqls == gsel,arr.ind=T)[2]
  north5 <- as.numeric(substr(MGRS, gll[1]+8,gll[1]+12))
  alln <- data.frame(hem = NA, gze = gze, gzl = gzl, gsnl = gsnl, ngsn = NA)
  fulln<-match_df(northdes,alln,on=c("gze","gzl","gsnl"))
  Northing <- 100000*fulln[1,5]+north5

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

  as.data.frame(rbind(c("MGRS"=MGRS,"Easting"=Easting,"Northing"=Northing,"Zone"=gzn,"Hemisphere"=hem,"Latitude"=round(lat,2),"Longitude"=round(lng,2))))
}
