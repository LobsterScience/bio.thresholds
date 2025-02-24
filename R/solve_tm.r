#' Finds Optimum (Tm) from Ratkowsky Model fits
#'
#' @param fit from NLS
#'
#' @return Estimate of Optimum (or maximum)
#' @export
#' @references Jonsson, B., Forseth, T., Jensen, A.J., and Naesje, T.F. 2001. Thermal performance of juvenile Atlantic salmon, Salmo salar L. Functional Ecology 15: 701-711.
#' @examples solve_tm(fit)
 solve_tm<-function(fit){
  g=coef(fit)[2]
  Tl=coef(fit)[3]
  Tu=coef(fit)[4]
  Tm=seq(Tl,Tu,0.01)
  X1<-log(1+g*(Tm-Tl))+g*(Tm-Tu) #equation 5
  Topt<-Tm[which(round(X1*100)/100==0)]
  if(length(Topt)>1) {
    Topt<-median(Topt)
  }
  return(Topt)
}
