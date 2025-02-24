#' Simulate Ratkowsky Model
#'
#' @param pars = vector of parameters (d,Tl,g,Tu)
#' @param x = variable to define parameters over
#'
#' @return Response Variable
#' @export
#' @references Ratkowsky, D.A., Lowry, R.K., McMeekin, T.A., Stokes, A.N., & Chandler, R.A. (1983). Model for bacterial culture growth rate through out the entire biokinetic temperature range. Journal of Bacteriology 154: 1222-1226.
#' @examples
#' theta <- c(d=0.1,Tl=2,g=0.2,Tu=20)
#' x=rep(2:20,each=1)
#' ratkowsky(theta,x)


ratkowsky = function(pars,x) {
  pars[1]*(x-pars[2])*(1-exp(pars[3]*(x-pars[4]))) #these function needs to be constrainted; anything beyond Tmin and Tmax is negative
  }
