#' Posterior Means and 95\% C.I.s of the NIE, NDE and TE
#'
#' Obtain posterior means and credible intervals of the effects.
#' @param dataTreatment The observed data under Z=1
#' @param dataControl The observed data under Z=0
#' @param prior a list giving the prior information
#' @param mcmc a list giving the MCMC parameters
#' @param state a list giving the current value of the parameters
#' @param status a logical variable indicating whether this run is new (TRUE) or the continuation of a previous analysis (FALSE)
#' @param na.action a function that indicates what should happen when the data contain NAs
#' @param q A dimension of the observed data, i.e., number of covariates plus 2
#' @param NN Number of samples drawn for each iteration from the joint distribution of the mediator and the covariates. Default is 10
#' @param n1 Number of observations under Z=1
#' @param n0 Number of observations under Z=0
#' @param extra.thin Giving the extra thinning interval
#' @param seed Value to be given to the seed
#'
#' @examples
#' \dontrun{
#' # Data
#' dataTreatment = c(11.163989,  13.357644,  12.264292,  16.848934,
#'  20.902562,  13.037647, 10.523266, 13.869733,  16.399047,  13.249197,
#'   11.284890,  11.450062,  10.817139, 15.245700, 9.920019,  11.573960,
#'    10.911009,  10.489587, 10.814474, 10.479123,  11.576099, 9.527664,
#'     19.139293,  10.419879, 11.477210, 9.974351,  11.661321,  11.107294,
#'      10.703919,  10.229351, 14.304395, 11.966633,  12.320991,  17.905459,
#'       9.661749, 10.385607, 20.406423,  10.777457, 9.825228,  13.691681,
#'        9.871743,  10.500229, 10.285037,  14.930915,  13.608770,  16.167101,
#'        13.602650,  14.786788, 10.867811, 10.512889,  11.412777,  11.982880,
#'         8.774465,  10.481118, 9.444334,  15.367506,  9.772886,  10.974271,
#'          9.893823,  20.102887,  11.444247,  10.500275, 9.847932,  11.732694,
#'           11.462192,  13.380679,  12.794651,  9.678496,  17.062716,  14.091250,
#'            9.289893,  10.719766,  10.636155,  9.798370,  10.298192,  10.618451,
#'            10.365415,  18.189078, 15.972311,  11.072824,  12.343840,  11.811951,
#'             10.639700,  15.114964, 11.880889,  11.335969, 11.388849,  10.878676,
#'              10.722346,  11.307921,  9.956277,  9.213929,  10.085210,  11.800195,
#'               10.658099,  10.718170,  11.422773,  10.623135,  11.142688,  10.719141 )
#' dataControl = c(9.060491, 7.566083, 8.081978, 6.450560, 7.252086,
#' 7.545289, 9.019265, 7.678507, 7.353213, 10.708195, 8.253596, 8.413139,
#'  9.588276, 10.589635 ,8.619169,9.088016, 7.063381, 6.403902,  8.413106,
#'   7.254758, 7.946906, 7.013008, 6.828835, 7.198124, 8.133756, 14.843925,
#'    8.043300, 9.978774, 11.630415, 6.912515, 7.354103, 6.634561, 9.672022,
#'    9.785791, 7.490891, 12.384268, 9.223109, 6.719683, 8.311750, 9.388895,
#'     7.773400, 14.815319, 20.050661, 6.417857, 10.674383, 9.293519, 14.880615,
#'      8.090496, 7.880399, 7.181893, 7.913207, 7.086220, 10.728450, 8.687074,
#'       6.958637, 14.135874, 7.273476, 10.325233, 7.736805, 14.173735, 7.773935,
#'        9.227020, 8.257979 , 11.986135, 7.875595, 7.341740, 6.036009, 8.976403,
#'        7.396193, 7.217195, 6.849920, 10.399304, 12.968451,7.879505, 7.816833,
#'         7.270345, 7.299595, 10.799175, 6.471268, 6.234942, 8.203609, 6.411400,
#'         10.768465, 8.162683, 7.262786, 8.473937, 9.886522, 7.018580, 6.836707,
#'          7.283106, 7.344599, 8.686345, 9.547393, 10.468727, 14.321967, 7.183868,
#'           10.743336, 6.865275, 7.725663, 7.168826)
#'
#' # Initial state
#' state <- NULL

#' MCMC parameters
#' nburn <- 500
#' nsave <- 1000
#' nskip <- 10
#' ndisplay <- 100
#' mcmc <- list(nburn=nburn,nsave=nsave,nskip=nskip,ndisplay=ndisplay)
#'
#' Prior information
#' s2 <- matrix(c(1000,0,0,1),ncol=2)
#' m2 <- c(20,3)
#' psiinv2 <- solve(matrix(c(1000,0,0,1),ncol=2))
#' prior <- list(a0=1,b0=1/5,nu1=4,nu2=4,s2=s2,m2=m2,psiinv2=psiinv2,tau1=0.01,tau2=0.01)

#' library(BNPMediation)
#' result <- bnpmediation(dataTreatment, dataControl, prior, mcmc, state, status = TRUE,
#'  na.action, q = 2, NN = 10, n1 = 10, n0 = 10, extra.thin = 0)
#'
#' # Use the plot_effects function from the package on the final output.
#'  plot_effects(result)
#' }
#'
#'
#' @return ENIE Posterior mean of the Natural Indirect Effect (NIE)
#' @return ENDE Posterior mean of the Natural Direct Effect (NDE)
#' @return ETE Posterior mean of the Total Effect (TE)
#' @return IE.c.i 95\% C.I. of the NIE
#' @return DE.c.i 95\% C.I. of the NDE
#' @return TE.c.i 95\% C.I. of the TE
#' @return Y11 Posterior samples of Y11
#' @return Y00 Posterior samples of Y00
#' @return Y10 Posterior samples of Y10
#'
#' @importFrom mnormt rmnorm dmnorm
#' @importFrom stats dnorm rnorm na.omit
#' @importFrom utils setTxtProgressBar txtProgressBar
#' @import DPpackage


#' @export
bnpmediation<-function(dataTreatment, dataControl, prior, mcmc, state, status=TRUE,na.action, q=2, NN=10, n1=10, n0=10, extra.thin=0, seed=12345)

{
  cat("***** Fitting observed data models via DPpackage::DPdensity()\n")
  obj1 = DPdensity(y=dataTreatment,prior=prior,mcmc=mcmc,state=state,status=TRUE, na.action=na.omit)
  obj0 = DPdensity(y=dataControl,prior=prior,mcmc=mcmc,state=state,status=TRUE, na.action=na.omit)

  cat("***** Running bnpmediation\n")

  obj1.dim <- dim(obj1$save.state$randsave)[2]-(q*(q+1)/2+2*q-1)
  obj0.dim <- dim(obj0$save.state$randsave)[2]-(q*(q+1)/2+2*q-1)

  Len.MCMC <- 1:dim(obj0$save.state$randsave)[1]
  if(extra.thin!=0){
    Len.MCMC <- Len.MCMC[seq(1, length(Len.MCMC), extra.thin)]
  }

  Ysamples<-OutSamples(obj1, obj0, q)
  Y11 <- Ysamples$Y1[Len.MCMC]
  Y00 <- Ysamples$Y0[Len.MCMC]

  set.seed(seed)

  mat.given.ij <- function(x, y) ifelse(x <= y, (q-1)*(x-1)+y-x*(x-1)/2, (q-1)*(y-1)+x-y*(y-1)/2)
  mat <- function(q) outer( 1:q, 1:q, mat.given.ij )

  pb <- txtProgressBar(min = 0, max = length(Len.MCMC), style = 3)

  Y10<-NULL

  index<-0
  for(j in Len.MCMC){
    index <- index + 1
    mu2 <- sapply(seq(2,obj0.dim, by=(q*(q+1)/2+q)), function(x)  obj0$save.state$randsave[j,x[1]:(x[1]+q-2)])
    sigma22 <- sapply(seq(q+q+1,obj0.dim, by=(q*(q+1)/2+q)), function(x)  obj0$save.state$randsave[j,x[1]:(x[1]+(q-1)*(q)/2-1)][mat(q-1)])
    if(q!=2){
      joint0 <- do.call("rbind", replicate(NN, data.frame(sapply(1:n0, function(x) rmnorm(1,mu2[,x],matrix(sigma22[,x],q-1,q-1,byrow=T) )))))
    }else{
      joint0 <- matrix(replicate(NN, sapply(1:n0, function(x) rnorm(1,mu2[x],sd=sqrt(sigma22[x]) )), simplify="array"), nrow=n0*NN)
    }
    unique.val <- unique(obj1$save.state$randsave[j,seq(1,obj1.dim,by=(q*(q+1)/2+q))])
    unique.ind <- NULL
    unique.prop <- NULL
    for(k in 1:length(unique.val)){
      unique.ind[k] <- which(obj1$save.state$randsave[j,seq(1,obj1.dim,by=(q*(q+1)/2+q))]==unique.val[k])[1]
      unique.prop[k] <- length(which(obj1$save.state$randsave[j,seq(1,obj1.dim,by=(q*(q+1)/2+q))]==unique.val[k]))/n1
    }
    b01 <- NULL
    Weight.num0 <- matrix(nrow=length(unique.val), ncol=n0*NN)
    B0 <- matrix(nrow=length(unique.val),ncol=n0*NN)

    t.ind<-0
    for(k in unique.ind){
      t.ind<-1+t.ind
      mu1<-obj1$save.state$randsave[j,(q*(q+1)/2+q)*k-(q*(q+1)/2+q)+1]
      mu2<-obj1$save.state$randsave[j,((q*(q+1)/2+q)*k-(q*(q+1)/2+q)+2):((q*(q+1)/2+q)*k-(q*(q+1)/2+q)+q)]
      sigma1<-obj1$save.state$randsave[j,(q*(q+1)/2+q)*k-(q*(q+1)/2+q)+q+1]
      sigma12<-obj1$save.state$randsave[j,(q*(q+1)/2+q)*k-(q*(q+1)/2+q)+((q+2):(2*q))]
      sigma22<-matrix(obj1$save.state$randsave[j,((q*(q+1)/2+q)*k-(q*(q+1)/2+q)+2*q+1):((q*(q+1)/2+q)*k)][mat(q-1)],q-1,q-1,byrow=TRUE)
      if(q!=2){
        Weight.num0[t.ind,1:(n0*NN)]<-unique.prop[t.ind]*dmnorm(joint0,mu2,sigma22)
      }else{
        Weight.num0[t.ind,1:(n0*NN)]<-unique.prop[t.ind]*dnorm(joint0,mu2,sd=sqrt(sigma22))
      }
      b01[t.ind]<-mu1-sigma12%*%solve(sigma22)%*%t(t(mu2))
      B0[t.ind,1:(n0*NN)]<-sigma12%*%solve(sigma22)%*%t(joint0)
    }
    Weight=apply(Weight.num0, 2, function(x) x/sum(x))
    test <- Weight*(b01+B0)
    Y10[index]<-mean(apply(test, 2, sum))
    Sys.sleep(0.05)
    setTxtProgressBar(pb, index)
  }

  z <- list(Y11=Y11,
            Y00=Y00,
            Y10=Y10,
            ENIE=mean(Y11-Y10),
            ENDE=mean(Y10-Y00),
            ETE=mean(Y11-Y00),
            TE.c.i=c(sort(Y11-Y00)[length(Len.MCMC)*0.025],sort(Y11-Y00)[length(Len.MCMC)*0.975]),
            IE.c.i=c(sort(Y11-Y10)[length(Len.MCMC)*0.025],sort(Y11-Y10)[length(Len.MCMC)*0.975]),
            DE.c.i=c(sort(Y10-Y00)[length(Len.MCMC)*0.025],sort(Y10-Y00)[length(Len.MCMC)*0.975]))
  z$call <- match.call()
  class(z) <- "bnpmediation"
  return(z)
}


