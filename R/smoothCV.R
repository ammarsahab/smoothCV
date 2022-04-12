
#' dma.dt
#'
#' Function wrapper for double moving averages.
#' @param trainset A set of univariate time series data. Can be a vector (type double) or a data.table.
#' @param m Number of most recent observations where the moving average will be taken at each point in time.
#' Should be less than or equal to the length of training set divided by 2.
#' @param nahead Number of observations to predict.
#' @return A list containing \code{smoothed}
#' (training set, moving averages, level and trend components, and predicted training values)
#' and \code{forc} (forecasted values to use against a test set).
#' @examples
#' dma.dt(crudenow$close,5,nahead=10)


dma.dt<-function(trainset,m,nahead=0){
  if(typeof(trainset)=="double"){
    trainset<-data.table("Data"=trainset)
  }else{
    setnames(trainset,names(trainset)[1],"Data")}
  smoothed<-trainset[,sma:=SMA(Data,n=m)
  ][,dma:=SMA(sma,n=m)
  ][,`:=`(At=2*sma-dma,
          Bt=2/(m-1)*(sma-dma))
  ][,forc:=shift(At+Bt,type="lag")]

  forecast<-last(trainset[,c("At","Bt")],1)[rep(1,nahead)
  ][,Bt:=Bt*seq(nahead)
  ][,forc:=At+Bt]
  return(list(smoothed,forecast,nahead))
}

#' sma.dt
#'
#' Function wrapper for single moving averages.
#' @param trainset A set of univariate time series data. Can be a vector (type double) or a data.table.
#' @param m Number of most recent observations where the moving average will be taken at each point in time.
#' Should be less than or equal to the length of training set.
#' @param nahead Number of observations to predict.
#' @return A list containing \code{smoothed}
#' (training set, moving averages, and predicted training values)
#' and \code{forc} (forecasted values to use against a test set).
#' @examples
#' sma.dt(crudenow$close,5,nahead=10)


sma.dt<-function(trainset,m,nahead=0){
  if(typeof(trainset)=="double"){
    trainset<-data.table("Data"=trainset)
  }else{
    setnames(trainset,names(trainset)[1],"Data")}
  smoothed<-trainset[,sma:=SMA(Data,n=m)
  ][,forc:=shift(sma,type="lag")]

  forecast<-last(trainset[,c("sma")],1)[rep(1,nahead)
  ]
  setnames(forecast,names(forecast)[1],"forc")
  return(list(smoothed,forecast,nahead))
}

#' esWrapper
#'
#' Function wrapper for exponential smoothing.
#' @param trainset A set of univariate time series data. Can be a vector (type double) or a data.table.
#' @param smoothing An object of class \code{HoltWinters}
#' @param nahead Number of observations to predict.
#' @return A list containing \code{smoothed}
#' (training set, moving averages, and predicted training values)
#' and \code{forc} (forecasted values to use against a test set).
#' @examples
#' esWrapper(crudenow$close, HoltWinters(crudenow$close, alpha=0.5, beta=0.5, gamma=F), nahead=10)


esWrapper<-function(trainset,smoothing,nahead=0){
  smoothed<-data.table("Data"     = trainset,
                       "Smoothed" = smoothing[[1]][1],
                       "Level"    = smoothing[[1]][2],
                       "Trend"    = ifelse(smoothing[[4]],
                                           smoothing[[1]][3],0),
                       "Seasonal" = ifelse(smoothing[[5]],
                                           smoothing[[1]][4],0))

  forecast<-data.table(forc = predict(smoothing,nahead))
  setnames(forecast,names(forecast)[1],"forc")
  return(list(smoothed,forecast,nahead))
}

#' smooth.acc
#'
#' Accuracy measures for results of certain smoothing methods.
#' @param smoothresult An object from the output of \code{dma.dt}, \code{sma.dt}, or \code{esWrapper}
#' @param againstself Whether to test forecast accuracy against the training set or a testing set
#' @param testset  A set of univariate time series data. Can be a vector (type double) or a data.table.
#' @return A data.table containing MSE and MAPE of the inputted smoothing method.
#' @examples
#' smooth.acc(sma.dt(crudenow$close,5,nahead=10), againstself=F, crudetest$close)
#' smooth.acc(es.Wrapper(crudenow$close, HoltWinters(crudenow$close, alpha=0.5, beta=0.5, gamma=F), nahead=10),againstself=T)


smooth.acc<-function(smoothresult,againstself=T,testset){
  if(againstself){
    metrics<-smoothresult[[1]][,.(
      "MSE"=mean((Data-forc)^2,na.rm=T),
      "MAPE"=mean(abs((Data-forc)/forc),na.rm=T)*100
    )]
  }else{
    bound<-smoothresult[[3]]
    metrics<-smoothresult[[2]][
      ,test:=testset[1:bound]][
        ,.("MSE"=mean((forc-test)^2),
           "MAPE"=mean(abs((forc-test)/test))*100
    )]
  }
  return(metrics)
}

#' MA.Grid
#'
#' Grid search for single and double moving averages to find optimal parameters that minimize errors
#' with respect to a certain test set.
#' @param trainset A set of univariate time series data. Can be a vector (type double) or a data.table.
#' @param testset Same as above.
#' @param start Starting value for m (see \code{sma.dt} and \code{dma.dt}).
#' @param end Maximum value for m. Must be less than or equal to the length of the training set for SMA,
#' or less than or equal to the length of the training set divided by two for DMA.
#' @param dist Distance between successive values for m.
#' The vector of m values will be constructed as \code{seq(start,end,dist)}.
#' @param nahead Number of observations to forecast.
#' @param type Type of moving average. Can be \code{"SMA"} or \code{"DMA"}
#' @return A data.table containing tested values of m, with MSE and MAPE for each value.
#' @examples
#' MA.Grid(crudenow$Close, crudetest$Close, start=1, end=30, dist=2, nahead=12, type="SMA")
#' MA.Grid(crudenow$Close, crudetest$Close, start=1, end=30, dist=2, nahead=12, type="DMA")

MA.Grid<-function(trainset,testset,start=2,end,dist=1,nahead,type){
    Ms=seq(start,end,dist)
    if(type=="SMA"){
    params<-data.table(M=Ms,
                       MSE=numeric(length(Ms)),
                       MAPE=numeric(length(Ms))
    )

    for(i in seq(1:(length(Ms)))){
      set(params,i,
          c("MSE","MAPE"),
          smooth.acc(smoothresult = sma.dt(trainset,start+dist*(i-1),nahead=nahead),
                     againstself  = F,
                     testset))

    }}else if (type=="DMA"){
      params<-data.table(M=Ms,
                         MSE=numeric(length(Ms)),
                         MAPE=numeric(length(Ms))
      )
      for(i in seq(1:(length(Ms)))){
        set(params,i,
            c("MSE","MAPE"),
            smooth.acc(smoothresult  = dma.dt(trainset,start+dist*(i-1),nahead=nahead),
                       againstself   = F,
                       testset))

      }
    }
return(params)
}

#' ES.Grid
#'
#' Grid search for exponential smoothing to find optimal parameters that minimizes errors
#' with respect to a certain test set.
#' @param type Type of moving average. Can be \code{"SES"} or \code{"DES"}
#' @param nahead Number of observations to forecast.
#' @param alphrange A vector of parameters for the level component. Can be created using \code{seq}, or by manually specifying a vector.
#' @param betarange A vector of parameters for the trend component. Can be created using \code{seq}, or by manually specifying a vector.
#' Not used if type="SES"
#' @param trainset A set of univariate time series data. Can be a vector (type double) or a data.table.
#' @param testset Same as above.
#' @return A data.table containing all combinations of alpha and beta with their respective MSE and MAPE values.
#' @examples
#' ES.Grid(type="SES", alphrange=seq(0.01,1,0.01), betarange=NA, nahead=12, crudenow$Close, crudetest$Close)
#' ES.Grid(type="DES", alphrange=seq(0.1,1,0.1), betarange=seq(0.1,1,0.1), nahead=12, crudenow$Close, crudetest$Close)

ES.grid<-function(type,alphrange,betarange,nahead,trainset,testset){
  if(type=="SES"){
    params<-data.table(alpha=alphrange)
    params[,`:=`(MSE=numeric(length(alphrange)),
                MAPE=numeric(length(alphrange)))]

    for(i in seq(1:length(alphrange))){
      set(params,i,
          c("MSE","MAPE"),
          smooth.acc(esWrapper(smoothing = HoltWinters(trainset,
                                                       alpha=params[[1]][i],
                                                       beta=F, gamma=F),
                               trainset=trainset,
                               nahead=nahead),
                     againstself=F, testset))
    }
  }else if (type=="DES"){
    params<-CJ(alphrange,betarange)
    params[,`:=`(MSE=numeric(nrow(params)),
                MAPE=numeric(nrow(params))
                )]

    for(i in seq(1:(length(alphrange)*length(betarange)))){
      set(params,i,
          c("MSE","MAPE"),
          smooth.acc(esWrapper(smoothing = HoltWinters(trainset,
                                                       alpha=params[[1]][i],
                                                       beta=params[[2]][i],
                                                       gamma=F),
                               trainset=trainset,
                               nahead=nahead),
                     againstself=F, testset))
    }

  }
return(params)
}

#' fcCV
#'
#' Forward chaining cross validation for SMA, DMA, SES, and DES.
#' The algorithm starts by picking a number of observations as the initial training set.
#' The remaining observations will be split into k folds. The first fold will be used as the initial test set.
#' Grid search is then committed using MA.Grid or ES.Grid.
#' The first fold is the incorporated to the training set. The second fold will be used as the next test set.
#' The algorithm repeats itself k times.
#' @param fullset  A set of univariate time series data. Can be a vector (type double) or a data.table.
#' @param initialn Number of observations in the initial training set.
#' @param folds Number of folds to divide the rest of the data into.
#' @param type Type of smoothing. Can be \code{"SMA"}, \code{"DMA"}, \code{"SES"}, \code{"DES"}.
#' @param start (if type="SMA" or type="DMA") Starting value for m (see \code{sma.dt} and \code{dma.dt}).
#' @param end (if type="SMA" or type="DMA") Maximum value for m. Must be less than or equal to the length of the training set for SMA,
#' or less than or equal to the length of the training set divided by two for DMA.
#' @param dist (if type="SMA" or type="DMA") Distance between successive values for m.
#' The vector of m values will be constructed as \code{seq(start,end,dist)}.
#' @param alphrange (if type="SES" or type="DES") A vector of parameters for the level component. Can be created using \code{seq}, or by manually specifying a vector.
#' @param betarange (if type="SES" or type="DES") A vector of parameters for the trend component. Can be created using \code{seq}, or by manually specifying a vector.
#' Not used if type="SES"
#' @return A data.table containing all parameter combinations during each iteration (fold) with their respective MSE and MAPE values.
#' Can be grouped by parameter values to obtain the mean error for each parameter value.
#' @examples
#' fcCV(fullset=crudenow$Close, initialn=36, folds=12, type="DES",
#' alphrange=seq(0.1,1,0.1), betarange=seq(0.1,1,0.1))
#' fcCV(fullset=crudenow$Close, initialn=36, folds=12, type="SMA",
#' start=2, end=30, dist=3)

fcCV<-function(fullset, initialn, folds, type,
               start, end, dist, alphrange, betarange){
  trainset<-fullset[1:initialn]

  if(typeof(fullset)=="double"){
    setlength<-length(fullset)
    lenfold<-(setlength-initialn)/folds

    fullset<-data.table("Data"=fullset)
  }else{
    setlength<-nrow(fullset)
    lenfold<-(setlength-initialn)/folds

    setnames(fullset,names(fullset)[1],"Data")}

  teststart<-initialn+1
  testend<-initialn+lenfold

  testset<-fullset[(teststart):(testend)]
  CVres <- vector(mode='list', length=folds)
  spoints <- rep(NA,folds)
  epoints <-rep(NA,folds)

  if(type=="SMA"){
    for(i in seq(1:folds)){
      spoints[i]<-teststart
      epoints[i]<-testend
      CVres[[i]]<-MA.Grid(trainset,testset,start,end,dist,nrow(testset),type="SMA")

      teststart<-teststart+lenfold
      if((testend+lenfold)<=setlength){
        testend<-testend+lenfold
      }else{
        testend<-setlength
      }
      trainset<-fullset[1:(teststart-1)]
      testset<-fullset[(teststart):(testend)]
    }
    npiter<-length(seq(start,end,dist))
  }else if(type=="DMA"){
    for(i in seq(1:folds)){
      spoints[i]<-teststart
      epoints[i]<-testend
      CVres[[i]]<-MA.Grid(trainset,testset,start,end,dist,nrow(testset),type="DMA")

      teststart<-teststart+lenfold
      if((testend+lenfold)<=setlength){
        testend<-testend+lenfold
      }else{
        testend<-setlength
      }
      trainset<-fullset[1:(teststart-1)]
      testset<-fullset[(teststart):(testend)]
    }
    npiter<-length(seq(start,end,dist))
  }else if(type=="SES"){
      for(i in seq(1:folds)){
        spoints[i]<-teststart
        epoints[i]<-testend
        CVres[[i]]<-ES.grid(type="SES",
                            alphrange=alphrange,
                            betarange=NA,
                            nrow(testset),
                            trainset,testset)

        teststart<-teststart+lenfold
        if((testend+lenfold)<=setlength){
          testend<-testend+lenfold
        }else{
          testend<-setlength
        }
        trainset<-fullset[1:(teststart-1)]
        testset<-fullset[(teststart):(testend)]
      }
    npiter<-length(alphrange)
    }else{
      for(i in seq(1:folds)){
        spoints[i]<-teststart
        epoints[i]<-testend
        CVres[[i]]<-ES.grid(type="DES",
                            alphrange=alphrange,
                            betarange=betarange,
                            nrow(testset),
                            trainset,testset)

        teststart<-teststart+lenfold
        if((testend+lenfold)<=setlength){
          testend<-testend+lenfold
        }else{
          testend<-setlength
        }
        trainset<-fullset[1:(teststart-1)]
        testset<-fullset[(teststart):(testend)]
      }
      npiter<-length(alphrange)*length(betarange)
    }
  CVres<-rbindlist(CVres)
  CVres[,iter:=rep(seq(1:folds),each=npiter)]
  return(list(spoints,epoints,CVres))
}
