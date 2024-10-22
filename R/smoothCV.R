
#' Function wrapper for double moving averages.
#'
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
  if(typeof(trainset)=="double" | typeof(trainset)=="integer"){
    trainset<-data.table("Data"=trainset)
  }else{
    setnames(trainset,names(trainset)[1],"Data")}
  if(!is.data.table(trainset)){
    trainset<-as.data.table(trainset)
  }
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

#' Function wrapper for single moving averages.
#'
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
  if(typeof(trainset)=="double" | typeof(trainset)=="integer"){
    trainset<-data.table("Data"=trainset)
  }else{
    setnames(trainset,names(trainset)[1],"Data")}
  if(!is.data.table(trainset)){
    trainset<-as.data.table(trainset)
  }
  smoothed<-trainset[,sma:=SMA(Data,n=m)
  ][,forc:=shift(sma,type="lag")]

  forecast<-last(trainset[,c("sma")],1)[rep(1,nahead)
  ]
  setnames(forecast,names(forecast)[1],"forc")
  return(list(smoothed,forecast,nahead))
}

#' Function wrapper for exponential smoothing.
#'
#' @param trainset A set of univariate time series data. Can be a vector (type double) or a data.table.
#' @param smoothing An object of class \code{HoltWinters}
#' @param nahead Number of observations to predict.
#' @return A list containing \code{smoothed}
#' (training set, moving averages, and predicted training values)
#' and \code{forc} (forecasted values to use against a test set).
#' @examples
#' esWrapper(crudenow$close, HoltWinters(crudenow$close, alpha=0.5, beta=0.5, gamma=F), nahead=10)


esWrapper<-function(trainset,smoothing,nahead=0){
  lagfactor<-ceiling(smoothing[[3]])+ceiling(smoothing[[4]])+ceiling(smoothing[[5]])
  smoothed<-data.table("Data"     = trainset,
                       "Smoothed" = c(rep(NA,lagfactor),smoothing[[1]][,1]),
                       "Level"    = c(rep(NA,lagfactor),smoothing[[1]][,2]),
                       "Trend"    = ifelse(smoothing[[4]],
                                           c(rep(NA,lagfactor),smoothing[[1]][,3]),
                                           0),
                       "Seasonal" = ifelse(smoothing[[5]],
                                           c(rep(NA,lagfactor),smoothing[[1]][,3]),
                                           0))

  forecast<-data.table(forc = predict(smoothing,nahead))
  setnames(forecast,names(forecast)[1],"forc")
  return(list(smoothed,forecast,nahead))
}

#' Accuracy measures for smoothing results.
#'
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
      "MAPE"=mean(abs((Data-forc)/forc),na.rm=T)*100,
      "MAE"=mean(abs(Data-forc))
    )]
  }else{
    bound<-smoothresult[[3]]
    metrics<-smoothresult[[2]][
      ,test:=testset[1:bound]][
        ,.("MSE"=mean((forc-test)^2),
           "MAPE"=mean(abs((forc-test)/test))*100,
           "MAE"=mean(abs(forc-test))
    )]
  }
  return(metrics)
}

#' Grid search for single and double moving average
#'
#' Grid search for single and double moving averages to find optimal window value that minimize errors
#' with respect to a certain test set. Iterates across given range of window values.
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
                       MAPE=numeric(length(Ms)),
                       MAE=numeric(length(Ms))
    )

    for(i in seq(1:(length(Ms)))){
      set(params,i,
          c("MSE","MAPE","MAE"),
          smooth.acc(smoothresult = sma.dt(trainset,start+dist*(i-1),nahead=nahead),
                     againstself  = F,
                     testset))

    }}else if (type=="DMA"){
      params<-data.table(M=Ms,
                         MSE=numeric(length(Ms)),
                         MAPE=numeric(length(Ms)),
                         MAE=numeric(length(Ms))
      )
      for(i in seq(1:(length(Ms)))){
        set(params,i,
            c("MSE","MAPE","MAE"),
            smooth.acc(smoothresult  = dma.dt(trainset,start+dist*(i-1),nahead=nahead),
                       againstself   = F,
                       testset))

      }
    }
return(params)
}

#' Grid search for exponential smoothing
#'
#' Grid search for exponential smoothing to find optimal parameters that minimizes errors
#' with respect to a certain test set. Only implemented for non-seasonal exponential smoothing
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
                MAPE=numeric(length(alphrange)),
                MAE=numeric(length(alphrange))
    )]

    for(i in seq(1:length(alphrange))){
      set(params,i,
          c("MSE","MAPE","MAE"),
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
                MAPE=numeric(nrow(params)),
                MAE=numeric(nrow(params))
                )]

    for(i in seq(1:(length(alphrange)*length(betarange)))){
      set(params,i,
          c("MSE","MAPE","MAE"),
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

#' Forward chaining cross validation
#'
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
#' @param localoptim Whether to optimize in each training set. If true, a range for paramter values will not need to be provided
#' @param alphrange (if type="SES" or type="DES") A vector of parameters for the level component. Can be created using \code{seq}, or by manually specifying a vector.
#' @param betarange (if type="DES") A vector of parameters for the trend component. Can be created using \code{seq}, or by manually specifying a vector.
#' Not used if type="SES"
#' @return A data.table containing all parameter combinations during each iteration (fold) with their respective MSE and MAPE values.
#' Can be grouped by parameter values to obtain the mean error for each parameter value.
#' @examples
#' fcCV(fullset=crudenow$Close, initialn=36, folds=12, type="DES", localoptim=F,
#' alphrange=seq(0.1,1,0.1), betarange=seq(0.1,1,0.1))
#' fcCV(fullset=crudenow$Close, initialn=36, folds=12, type="SMA", localoptim=F,
#' start=2, end=30, dist=3)

fcCV<-function(fullset, initialn, folds, type,
               start, end, dist, localoptim=F,
               alphrange, betarange){
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

    if(localoptim==F){

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
    }else if(localoptim==T){

      for(i in seq(1:folds)){
        spoints[i]<-teststart
        epoints[i]<-testend
        smoothing<- HoltWinters(trainset,beta=F,gamma=F)
        CVres[[i]]<-data.table(alpha=smoothing[[3]],
                               MSE=numeric(1),
                                MAPE=numeric(1),
                                MAE=numeric(1))
        set(CVres[[i]],1L,c("MSE","MAPE","MAE"),
            smooth.acc(esWrapper(smoothing = smoothing,trainset=trainset,
                                 nahead= nrow(testset)),
                       againstself=F, testset))

        teststart<-teststart+lenfold
        if((testend+lenfold)<=setlength){
          testend<-testend+lenfold
        }else{
          testend<-setlength
        }
        trainset<-fullset[1:(teststart-1)]
        testset<-fullset[(teststart):(testend)]
      }
      npiter<-1
    }
  }else{

      if(localoptim==F){

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
    }else if(localoptim==T){

        for(i in seq(1:folds)){
        spoints[i]<-teststart
        epoints[i]<-testend
        smoothing<- HoltWinters(trainset,gamma=F)
        CVres[[i]]<-data.table(alpha=smoothing[[3]],
                               beta=smoothing[[4]],
                               MSE=numeric(1),
                               MAPE=numeric(1),
                               MAE=numeric(1))
        set(CVres[[i]],1L,c("MSE","MAPE","MAE"),
            smooth.acc(esWrapper(smoothing = smoothing,trainset=trainset,
                                 nahead=nrow(testset)),
                       againstself=F, testset))

        teststart<-teststart+lenfold
        if((testend+lenfold)<=setlength){
          testend<-testend+lenfold
        }else{
          testend<-setlength
        }
        trainset<-fullset[1:(teststart-1)]
        testset<-fullset[(teststart):(testend)]
      }
      npiter<-1
    }
  }
  CVres<-rbindlist(CVres)
  CVres[,iter:=rep(seq(1:folds),each=npiter)]
  return(list(spoints,epoints,CVres))
}


#' Rolling block cross validation
#'
#' The algorithm starts by splitting all observations into k folds.
#' Specified portions of the fold will be split into training and test sets.
#' Grid search is then committed using MA.Grid or ES.Grid.
#' This occurs in all k folds.
#' @param fullset  A set of univariate time series data. Can be a vector (type double) or a data.table.
#' @param folds Number of folds to divide the rest of the data into.
#' @param trainsplit Proportion of observations to be used as training set in each fold.
#' The remaining observations will automatically be allocated to the test set.
#' @param type Type of smoothing. Can be \code{"SMA"}, \code{"DMA"}, \code{"SES"}, \code{"DES"}.
#' @param start (if type="SMA" or type="DMA") Starting value for m (see \code{sma.dt} and \code{dma.dt}).
#' @param end (if type="SMA" or type="DMA") Maximum value for m.
#' Must be less than or equal to the length of the training set for SMA,
#' or less than or equal to the length of the training set divided by two for DMA.
#' @param dist (if type="SMA" or type="DMA") Distance between successive values for m.
#' The vector of m values will be constructed as \code{seq(start,end,dist)}.
#' @param localoptim Whether to optimize in each training set.
#' If true, there is no need to provide values for alpha and beta..
#' @param alphrange (if type="SES" or type="DES") A vector of parameters for the level component.
#' Can be created using \code{seq}, or by manually specifying a vector.
#' @param betarange (if type="DES") A vector of parameters for the trend component.
#' Can be created using \code{seq}, or by manually specifying a vector.
#' Not used if type="SES"
#' @return A data.table containing all parameter combinations during each iteration (fold) with their respective MSE and MAPE values.
#' Can be grouped by parameter values to obtain the mean error for each parameter value.
#' @examples
#' fcCV(fullset=crudenow$Close, initialn=36, folds=12, type="DES", localoptim=F,
#' alphrange=seq(0.1,1,0.1), betarange=seq(0.1,1,0.1))
#' fcCV(fullset=crudenow$Close, initialn=36, folds=12, type="SMA", localoptim=F,
#' start=2, end=30, dist=3)

rbCV<-function(fullset, trainsplit, folds, type,
               start, end, dist, localoptim=F,
               alphrange, betarange){
  if(typeof(fullset)=="double"){
    setlength<-length(fullset)

    fullset<-data.table("Data"=fullset)
  }else{
    setlength<-nrow(fullset)

    setnames(fullset,names(fullset)[1],"Data")}

  foldsize<-ceiling(nrow(fullset)/folds)

  foldstart<-1
  foldend<-1+foldsize

  CVres <- vector(mode='list', length = folds)
  foldloc <- vector(mode='list', length = folds)
  trainloc <-vector(mode='list', length = folds)
  testloc <- vector(mode='list', length = folds)

  if(type=="SMA"){

    for(i in seq(1:folds)){
      trainend<- foldstart+ceiling(trainsplit*(foldend-foldstart))
      trainset<-fullset[foldstart:trainend]
      testset<- fullset[(trainend+1):foldend]

      foldloc[[i]]<-c(foldstart,foldend)
      trainloc[[i]]<-c(foldstart,trainend)
      testloc[[i]]<-c((trainend+1),foldend)
      CVres[[i]]<-MA.Grid(trainset,testset,start,end,dist,nrow(testset),type="SMA")

      foldstart<-foldend+1
      if((foldstart+foldsize)<=setlength){
        foldend<-foldstart+foldsize

      }else{
        foldend<-setlength

      }
    }
    npiter<-length(seq(start,end,dist))
  }else if(type=="DMA"){

    for(i in seq(1:folds)){
      trainend<- foldstart+ceiling(trainsplit*(foldend-foldstart))
      trainset<-fullset[foldstart:trainend]
      testset<- fullset[(trainend+1):foldend]

      foldloc[[i]]<-c(foldstart,foldend)
      trainloc[[i]]<-c(foldstart,trainend)
      testloc[[i]]<-c((trainend+1),foldend)
      CVres[[i]]<-MA.Grid(trainset,testset,start,end,dist,nrow(testset),type="DMA")

      foldstart<-foldend+1
      if((foldstart+foldsize)<=setlength){
        foldend<-foldstart+foldsize

      }else{
        foldend<-setlength

      }
    }
    npiter<-length(seq(start,end,dist))
  }else if(type=="SES"){

    if(localoptim==F){

      for(i in seq(1:folds)){
        trainend<- foldstart+ceiling(trainsplit*(foldend-foldstart))
        trainset<-fullset[foldstart:trainend]
        testset<- fullset[(trainend+1):foldend]

        foldloc[[i]]<-c(foldstart,foldend)
        trainloc[[i]]<-c(foldstart,trainend)
        testloc[[i]]<-c((trainend+1),foldend)

        CVres[[i]]<-ES.grid(type="SES",
                            alphrange=alphrange,
                            betarange=NA,
                            nrow(testset),
                            trainset,testset)


        foldstart<-foldend+1
        if((foldstart+foldsize)<=setlength){
          foldend<-foldstart+foldsize

        }else{
          foldend<-setlength

        }
      }
      npiter<-length(alphrange)
    }else if(localoptim==T){

      for(i in seq(1:folds)){
        trainend<- foldstart+ceiling(trainsplit*(foldend-foldstart))
        trainset<-fullset[foldstart:trainend]
        testset<- fullset[(trainend+1):foldend]

        foldloc[[i]]<-c(foldstart,foldend)
        trainloc[[i]]<-c(foldstart,trainend)
        testloc[[i]]<-c((trainend+1),foldend)

        smoothing<- HoltWinters(trainset,beta=F,gamma=F)
        CVres[[i]]<-data.table(alpha=smoothing[[3]],
                               MSE=numeric(1),
                               MAPE=numeric(1),
                               MAE=numeric(1))

        set(CVres[[i]],1L,c("MSE","MAPE","MAE"),
            smooth.acc(esWrapper(smoothing = smoothing,trainset=trainset,
                                 nahead= nrow(testset)),
                       againstself=F, testset))

        foldstart<-foldend+1
        if((foldstart+foldsize)<=setlength){
          foldend<-foldstart+foldsize

        }else{
          foldend<-setlength

        }
      }
      npiter<-1
    }
  }else{

    if(localoptim==F){

      for(i in seq(1:folds)){
        trainend<- foldstart+ceiling(trainsplit*(foldend-foldstart))
        trainset<-fullset[foldstart:trainend]
        testset<- fullset[(trainend+1):foldend]

        foldloc[[i]]<-c(foldstart,foldend)
        trainloc[[i]]<-c(foldstart,trainend)
        testloc[[i]]<-c((trainend+1),foldend)

        CVres[[i]]<-ES.grid(type="DES",
                            alphrange=alphrange,
                            betarange=betarange,
                            nrow(testset),
                            trainset,testset)

        foldstart<-foldend+1
        if((foldstart+foldsize)<=setlength){
          foldend<-foldstart+foldsize

        }else{
          foldend<-setlength

        }
      }
      npiter<-length(alphrange)*length(betarange)
    }else if(localoptim==T){

      for(i in seq(1:folds)){
        trainend<- foldstart+ceiling(trainsplit*(foldend-foldstart))
        trainset<-fullset[foldstart:trainend]
        testset<- fullset[(trainend+1):foldend]

        foldloc[[i]]<-c(foldstart,foldend)
        trainloc[[i]]<-c(foldstart,trainend)
        testloc[[i]]<-c((trainend+1),foldend)

        smoothing<- HoltWinters(trainset,gamma=F)
        CVres[[i]]<-data.table(alpha=smoothing[[3]],
                               beta=smoothing[[4]],
                               MSE=numeric(1),
                               MAPE=numeric(1),
                               MAE=numeric(1))

        set(CVres[[i]],1L,c("MSE","MAPE","MAE"),
            smooth.acc(esWrapper(smoothing = smoothing,trainset=trainset,
                                 nahead=nrow(testset)),
                       againstself=F, testset))

        foldstart<-foldend+1
        if((foldstart+foldsize)<=setlength){
          foldend<-foldstart+foldsize

        }else{
          foldend<-setlength

        }
      }
      npiter<-1
    }
  }
  CVres<-rbindlist(CVres)
  CVres[,iter:=rep(seq(1:folds),each=npiter)]
  return(list(foldloc,trainloc,testloc,CVres))
}

#' Repeated holdout cross-validation
#'
#' The method splits the data into a training set, a testing set, and a set which will be split into
#' training and testing portions based on a random number generator.
#' Grid search, or optimization, is then committed for k iterations
#' @param fullset  A set of univariate time series data. Can be a vector (type double) or a data.table.
#' @param iter Number of repetitions
#' @param trainsplit Proportion of observations to be used as training set
#' @param randsplit Proportion of observations in the set with uncertain membership.
#' Will be split into training and testing proportions.
#' @param type Type of smoothing. Can be \code{"SMA"}, \code{"DMA"}, \code{"SES"}, \code{"DES"}.
#' @param start (if type="SMA" or type="DMA") Starting value for m (see \code{sma.dt} and \code{dma.dt}).
#' @param end (if type="SMA" or type="DMA") Maximum value for m.
#' Must be less than or equal to the length of the training set for SMA,
#' or less than or equal to the length of the training set divided by two for DMA.
#' @param dist (if type="SMA" or type="DMA") Distance between successive values for m.
#' The vector of m values will be constructed as \code{seq(start,end,dist)}.
#' @param localoptim Whether to optimize in each training set.
#' If true, there is no need to provide values for alpha and beta..
#' @param alphrange (if type="SES" or type="DES") A vector of parameters for the level component.
#' Can be created using \code{seq}, or by manually specifying a vector.
#' @param betarange (if type="DES") A vector of parameters for the trend component.
#' Can be created using \code{seq}, or by manually specifying a vector.
#' Not used if type="SES"
#' @return A data.table containing all parameter combinations during each iteration (fold) with their respective MSE and MAPE values.
#' Can be grouped by parameter values to obtain the mean error for each parameter value.
#' @examples
#' fcCV(fullset=crudenow$Close, initialn=36, folds=12, type="DES", localoptim=F,
#' alphrange=seq(0.1,1,0.1), betarange=seq(0.1,1,0.1))
#' fcCV(fullset=crudenow$Close, initialn=36, folds=12, type="SMA", localoptim=F,
#' start=2, end=30, dist=3)

rhCV<-function(fullset, trainsplit, randsplit, iter, type,
               start, end, dist, localoptim=F,
               alphrange, betarange){

  if(typeof(fullset)=="double"){
    setlength<-length(fullset)

    fullset<-data.table("Data"=fullset)
  }else{
    setlength<-nrow(fullset)

    setnames(fullset,names(fullset)[1],"Data")}
  CVres <- vector(mode='list', length = folds)
  foldloc <- vector(mode='list', length = folds)

  randstart<-ceiling(trainsplit*setlength)
  randend<-ceiling((trainsplit+randsplit)*setlength)

  if(type=="SMA"){

    for(i in seq(1:iter)){
      trainend<-round(runif(1,min=randstart, max=randend))
      trainset<-fullset[1:trainend]
      testset<- fullset[(trainend+1):setlength]

      foldloc[[i]]<-c(1,trainend,trainend+1,setlength)
      CVres[[i]]<-MA.Grid(trainset,testset,start,end,dist,nrow(testset),type="SMA")
      }
    npiter<-length(seq(start,end,dist))
  }else if(type=="DMA"){

    for(i in seq(1:iter)){
      trainend<-round(runif(1,min=randstart, max=randend))
      trainset<-fullset[1:trainend]
      testset<- fullset[(trainend+1):setlength]

      foldloc[[i]]<-c(1,trainend,trainend+1,setlength)

      CVres[[i]]<-MA.Grid(trainset,testset,start,end,dist,nrow(testset),type="DMA")
    }
    npiter<-length(seq(start,end,dist))

  }else if(type=="SES"){

    if(localoptim==F){

      for(i in seq(1:iter)){
        trainend<-round(runif(1,min=randstart, max=randend))
        trainset<-fullset[1:trainend]
        testset<- fullset[(trainend+1):setlength]

        foldloc[[i]]<-c(1,trainend,trainend+1,setlength)

        CVres[[i]]<-ES.grid(type="SES",
                            alphrange=alphrange,
                            betarange=NA,
                            nrow(testset),
                            trainset,testset)
      }
      npiter<-length(alphrange)
    }else if(localoptim==T){

      for(i in seq(1:iter)){
        trainend<-round(runif(1,min=randstart, max=randend))
        trainset<-fullset[1:trainend]
        testset<- fullset[(trainend+1):setlength]

        foldloc[[i]]<-c(1,trainend,trainend+1,setlength)

        smoothing<- HoltWinters(trainset,beta=F,gamma=F)
        CVres[[i]]<-data.table(alpha=smoothing[[3]],
                               MSE=numeric(1),
                               MAPE=numeric(1),
                               MAE=numeric(1))

        set(CVres[[i]],1L,c("MSE","MAPE","MAE"),
            smooth.acc(esWrapper(smoothing = smoothing,trainset=trainset,
                                 nahead= nrow(testset)),
                       againstself=F, testset))

      }
      npiter<-1
    }
  }else{

    if(localoptim==F){

      for(i in seq(1:iter)){
        trainend<-round(runif(1,min=randstart, max=randend))
        trainset<-fullset[1:trainend]
        testset<- fullset[(trainend+1):setlength]

        foldloc[[i]]<-c(1,trainend,trainend+1,setlength)

        CVres[[i]]<-ES.grid(type="DES",
                            alphrange=alphrange,
                            betarange=betarange,
                            nrow(testset),
                            trainset,testset)

      }
      npiter<-length(alphrange)*length(betarange)
    }else if(localoptim==T){

      for(i in seq(1:iter)){
        trainend<-round(runif(1,min=randstart, max=randend))
        trainset<-fullset[1:trainend]
        testset<- fullset[(trainend+1):setlength]

        foldloc[[i]]<-c(1,trainend,trainend+1,setlength)

        smoothing<- HoltWinters(trainset,gamma=F)
        CVres[[i]]<-data.table(alpha=smoothing[[3]],
                               beta=smoothing[[4]],
                               MSE=numeric(1),
                               MAPE=numeric(1),
                               MAE=numeric(1))

        set(CVres[[i]],1L,c("MSE","MAPE","MAE"),
            smooth.acc(esWrapper(smoothing = smoothing,trainset=trainset,
                                 nahead=nrow(testset)),
                       againstself=F, testset))

      }
      npiter<-1
    }
  }
  CVres<-rbindlist(CVres)
  CVres[,iter:=rep(seq(1:iter),each=npiter)]
  return(list(foldloc,CVres))
}
