### Other Functions


CalcNicheFitDiffs_fromSTANfit<-function(A_draws, R_draws){
  
  if(!is.matrix(A_draws)){
    stop('Expects matrix. NB, changed function from taking a STAN fit object to allow pre-filtering')
    
  }
  
  r <-R_draws
  
  ### Arranging 
  A <- array(t(A_draws), dim = c(6,6, nrow(A_draws)))
  
  
  R<- r-1 # Shift to intrinsic growth rate from reproductive value Rows = draws, cols = species
  
  n = nrow(r)
  species <- c('BIR','PAL', 'PAN', 'PSA', 'SIM', 'SUL')
  
  
  # #################################
  # ### Checking this is right
  #  array(7:42, dim = c(6,6))  ## Old version need to replicate
  # Orig <- matrix(1:(36*4), ncol = 36, byrow=T)
  # Atmp <- array(  t(Orig), dim = c(6,6, 4))
  # Atmp[,,1]
  # ##################
  
  NicheOverlap = array(NA, c(6,6,n  ))
  FitnessDifference = array(NA, c(6,6,n  ))
  
  ### Vectorised over draws (3rd dimension)
  for(i in 1:6){
    for(j in 1:6){
      if(i ==j){ 
        NicheOverlap[i,j,]<- NA
        FitnessDifference[i,j,]<- NA
      }else{
        NicheOverlap[i,j,]<- sqrt(  (A[i,j,]/A[j,j,]) * (A[j,i,]/A[i,i,])   )
        FitnessDifference[i,j,] <-   (R[,i]/R[,j])   *  sqrt(  (A[j,j,]*A[j,i,]) / (A[i,i,]*A[i,j,])) 
      }
    }
  }
  
  NicheDiff = 1-NicheOverlap 
  # NicheDiff[NicheDiff<0]<-0  ##### Possibly remove this? or make NA
  
  colnames(NicheDiff)<- species
  rownames(NicheDiff)<- species
  
  colnames(FitnessDifference)<- species
  rownames(FitnessDifference)<- species
  
  return(list('FitDiff'=FitnessDifference, 'NicheDiff'=NicheDiff))
} 


ExtractDiffs<-function(i,  Ni_FiDraws ){
  FitnessDifference<-  Ni_FiDraws$FitDiff
  NicheDiff        <-  Ni_FiDraws$NicheDiff
  ToSelect<-which(lower.tri(FitnessDifference[,,1]), arr.ind = TRUE)
  
  df <- data.frame(FitnessDiff = FitnessDifference[ToSelect[i,1],ToSelect[i,2],] ,
                   NicheDiff =         NicheDiff[ToSelect[i,1],ToSelect[i,2],])
  
  sp<- c('BIR','PAL', 'PAN', 'PSA', 'SIM', 'SUL')
  df$Combination <- paste0(sp[ToSelect[i,1]],'_',sp[ToSelect[i,2]])
  
  return(df)
}



PostMean_NicheFitDiffs_fromSTANfit<-function(fit_object){
  species <- c('BIR','PAL', 'PAN', 'PSA', 'SIM', 'SUL')
  
  FIT1_Draws<-  rstan::extract(fit_object )
  x<- summary( fit_object)$summary[,1]
  r <-x[1:6]
  
  R<- r-1 # Shift to intrinsic growth rate from reproductive value
  
  A <- array(t(x[7:42]), dim = c(6,6))
  
  NicheOverlap = array(NA, c(6,6))
  FitnessDifference = array(NA, c(6,6))
  
  for(i in 1:6){
    for(j in 1:6){
      if(i ==j){ 
        NicheOverlap[i,j]<- NA
        FitnessDifference[i,j]<- NA
      }else{
        NicheOverlap[i,j]<- sqrt(  (A[i,j]/A[j,j]) * (A[j,i]/A[i,i])   )
        FitnessDifference[i,j] <-   (R[i]/R[j])   *  sqrt(  (A[j,j]*A[j,i]) / (A[i,i]*A[i,j])) 
        
      }
    }
  }
  
  NicheDiff = 1-NicheOverlap 
  # NicheDiff[NicheDiff<0]<-0  ##### Possibly remove this? or make NA
  
  colnames(NicheDiff)<- species
  rownames(NicheDiff)<- species
  
  colnames(FitnessDifference)<- species
  rownames(FitnessDifference)<- species
  
  return(list('FitDiff'=FitnessDifference, 'NicheDiff'=NicheDiff))
} 


Calc_FitDiff_components_fromSTANfit<- function(A_draws, R_draws){
  
  r <-R_draws
  A <- array(t(A_draws), dim = c(6,6, nrow(A_draws)))
  R<- r-1 # Shift to intrinsic growth rate from reproductive value Rows = draws, cols = species
  n = nrow(r)
  species <- c('BIR','PAL', 'PAN', 'PSA', 'SIM', 'SUL')
  
  R_Part = array(NA, c(6,6,n  ))
  A_Part = array(NA, c(6,6,n  ))
  
  for(i in 1:6){
    for(j in 1:6){
      if(i ==j){ 
        R_Part[i,j,]<- NA
        A_Part[i,j,]<- NA
      }else{
        R_Part[i,j,]<-   (R[,i]/R[,j]) 
        A_Part[i,j,] <- sqrt(  (A[j,j,]*A[j,i,]) / (A[i,i,]*A[i,j,])) 
      }
    }
  }
  colnames(A_Part)<- species
  rownames(A_Part)<- species
  colnames(R_Part)<- species
  rownames(R_Part)<- species
  
  
  LongFormDF <-  map_df(1:15, function(i){
    ToSelect<-which(lower.tri(A_Part[,,1]), arr.ind = TRUE)
    
    df <- data.frame(A_Part = A_Part[ToSelect[i,1],ToSelect[i,2],] ,
                     R_Part =         R_Part[ToSelect[i,1],ToSelect[i,2],])
    df$Combination <- paste0(species[ToSelect[i,1]],'_',species[ToSelect[i,2]])
    return(df)
  })
  return(LongFormDF)
}
