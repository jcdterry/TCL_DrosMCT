
Write_Stan_NegB_2R_2A<-function(d,
                                Name,
                                AlphaLowerBound = 0,
                                AlphaSD= 1, 
                                CalcGQ = TRUE){
  
  # Function to generate a Stan model
  # uses raw count data and a Negative-Binomial error distrubtion
  # Fits separate A's to the two treatments 
  
  #############
  ### Core Start
  #########
  
  Base<- paste0('data {
  int<lower=0> N;
  int y[N];// Observed Count of emergences
}
parameters {

//  parameters for first treatment

vector<lower=0>[6] r1;   // 6 growth terms Fixed Lower bound
vector<lower=',AlphaLowerBound,'>[36] A1;  //36 A terms. Fixed Lower bound
real<lower=0> phi1;  // divergence

//  parameters for second treatment

vector<lower=0>[6] r2;   // 6 growth terms Fixed Lower bound
vector<lower=',AlphaLowerBound,'>[36] A2;  //36 A terms. Fixed Lower bound
real<lower=0> phi2;  // divergence


}
model {

r1[1:6] ~ normal(10, 10);  // this "r" should be larger than the "r" using growthRate
A1[1:36] ~ normal(0, ',AlphaSD,');
r2[1:6] ~ normal(10, 10);
A2[1:36] ~ normal(0, ',AlphaSD,');
phi1 ~ cauchy(0,10);  // prior of phi (overdispersion controlling parameter)
phi2 ~ cauchy(0,10);
') 
  
  
  GQ<-'generated quantities {
  vector[N] log_lik; 
  vector[N] y_sim;\n'          # Start of generated quantitities section
  
  
  pars = 1:42 # For identifying
  
  #######
  rs <- pars[1:6]
  A <- array(pars[7:42], dim = c(6,6))
  A0 <- cbind(A, matrix(0,6,1))
  n=nrow(d)
  
  N0s<-  as.matrix(mutate(select(d, BIR_0:SUL_0), Void = 0))*2  ##### Using total counts, not pairs (BIR_0:SUL_0 are pairs, so need doubling)
  r_sel <- d$i  ## i is 1-6, denoting the index of the focal species
  
  ## Dataframe of indexes to use to select which parameters to apply for each data point
  X<-  data.frame(N0_i = N0s[as.matrix(data.frame(k=1:n,i=d$i))],     # Aii_selection1
                  N0_j = N0s[as.matrix(data.frame(k=1:n,j=d$j))],     # Aij_selection1
                  Aii_2  = A[as.matrix(data.frame(i1=d$i,i2=d$i))]-6 , # Aii_selection2
                  Aij_2 = A0[as.matrix(data.frame(i=d$i,j=d$j))]-6)    # Aij_selection2
  
  for(i in 1:n){
    if( d$Parasitism[i] == 'Safe'){
      GrowthCore = paste0( 'r1[',r_sel[i],']/(1+',X$N0_i[i],'*A1[',X$Aii_2[i],']')  # NB needs brackets of denominator closing!
      if(X$Aij_2[i] >0){
        GrowthCore <- paste0(GrowthCore, '+',X$N0_j[i],'*A1[',X$Aij_2[i],  ']')    # Add on other species  
      }
      GrowthCore <- paste0(X$N0_i[i],'*(',GrowthCore,')), phi1);\n') 
    }else{
      GrowthCore = paste0( 'r2[',r_sel[i],']/(1+',X$N0_i[i],'*A2[',X$Aii_2[i],']')  # NB needs brackets of denominator closing!
      if(X$Aij_2[i] >0){
        GrowthCore <- paste0(GrowthCore, '+',X$N0_j[i],'*A2[',X$Aij_2[i],  ']')    # Add on other species  
      }
      GrowthCore <- paste0(X$N0_i[i],'*(',GrowthCore,')), phi2);\n') 
    }
    Base <- paste0(Base, 'y[',i,'] ~  neg_binomial_2(',GrowthCore)
    GQ <- paste0(GQ, 'y_sim[',i,'] = neg_binomial_2_rng(',GrowthCore,
                 'log_lik[',i,'] = neg_binomial_2_lpmf(y[',i,']|', GrowthCore)
  }
  
  Both <- paste0(Base,'\n} \n', ifelse(CalcGQ, paste0(GQ,'\n}'), ''))
  writeLines(Both, paste0('StanModels/BuiltModel_NegB_',Name,'.stan'))
  return(paste0('Function Finished. Output: <BuiltModel_NegB_',Name,'.stan>'))
}


Write_Stan_NegB_2R_1A<-function(d,
                                Name,
                                AlphaLowerBound = 0,
                                AlphaSD= 1, 
                                CalcGQ = TRUE){
  
  # Function to generate a Stan model
  # uses raw count data and a Negative-Binomial error distrubtion
  # Fits one set of A's to the two treatments 
  
  #############
  ### Core Start
  #########
  
  Base<- paste0('data {
  int<lower=0> N;
  int y[N];// Observed Count of emergences
}
parameters {


vector<lower=',AlphaLowerBound,'>[36] A;  //36 A terms. Fixed Lower bound

//  parameters for first treatment

vector<lower=0>[6] r1;   // 6 growth terms Fixed Lower bound
real<lower=0> phi1;  // divergence

//  parameters for second treatment

vector<lower=0>[6] r2;   // 6 growth terms Fixed Lower bound
real<lower=0> phi2;  // divergence


}
model {

A[1:36] ~ normal(0, ',AlphaSD,');
r1[1:6] ~ normal(10, 10);  // this "r" should be larger than the "r" using growthRate
r2[1:6] ~ normal(10, 10);
phi1 ~ cauchy(0,10);  // prior of phi (overdispersion controlling parameter)
phi2 ~ cauchy(0,10);
') 
  
  
  GQ<-'generated quantities {
  vector[N] log_lik; 
  vector[N] y_sim;\n'          # Start of generated quantitities section
  
  
  pars = 1:42 # For identifying
  
  #######
  rs <- pars[1:6]
  A <- array(pars[7:42], dim = c(6,6))
  A0 <- cbind(A, matrix(0,6,1))
  n=nrow(d)
  
  N0s<-  as.matrix(mutate(select(d, BIR_0:SUL_0), Void = 0))*2  ##### Using total counts, not pairs (BIR_0:SUL_0 are pairs, so need doubling)
  r_sel <- d$i  ## i is 1-6, denoting the index of the focal species
  
  ## Dataframe of indexes to use to select which parameters to apply for each data point
  X<-  data.frame(N0_i = N0s[as.matrix(data.frame(k=1:n,i=d$i))],     # Aii_selection1
                  N0_j = N0s[as.matrix(data.frame(k=1:n,j=d$j))],     # Aij_selection1
                  Aii_2  = A[as.matrix(data.frame(i1=d$i,i2=d$i))]-6 , # Aii_selection2
                  Aij_2 = A0[as.matrix(data.frame(i=d$i,j=d$j))]-6)    # Aij_selection2
  
  for(i in 1:n){
    if( d$Parasitism[i] == 'Safe'){
      GrowthCore = paste0( 'r1[',r_sel[i],']/(1+',X$N0_i[i],'*A[',X$Aii_2[i],']')  # NB needs brackets of denominator closing!
      if(X$Aij_2[i] >0){
        GrowthCore <- paste0(GrowthCore, '+',X$N0_j[i],'*A[',X$Aij_2[i],  ']')    # Add on other species  
      }
      GrowthCore <- paste0(X$N0_i[i],'*(',GrowthCore,')), phi1);\n') 
    }else{
      GrowthCore = paste0( 'r2[',r_sel[i],']/(1+',X$N0_i[i],'*A[',X$Aii_2[i],']')  # NB needs brackets of denominator closing!
      if(X$Aij_2[i] >0){
        GrowthCore <- paste0(GrowthCore, '+',X$N0_j[i],'*A[',X$Aij_2[i],  ']')    # Add on other species  
      }
      GrowthCore <- paste0(X$N0_i[i],'*(',GrowthCore,')), phi2);\n') 
    }
    Base <- paste0(Base, 'y[',i,'] ~  neg_binomial_2(',GrowthCore)
    GQ <- paste0(GQ, 'y_sim[',i,'] = neg_binomial_2_rng(',GrowthCore,
                 'log_lik[',i,'] = neg_binomial_2_lpmf(y[',i,']|', GrowthCore)
  }
  
  Both <- paste0(Base,'\n} \n', ifelse(CalcGQ, paste0(GQ,'\n}'), ''))
  writeLines(Both, paste0('StanModels/BuiltModel_NegB_',Name,'.stan'))
  return(paste0('Function Finished. Output: <BuiltModel_NegB_',Name,'.stan>'))
}



Write_Stan_NegB_1R_1A<-function(d,
                                Name,
                                AlphaLowerBound = 0,
                                AlphaSD= 1, 
                                CalcGQ = TRUE){
  
  # Function to generate a Stan model
  # uses raw count data and a Negative-Binomial error distribution
  # Fits one set of A's and R's to the two treatments 
  
  #############
  ### Core Start
  #########
  
  Base<- paste0('data {
  int<lower=0> N;
  int y[N];// Observed Count of emergences
}
parameters {


vector<lower=',AlphaLowerBound,'>[36] A;  //36 A terms. Fixed Lower bound
vector<lower=0>[6] r;   // 6 growth terms Fixed Lower bound

//  parameters for first treatment
real<lower=0> phi1;  // divergence

//  parameters for second treatment
real<lower=0> phi2;  // divergence


}
model {

A[1:36] ~ normal(0, ',AlphaSD,');
r[1:6] ~ normal(10, 10);  // this "r" should be larger than the "r" using growthRate
phi1 ~ cauchy(0,10);  // prior of phi (overdispersion controlling parameter)
phi2 ~ cauchy(0,10);
') 
  
  GQ<-'generated quantities {
  vector[N] log_lik; 
  vector[N] y_sim;\n'          # Start of generated quantitities section
  
  
  pars = 1:42 # For identifying
  
  #######
  rs <- pars[1:6]
  A <- array(pars[7:42], dim = c(6,6))
  A0 <- cbind(A, matrix(0,6,1))
  n=nrow(d)
  
  N0s<-  as.matrix(mutate(select(d, BIR_0:SUL_0), Void = 0))*2  ##### Using total counts, not pairs (BIR_0:SUL_0 are pairs, so need doubling)
  r_sel <- d$i  ## i is 1-6, denoting the index of the focal species
  
  ## Dataframe of indexes to use to select which parameters to apply for each data point
  X<-  data.frame(N0_i = N0s[as.matrix(data.frame(k=1:n,i=d$i))],     # Aii_selection1
                  N0_j = N0s[as.matrix(data.frame(k=1:n,j=d$j))],     # Aij_selection1
                  Aii_2  = A[as.matrix(data.frame(i1=d$i,i2=d$i))]-6 , # Aii_selection2
                  Aij_2 = A0[as.matrix(data.frame(i=d$i,j=d$j))]-6)    # Aij_selection2
  
  for(i in 1:n){
    if( d$Parasitism[i] == 'Safe'){
      GrowthCore = paste0( 'r[',r_sel[i],']/(1+',X$N0_i[i],'*A[',X$Aii_2[i],']')  # NB needs brackets of denominator closing!
      if(X$Aij_2[i] >0){
        GrowthCore <- paste0(GrowthCore, '+',X$N0_j[i],'*A[',X$Aij_2[i],  ']')    # Add on other species  
      }
      GrowthCore <- paste0(X$N0_i[i],'*(',GrowthCore,')), phi1);\n') 
    }else{
      GrowthCore = paste0( 'r[',r_sel[i],']/(1+',X$N0_i[i],'*A[',X$Aii_2[i],']')  # NB needs brackets of denominator closing!
      if(X$Aij_2[i] >0){
        GrowthCore <- paste0(GrowthCore, '+',X$N0_j[i],'*A[',X$Aij_2[i],  ']')    # Add on other species  
      }
      GrowthCore <- paste0(X$N0_i[i],'*(',GrowthCore,')), phi2);\n') 
    }
    Base <- paste0(Base, 'y[',i,'] ~  neg_binomial_2(',GrowthCore)
    GQ <- paste0(GQ, 'y_sim[',i,'] = neg_binomial_2_rng(',GrowthCore,
                 'log_lik[',i,'] = neg_binomial_2_lpmf(y[',i,']|', GrowthCore)
  }
  
  Both <- paste0(Base,'\n} \n', ifelse(CalcGQ, paste0(GQ,'\n}'), ''))
  writeLines(Both, paste0('StanModels/BuiltModel_NegB_',Name,'.stan'))
  return(paste0('Function Finished. Output: <BuiltModel_NegB_',Name,'.stan>'))
}

