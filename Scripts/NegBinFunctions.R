
Write_Stan_NegB_DoubleA<-function(PreparedData1,PreparedData2, Name,
                                  AlphaLowerBound = TRUE,
                                  AlphaSD= 1){
  
  # Function to Fit a Stan, using raw count data and a Negative-Binomial error distrubtion
  ### This one fits A's to the two treatments totally seperately
  
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
  vector',ifelse(AlphaLowerBound,'<lower=0>',''),'[36] A1;  //36 A terms. Fixed Lower bound
  real<lower=0> phi1;  // divergence
  
//  parameters for second treatment

  vector<lower=0>[6] r2;   // 6 growth terms Fixed Lower bound
  vector',ifelse(AlphaLowerBound,'<lower=0>',''),'[36] A2;  //36 A terms. Fixed Lower bound
  real<lower=0> phi2;  // divergence

}
model {
r1[1:6] ~ normal(10, 10);  // this "r" should be larger than the "r" using growthRate
A1[1:36] ~ normal(0, ',AlphaSD,');
r2[1:6] ~ normal(10, 10);
A2[1:36] ~ normal(0, ',AlphaSD,');
phi1 ~ cauchy(0,10);  // prior of phi
phi2 ~ cauchy(0,10);
') 
  
  
  GQ<-'generated quantities {
  vector[N] log_lik;
  vector[N] y_sim;\n'          # Start off generated quantitities section
  
  
  pars = 1:42 # For identifying
  
  #######
  rs <- pars[1:6]
  A <- array(pars[7:42], dim = c(6,6))
  A0 <- cbind(A, matrix(0,6,1))
  
  
  ############# Do first treatment ######
  #######################################
  
  d=PreparedData1
  n=nrow(d)
  
  n_firstset <- n
  
  N0s<-  as.matrix(mutate(select(d, BIR_0:SUL_0), Void = 0))*2  ##### Using total counts, not pairs (BIR_0:SUL_0 are pairs)
  
  Aii_selection1<- as.matrix(data.frame(k=1:n,i=d$i))
  Aii_selection2<- as.matrix(data.frame(i1=d$i,i2=d$i))
  
  Aij_selection1<- as.matrix(data.frame(k=1:n,j=d$j))
  Aij_selection2<- as.matrix(data.frame(i=d$i,j=d$j))
  
  r_sel <- d$i
  
  R<-  data.frame(N0_i = N0s[Aii_selection1],
                  N0_j = N0s[Aij_selection1],
                  Aii_2 = A[Aii_selection2]-6 ,
                  Aij_2 = A0[Aij_selection2]-6)
  
  
  
  for(i in 1:n){
    if(R$Aij_2[i] >0){
      Base <- paste0(Base,   'target += neg_binomial_2_lpmf(y[',i,']|    (', R$N0_i[i],'*r1[',r_sel[i],'])/(1+',  R$N0_i[i],'*A1[',R$Aii_2[i],']+', R$N0_j[i],'*A1[',R$Aij_2[i],']), phi1);\n'  )
      GQ <- paste0(GQ, 'log_lik[',i,']= neg_binomial_2_lpmf(y[',i,']|    (', R$N0_i[i],'*r1[',r_sel[i],'])/(1+',  R$N0_i[i],'*A1[',R$Aii_2[i],']+', R$N0_j[i],'*A1[',R$Aij_2[i],']), phi1);\n',
                   'y_sim[',i,'] = neg_binomial_2_rng((', R$N0_i[i],'*r1[',r_sel[i],'])/(1+',  R$N0_i[i],'*A1[',R$Aii_2[i],']+', R$N0_j[i],'*A1[',R$Aij_2[i],']), phi1);\n')
    }else{
      Base <-paste0(Base,     'target += neg_binomial_2_lpmf(y[',i,']|(', R$N0_i[i],'*r1[',r_sel[i],'])/(1+',  R$N0_i[i],'*A1[',R$Aii_2[i],']), phi1);\n'  )
      GQ <- paste0(GQ, 'log_lik[',i,'] = neg_binomial_2_lpmf(y[',i,']|(', R$N0_i[i],'*r1[',r_sel[i],'])/(1+',  R$N0_i[i],'*A1[',R$Aii_2[i],']), phi1);\n',
                   'y_sim[',i,'] = neg_binomial_2_rng((', R$N0_i[i],'*r1[',r_sel[i],'])/(1+',  R$N0_i[i],'*A1[',R$Aii_2[i],']), phi1);\n')
    }
  }
  
  Base <-paste0(Base,'\n // Start of second treatment  \n ')
  
  ############# Do second treatment ######
  #######################################
  
  d=PreparedData2
  
  n=nrow(d)
  
  N0s<-  as.matrix(mutate(select(d, BIR_0:SUL_0), Void = 0))*2  ##### Using total counts, not pairs (BIR_0:SUL_0 are pairs)
  
  Aii_selection1<- as.matrix(data.frame(k=1:n,i=d$i))
  Aii_selection2<- as.matrix(data.frame(i1=d$i,i2=d$i))
  
  Aij_selection1<- as.matrix(data.frame(k=1:n,j=d$j))
  Aij_selection2<- as.matrix(data.frame(i=d$i,j=d$j))
  
  r_sel <- d$i
  
  R<-  data.frame(N0_i = N0s[Aii_selection1],
                  N0_j = N0s[Aij_selection1],
                  Aii_2 = A[Aii_selection2]-6 ,
                  Aij_2 = A0[Aij_selection2]-6)
  
  for(i in 1:n){
    if(R$Aij_2[i] >0){
      Base <- paste0(Base, 'target += neg_binomial_2_lpmf(y[',i+n_firstset,']| (', R$N0_i[i],'*r2[',r_sel[i],'])/(1+',  R$N0_i[i],'*A2[',R$Aii_2[i],']+', R$N0_j[i],'*A2[',R$Aij_2[i],']), phi2);\n'  )
      GQ <- paste0(GQ,'log_lik[',i+n_firstset,'] = neg_binomial_2_lpmf(y[',i+n_firstset,']| (', R$N0_i[i],'*r2[',r_sel[i],'])/(1+',  R$N0_i[i],'*A2[',R$Aii_2[i],']+', R$N0_j[i],'*A2[',R$Aij_2[i],']), phi2);\n',
                   'y_sim[',i+n_firstset,'] = neg_binomial_2_rng((',R$N0_i[i],'*r2[',r_sel[i],'])/(1+',  R$N0_i[i],'*A2[',R$Aii_2[i],']+', R$N0_j[i],'*A2[',R$Aij_2[i],']), phi2);\n')
    }else{
      Base <-paste0(Base,                   'target += neg_binomial_2_lpmf(y[',i+n_firstset,']|(', R$N0_i[i],'*r2[',r_sel[i],'])/(1+',  R$N0_i[i],'*A2[',R$Aii_2[i],']), phi2);\n'  )
      GQ <- paste0(GQ,   'log_lik[',i+n_firstset,'] =  neg_binomial_2_lpmf(y[',i+n_firstset,']|(', R$N0_i[i],'*r2[',r_sel[i],'])/(1+',  R$N0_i[i],'*A2[',R$Aii_2[i],']), phi2);\n',
                   'y_sim[',i+n_firstset,'] = neg_binomial_2_rng((',R$N0_i[i],'*r2[',r_sel[i],'])/(1+',  R$N0_i[i],'*A2[',R$Aii_2[i],']), phi2);\n')
    }
  }
  
  Both <- paste0(Base,'\n} \n', GQ, '\n}')
  writeLines(Both, paste0('StanModels/BuiltModel_',Name,'.stan'))
  return(paste0('Function Finished. Output: <BuiltModel_',Name,'.stan>'))
}


Write_Stan_NegB_SingleA<-function(PreparedData1,PreparedData2, Name,
                                  AlphaLowerBound = TRUE,
                                  AlphaSD= 1){
  
  # Function to Fit a Stan 
  ### This one fits one set of A's to the two treatments 
  
  #############
  ### Core Start
  #########
  
  Base<- paste0('data {
  int<lower=0> N;
  int y[N]; // Observed Count of emergences
}
parameters {

//  parameters for first treatment

  vector<lower=0>[6] r1;   // 6 growth terms Fixed Lower bound
  real<lower=0> phi1;
  
//  parameters for second treatment

  vector<lower=0>[6] r2;   // 6 growth terms Fixed Lower bound
  real<lower=0> phi2;

// For Both  
  
  vector',ifelse(AlphaLowerBound,'<lower=0>',''),'[36] A;  //36 A terms. Fixed Lower bound

}
model {
r1[1:6] ~ normal(10, 10);
r2[1:6] ~ normal(10, 10);
phi1 ~ cauchy(0,10);
phi2 ~ cauchy(0,10);
A[1:36] ~ normal(0, ',AlphaSD,');
' ) 
  
  
  GQ<-'generated quantities {
  vector[N] log_lik;
  vector[N] y_sim;\n'          # Start off generated quantitities section
  
  pars = 1:42 # For identifying
  
  #######
  rs <- pars[1:6]
  A <- array(pars[7:42], dim = c(6,6))
  A0 <- cbind(A, matrix(0,6,1))
  
  
  ############# Do first treatment ######
  #######################################
  
  d=PreparedData1
  n=nrow(d)
  
  n_firstset <- n
  
  N0s<-  as.matrix(mutate(select(d, BIR_0:SUL_0), Void = 0))*2  ##### Using total counts, not pairs (BIR_0:SUL_0 are pairs)
  
  Aii_selection1<- as.matrix(data.frame(k=1:n,i=d$i))
  Aii_selection2<- as.matrix(data.frame(i1=d$i,i2=d$i))
  
  Aij_selection1<- as.matrix(data.frame(k=1:n,j=d$j))
  Aij_selection2<- as.matrix(data.frame(i=d$i,j=d$j))
  
  r_sel <- d$i
  
  R<-  data.frame(N0_i = N0s[Aii_selection1],
                  N0_j = N0s[Aij_selection1],
                  Aii_2 = A[Aii_selection2]-6 ,
                  Aij_2 = A0[Aij_selection2]-6)
  
  for(i in 1:n){
    if(R$Aij_2[i] >0){
      Base <- paste0(Base,    'target += neg_binomial_2_lpmf(y[',i,']|    (', R$N0_i[i],'*r1[',r_sel[i],'])/(1+',  R$N0_i[i],'*A[',R$Aii_2[i],']+', R$N0_j[i],'*A[',R$Aij_2[i],']), phi1);\n'  )
      GQ <- paste0(GQ, 'log_lik[',i,'] = neg_binomial_2_lpmf(y[',i,']|    (', R$N0_i[i],'*r1[',r_sel[i],'])/(1+',  R$N0_i[i],'*A[',R$Aii_2[i],']+', R$N0_j[i],'*A[',R$Aij_2[i],']), phi1);\n',
                   'y_sim[',i,'] = neg_binomial_2_rng((', R$N0_i[i],'*r1[',r_sel[i],'])/(1+',  R$N0_i[i],'*A[',R$Aii_2[i],']+', R$N0_j[i],'*A[',R$Aij_2[i],']), phi1);\n')
      
    }else{
      Base <-paste0(Base,     'target += neg_binomial_2_lpmf(y[',i,']| (', R$N0_i[i],'*r1[',r_sel[i],'])/(1+',  R$N0_i[i],'*A[',R$Aii_2[i],']), phi1);\n'  )
      GQ <- paste0(GQ, 'log_lik[',i,'] = neg_binomial_2_lpmf(y[',i,']| (', R$N0_i[i],'*r1[',r_sel[i],'])/(1+',  R$N0_i[i],'*A[',R$Aii_2[i],']), phi1);\n',
                   'y_sim[',i,'] = neg_binomial_2_rng((',R$N0_i[i],'*r1[',r_sel[i],'])/(1+',  R$N0_i[i],'*A[',R$Aii_2[i],']), phi1);\n')
    }
  }
  
  Base <-paste0(Base,'\n // Start of second treatment \n ')
  
  ############# Do second treatment ######
  #######################################
  
  d=PreparedData2
  
  n=nrow(d)
  
  N0s<-  as.matrix(mutate(select(d, BIR_0:SUL_0), Void = 0))*2  ##### Using total counts, not pairs (BIR_0:SUL_0 are pairs)
  
  Aii_selection1<- as.matrix(data.frame(k=1:n,i=d$i))
  Aii_selection2<- as.matrix(data.frame(i1=d$i,i2=d$i))
  
  Aij_selection1<- as.matrix(data.frame(k=1:n,j=d$j))
  Aij_selection2<- as.matrix(data.frame(i=d$i,j=d$j))
  
  r_sel <- d$i
  
  R<-  data.frame(N0_i = N0s[Aii_selection1],
                  N0_j = N0s[Aij_selection1],
                  Aii_2 = A[Aii_selection2]-6 ,
                  Aij_2 = A0[Aij_selection2]-6)
  
  for(i in 1:n){
    if(R$Aij_2[i] >0){
      Base <- paste0(Base,                 'target += neg_binomial_2_lpmf(y[',i+n_firstset,']|   (', R$N0_i[i],'*r2[',r_sel[i],'])/(1+',  R$N0_i[i],'*A[',R$Aii_2[i],']+', R$N0_j[i],'*A[',R$Aij_2[i],']), phi2);\n'  )
      GQ <- paste0(GQ,   'log_lik[',i+n_firstset,'] = neg_binomial_2_lpmf(y[',i+n_firstset,']|   (', R$N0_i[i],'*r2[',r_sel[i],'])/(1+',  R$N0_i[i],'*A[',R$Aii_2[i],']+', R$N0_j[i],'*A[',R$Aij_2[i],']), phi2);\n',
                   'y_sim[',i+n_firstset,'] = neg_binomial_2_rng((',R$N0_i[i],'*r2[',r_sel[i],'])/(1+',  R$N0_i[i],'*A[',R$Aii_2[i],']+', R$N0_j[i],'*A[',R$Aij_2[i],']), phi2);\n')
    }else{
      Base <-paste0(Base,                  'target += neg_binomial_2_lpmf(y[',i+n_firstset,']| (', R$N0_i[i],'*r2[',r_sel[i],'])/(1+',  R$N0_i[i],'*A[',R$Aii_2[i],']), phi2);\n'  )
      GQ <- paste0(GQ,   'log_lik[',i+n_firstset,'] = neg_binomial_2_lpmf(y[',i+n_firstset,']| (', R$N0_i[i],'*r2[',r_sel[i],'])/(1+',  R$N0_i[i],'*A[',R$Aii_2[i],']), phi2);\n',
                   'y_sim[',i+n_firstset,'] = neg_binomial_2_rng((',R$N0_i[i],'*r2[',r_sel[i],'])/(1+',  R$N0_i[i],'*A[',R$Aii_2[i],']), phi2);\n')
    }
  }
  
  #######
  ### LOO needs log lik in generated quantitites. 
  ###### 
  
  
  Both <- paste0(Base,'\n} \n', GQ, '\n}')
  
  
  writeLines(Both, paste0('StanModels/BuiltModel_',Name,'.stan'))
  return(paste0('Function Finished. Output: <BuiltModel_',Name,'.stan>'))
}