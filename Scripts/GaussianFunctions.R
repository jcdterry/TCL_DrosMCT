## Functions for Fitting MCT fly data

WriteStanFunction<-function(PreparedData, Name,
                            AlphaLowerBound = TRUE,
                            AlphaSD= 1, 
                            log_alpha= FALSE){
  
  
  d=PreparedData
  
  n=nrow(d)
  
  print(paste0(Name,' seems to have ',n,' entries'))
  
  N0s<-  as.matrix(mutate(select(d, BIR_0:SUL_0), Void = 0))
  
  Aii_selection1<- as.matrix(data.frame(k=1:n,i=d$i))
  Aii_selection2<- as.matrix(data.frame(i1=d$i,i2=d$i))
  
  Aij_selection1<- as.matrix(data.frame(k=1:n,j=d$j))
  Aij_selection2<- as.matrix(data.frame(i=d$i,j=d$j))
  
  r_sel <- d$i
  
  ### ^^^ These values are all fixed, and could be hardcoded into the function
  
  pars = 1:42 # For identifying
  
  #######
  rs <- pars[1:6]
  A <- array(pars[7:42], dim = c(6,6))
  A0 <- cbind(A, matrix(0,6,1))
  
  R<-  data.frame(N0_i = N0s[Aii_selection1],
                  N0_j = N0s[Aij_selection1],
                  Aii_2 = A[Aii_selection2]-6 ,
                  Aij_2 = A0[Aij_selection2]-6)
  
  Base<- paste0('data {
  int<lower=0> N;
  vector[N] y; // Observed Growth Rate
}
parameters {
  vector<lower=0>[6] r;   // 6 growth terms Fixed Lower bound
  vector',ifelse(AlphaLowerBound,'<lower=0>',''),'[36] A;  //36 A terms. Fixed Lower bound
  real<lower=0> sigma;   // Overall error term
}
model {
r[1:6] ~ normal(1, 1);
A[1:36] ~ normal(0, ',AlphaSD,');
' ) 

  if(log_alpha){
    
    for(i in 1:n){
      if(R$Aij_2[i] >0){
        Base <- paste0(Base, 'target += normal_lpdf(y[',i,']|    r[',r_sel[i],']/(1+',  R$N0_i[i],'*exp(A[',R$Aii_2[i],'])+', R$N0_j[i],'*exp(A[',R$Aij_2[i],'])),  sigma );\n'  )
      }else{
        Base <-paste0(Base,'target += normal_lpdf(y[',i,']| r[',r_sel[i],']/(1+',  R$N0_i[i],'*exp(A[',R$Aii_2[i],'])),  sigma );\n'  )
      }
    } 
  } else{
    
    for(i in 1:n){
      if(R$Aij_2[i] >0){
        Base <- paste0(Base, 'target += normal_lpdf(y[',i,']|    r[',r_sel[i],']/(1+',  R$N0_i[i],'*A[',R$Aii_2[i],']+', R$N0_j[i],'*A[',R$Aij_2[i],']),  sigma );\n'  )
      }else{
        Base <-paste0(Base,'target += normal_lpdf(y[',i,']| r[',r_sel[i],']/(1+',  R$N0_i[i],'*A[',R$Aii_2[i],']),  sigma );\n'  )
      }
    }
  }
Base <- paste0(Base, '}\n')

Base <- paste0(Base, "generated quantities{
real y_sim[N];\n")
for(i in 1:n){
  if(R$Aij_2[i] >0){
    Base <- paste0(Base, 'y_sim[', i, '] = normal_rng(r[',r_sel[i],']/(1+',  R$N0_i[i],'*(A[',R$Aii_2[i],'])+', R$N0_j[i],'*(A[',R$Aij_2[i],'])),  sigma );\n')
  }else{
    Base <- paste0(Base, 'y_sim[', i, '] = normal_rng(r[',r_sel[i],']/(1+',  R$N0_i[i],'*(A[',R$Aii_2[i],'])),  sigma );\n')
  }
}
Base <- paste0(Base, '}\n')
writeLines(Base, paste0('StanModels/BuiltModel_',Name,'.stan'))


return(paste0('Function Finished. Output: <BuiltModel_',Name,'.stan>'))
}






Write_BothTreatmentsStanFunction_SingleA<-function(PreparedData1,PreparedData2, Name,
                                                   AlphaLowerBound = TRUE,
                                                   AlphaSD= 1){
  
  # Function to Fit a Stan 
  ### This one fits one set of A's to the two treatments 
  ### Seperate error term  for the sets is maintained (need to think about how to handle this)
  
  
  #############
  ### Core Start
  #########
  
  Base<- paste0('data {
  int<lower=0> N;
  vector[N] y; // Observed Growth Rate
}
parameters {

//  parameters for first treatment

  vector<lower=0>[6] r1;   // 6 growth terms Fixed Lower bound
  real<lower=0> sigma1;   // Overall error term
  
  
//  parameters for second treatment

  vector<lower=0>[6] r2;   // 6 growth terms Fixed Lower bound
  real<lower=0> sigma2;   // Overall error term
  
// For Both  
  
  vector',ifelse(AlphaLowerBound,'<lower=0>',''),'[36] A;  //36 A terms. Fixed Lower bound

}
model {
r1[1:6] ~ normal(1, 1);
r2[1:6] ~ normal(1, 1);

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

N0s<-  as.matrix(mutate(select(d, BIR_0:SUL_0), Void = 0))

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
    Base <- paste0(Base,    'target += normal_lpdf(y[',i,']|r1[',r_sel[i],']/(1+',  R$N0_i[i],'*A[',R$Aii_2[i],']+', R$N0_j[i],'*A[',R$Aij_2[i],']),  sigma1 );\n'  )
    GQ <- paste0(GQ, 'log_lik[',i,'] = normal_lpdf(y[',i,']|r1[',r_sel[i],']/(1+',  R$N0_i[i],'*A[',R$Aii_2[i],']+', R$N0_j[i],'*A[',R$Aij_2[i],']),  sigma1 );\n',
                     'y_sim[', i, '] = normal_rng(          r1[',r_sel[i], ']/(1+', R$N0_i[i],'*A[',R$Aii_2[i],']+', R$N0_j[i],'*A[',R$Aij_2[i],']),  sigma1 );\n')
    
    }else{
    Base <-paste0(Base,     'target += normal_lpdf(y[',i,']| r1[',r_sel[i],']/(1+',  R$N0_i[i],'*A[',R$Aii_2[i],']),  sigma1 );\n'  )
    GQ <- paste0(GQ, 'log_lik[',i,'] = normal_lpdf(y[',i,']| r1[',r_sel[i],']/(1+',  R$N0_i[i],'*A[',R$Aii_2[i],']),  sigma1 );\n',
                    'y_sim[', i, '] = normal_rng(           r1[', r_sel[i],']/(1+',  R$N0_i[i],'*A[',R$Aii_2[i],']),  sigma1 );\n')
  }
}

Base <-paste0(Base,'\n // Start of second treatment \n ')

############# Do second treatment ######
#######################################

d=PreparedData2

n=nrow(d)

N0s<-  as.matrix(mutate(select(d, BIR_0:SUL_0), Void = 0))

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
    Base <- paste0(Base, 'target += normal_lpdf(y[',i+n_firstset,']|                r2[',r_sel[i],']/(1+',  R$N0_i[i],'*A[',R$Aii_2[i],']+', R$N0_j[i],'*A[',R$Aij_2[i],']),  sigma2 );\n'  )
    GQ <- paste0(GQ,   'log_lik[',i+n_firstset,'] = normal_lpdf(y[',i+n_firstset,']|r2[',r_sel[i],']/(1+',  R$N0_i[i],'*A[',R$Aii_2[i],']+', R$N0_j[i],'*A[',R$Aij_2[i],']),  sigma2 );\n',
                       'y_sim[', i+n_firstset, '] = normal_rng(                     r2[', r_sel[i],']/(1+', R$N0_i[i],'*A[',R$Aii_2[i],']+', R$N0_j[i],'*A[',R$Aij_2[i],']),  sigma2 );\n')
  }else{
    Base <-paste0(Base,'target += normal_lpdf(y[',i+n_firstset,']|                  r2[',r_sel[i],']/(1+',  R$N0_i[i],'*A[',R$Aii_2[i],']),  sigma2 );\n'  )
    GQ <- paste0(GQ,   'log_lik[',i+n_firstset,'] = normal_lpdf(y[',i+n_firstset,']|r2[',r_sel[i],']/(1+',  R$N0_i[i],'*A[',R$Aii_2[i],']),  sigma2 );\n',
                       'y_sim[', i+n_firstset, '] = normal_rng(                     r2[',r_sel[i],']/(1+',  R$N0_i[i],'*A[',R$Aii_2[i],']),  sigma2 );\n')
  }
}

#######
### LOO needs log lik in generated quantitites. 
###### 


Both <- paste0(Base,'\n} \n', GQ, '\n}')


writeLines(Both, paste0('StanModels/BuiltModel_',Name,'.stan'))
return(paste0('Function Finished. Output: <BuiltModel_',Name,'.stan>'))
}

Write_BothTreatmentsStanFunction_DoubleA<-function(PreparedData1,PreparedData2, Name,
                                                   AlphaLowerBound = TRUE,
                                                   AlphaSD= 1){
  
  # Function to Fit a Stan 
  ### This one fits A's to the two treatments totally seperately
  ### Seperate error term too (need to think about how to handle this)
  
  
  #############
  ### Core Start
  #########
  
  Base<- paste0('data {
  int<lower=0> N;
  vector[N] y; // Observed Growth Rate
}
parameters {

//  parameters for first treatment

  vector<lower=0>[6] r1;   // 6 growth terms Fixed Lower bound
  vector',ifelse(AlphaLowerBound,'<lower=0>',''),'[36] A1;  //36 A terms. Fixed Lower bound
  real<lower=0> sigma1;   // Overall error term
  
  
//  parameters for second treatment

  vector<lower=0>[6] r2;   // 6 growth terms Fixed Lower bound
  vector',ifelse(AlphaLowerBound,'<lower=0>',''),'[36] A2;  //36 A terms. Fixed Lower bound
  real<lower=0> sigma2;   // Overall error term
  

}
model {
r1[1:6] ~ normal(1, 1);
A1[1:36] ~ normal(0, ',AlphaSD,');
r2[1:6] ~ normal(1, 1);
A2[1:36] ~ normal(0, ',AlphaSD,');
// why does sigma have no prior?
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

N0s<-  as.matrix(mutate(select(d, BIR_0:SUL_0), Void = 0))

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
    Base <- paste0(Base,   'target += normal_lpdf(y[',i,']|    r1[',r_sel[i],']/(1+',  R$N0_i[i],'*A1[',R$Aii_2[i],']+', R$N0_j[i],'*A1[',R$Aij_2[i],']),  sigma1 );\n'  )
    GQ <- paste0(GQ, 'log_lik[',i,']= normal_lpdf(y[',i,']|    r1[',r_sel[i],']/(1+',  R$N0_i[i],'*A1[',R$Aii_2[i],']+', R$N0_j[i],'*A1[',R$Aij_2[i],']),  sigma1 );\n',
                     'y_sim[', i, "] = normal_rng(             r1[",r_sel[i],']/(1+',  R$N0_i[i],'*A1[',R$Aii_2[i],']+', R$N0_j[i],'*A1[',R$Aij_2[i],']),  sigma1 );\n')
  }else{
    Base <-paste0(Base,     'target += normal_lpdf(y[',i,']| r1[',r_sel[i],']/(1+',  R$N0_i[i],'*A1[',R$Aii_2[i],']),  sigma1 );\n'  )
    GQ <- paste0(GQ, 'log_lik[',i,'] = normal_lpdf(y[',i,']| r1[',r_sel[i],']/(1+',  R$N0_i[i],'*A1[',R$Aii_2[i],']),  sigma1 );\n',
                     "y_sim[", i, "] = normal_rng(           r1[",r_sel[i],']/(1+',  R$N0_i[i],'*A1[',R$Aii_2[i],']),  sigma1 );\n')
  }
}

Base <-paste0(Base,'\n // Start of second treatment  \n ')

############# Do second treatment ######
#######################################

d=PreparedData2

n=nrow(d)

N0s<-  as.matrix(mutate(select(d, BIR_0:SUL_0), Void = 0))

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
    Base <- paste0(Base,'target += normal_lpdf(y[',i+n_firstset,']|                     r2[',r_sel[i],']/(1+',  R$N0_i[i],'*A2[',R$Aii_2[i],']+', R$N0_j[i],'*A2[',R$Aij_2[i],']),  sigma2 );\n'  )
    GQ <- paste0(GQ,   'log_lik[',i+n_firstset,'] = normal_lpdf(y[',i+n_firstset,']|    r2[',r_sel[i],']/(1+',  R$N0_i[i],'*A2[',R$Aii_2[i],']+', R$N0_j[i],'*A2[',R$Aij_2[i],']),  sigma2 );\n',
                         "y_sim[", i+n_firstset, "] = normal_rng(                       r2[",r_sel[i],']/(1+',  R$N0_i[i],'*A2[',R$Aii_2[i],']+', R$N0_j[i],'*A2[',R$Aij_2[i],']),  sigma2 );\n')
  }else{
    Base <-paste0(Base,'target += normal_lpdf(y[',i+n_firstset,']|                    r2[',r_sel[i],']/(1+',  R$N0_i[i],'*A2[',R$Aii_2[i],']),  sigma2 );\n'  )
    GQ <- paste0(GQ,   'log_lik[',i+n_firstset,'] =  normal_lpdf(y[',i+n_firstset,']| r2[',r_sel[i],']/(1+',  R$N0_i[i],'*A2[',R$Aii_2[i],']),  sigma2 );\n',
                        "y_sim[", i+n_firstset,"] = normal_rng(                       r2[",r_sel[i],']/(1+',  R$N0_i[i],'*A2[',R$Aii_2[i],']),  sigma2 );\n')
  }
}

Both <- paste0(Base,'\n} \n', GQ, '\n}')
writeLines(Both, paste0('StanModels/BuiltModel_',Name,'.stan'))
return(paste0('Function Finished. Output: <BuiltModel_',Name,'.stan>'))
}




Write_BothTreatmentsStanFunction_SingleSetRandA<-function(PreparedData1,PreparedData2, Name,
                                                          AlphaLowerBound = TRUE,
                                                          AlphaSD= 1){
  
  # Function to Fit a Stan 
  ### This one fits one set of A's to the two treatments  and also one set of R's
  ### Seperate error term  for the sets is maintained (need to think about how to handle this)
  
  
  #############
  ### Core Start
  #########
  
  Base<- paste0('data {
  int<lower=0> N;
  vector[N] y; // Observed Growth Rate
}
parameters {

//  parameters for first treatment

  real<lower=0> sigma1;   // Overall error term
  
  
//  parameters for second treatment

  real<lower=0> sigma2;   // Overall error term
  
// For Both  
  
  vector',ifelse(AlphaLowerBound,'<lower=0>',''),'[36] A;  //36 A terms. Fixed Lower bound
  vector<lower=0>[6] r;   // 6 growth terms Fixed Lower bound

}
model {
r[1:6] ~ normal(1, 1);
A[1:36] ~ normal(0, ',AlphaSD,');
')  

GQ<-'generated quantities {
  vector[N] log_lik; 
  vector[N] y_sim;  \n'          # Start off generated quantitities section

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

N0s<-  as.matrix(mutate(select(d, BIR_0:SUL_0), Void = 0))

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
    Base <- paste0(Base,    'target += normal_lpdf(y[',i,']|    r[',r_sel[i],']/(1+',  R$N0_i[i],'*A[',R$Aii_2[i],']+', R$N0_j[i],'*A[',R$Aij_2[i],']),  sigma1 );\n'  )
    GQ <- paste0(GQ, 'log_lik[',i,'] = normal_lpdf(y[',i,']|    r[',r_sel[i],']/(1+',  R$N0_i[i],'*A[',R$Aii_2[i],']+', R$N0_j[i],'*A[',R$Aij_2[i],']),  sigma1 );\n',
                       'y_sim[',i,'] = normal_rng(              r[',r_sel[i],']/(1+',  R$N0_i[i],'*A[',R$Aii_2[i],']+', R$N0_j[i],'*A[',R$Aij_2[i],']),  sigma1 );\n')
    
  }else{
    Base <-paste0(Base,     'target += normal_lpdf(y[',i,']| r[',r_sel[i],']/(1+',  R$N0_i[i],'*A[',R$Aii_2[i],']),  sigma1 );\n'  )
    GQ <- paste0(GQ, 'log_lik[',i,'] = normal_lpdf(y[',i,']| r[',r_sel[i],']/(1+',  R$N0_i[i],'*A[',R$Aii_2[i],']),  sigma1 );\n' ,
                       'y_sim[',i,'] = normal_rng(           r[',r_sel[i],']/(1+',  R$N0_i[i],'*A[',R$Aii_2[i],']),  sigma1 );\n' )
  }
}

Base <-paste0(Base,'\n // Start of second treatment \n ')

############# Do second treatment ######
#######################################

d=PreparedData2

n=nrow(d)

N0s<-  as.matrix(mutate(select(d, BIR_0:SUL_0), Void = 0))

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
    Base <- paste0(Base, 'target += normal_lpdf(y[',i+n_firstset,']|    r[',r_sel[i],']/(1+',  R$N0_i[i],'*A[',R$Aii_2[i],']+', R$N0_j[i],'*A[',R$Aij_2[i],']),  sigma2 );\n'  )
    GQ <- paste0(GQ,   'log_lik[',i+n_firstset,']  = normal_lpdf(y[',i+n_firstset,']|    r[',r_sel[i],']/(1+',  R$N0_i[i],'*A[',R$Aii_2[i],']+', R$N0_j[i],'*A[',R$Aij_2[i],']),  sigma2 );\n',
                        "y_sim[", i+n_firstset, "] = normal_rng(                         r[",r_sel[i],']/(1+',  R$N0_i[i],'*A[',R$Aii_2[i],']+', R$N0_j[i],'*A[',R$Aij_2[i],']),  sigma2 );\n' )
  }else{
    Base <-paste0(Base,'target += normal_lpdf(y[',i+n_firstset,']| r[',r_sel[i],']/(1+',  R$N0_i[i],'*A[',R$Aii_2[i],']),  sigma2 );\n'  )
    GQ <- paste0(GQ,   'log_lik[',i+n_firstset,'] = normal_lpdf(y[',i+n_firstset,']|    r[',r_sel[i],']/(1+',  R$N0_i[i],'*A[',R$Aii_2[i],']),  sigma2 );\n' ,
                       "y_sim[", i+n_firstset, '] = normal_rng(                         r[',r_sel[i],']/(1+',  R$N0_i[i],'*A[',R$Aii_2[i],']),  sigma2 );\n')
  }
}

#######
### LOO needs log lik in generated quantitites. 
###### 


Both <- paste0(Base,'\n} \n', GQ, '\n}')


writeLines(Both, paste0('StanModels/BuiltModel_',Name,'.stan'))
return(paste0('Function Finished. Output: <BuiltModel_',Name,'.stan>'))
}
