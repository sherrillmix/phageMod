data{
  int <lower=0> nObs;
  int <lower=0, upper=1> isPhage[nObs];
  int <lower=0> counts[nObs];
  real <lower=0> expected[nObs];
  //real <lower=0> isModified[nObs];
}

parameters{
  real betaIntercept;
  real betaPhage;
  //real betaMod;
  //real betaPhageMod;
  real <lower=0, upper=1> theta[2];
}

transformed parameters {
  real <lower=0> lambda[nObs];
  for(ii in 1:nObs){
    //+betaMod*isModified[ii]+betaPhageMod*isModified[ii]*isPhage[ii])
    lambda[ii]<-exp(betaIntercept+betaPhage*isPhage[ii]);
  }
}

model{
  //flat priors

  for(ii in 1:nObs){
    if(counts[ii]==0){
      increment_log_prob(log_sum_exp(bernoulli_log(1,theta[isPhage[ii]+1]^expected[ii]), bernoulli_log(0,theta[isPhage[ii]+1]^expected[ii])+ poisson_log(counts[ii],expected[ii]*lambda[ii])));
    }else{
      increment_log_prob(bernoulli_log(0,theta[isPhage[ii]+1]^expected[ii])+ poisson_log(counts[ii],expected[ii]*lambda[ii]));
    }
  }
}
