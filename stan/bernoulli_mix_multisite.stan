functions {
  // functions
  real bern_mix_lpmf(int x, real p, real q){
    real l;
    l = 0;
    
    // for(i in x){
    //   if (i == -1)
    //     // l += p * (1-q);
    //     l += log(p) + log(1-q);
    //     else if (i == 0)
    //       // l += 1 - p;
    //       l += log(1-p);
    //       else if (i == 1)
    //         // l += p * q;
    //         l += log(p) + log(q);
    // }
    if (x == -1)
      l += log(p) + log(1-q);
    else if (x == 0)
      l += log(1-p);
    else if (x == 1)
      l += log(p) + log(q);
  
    return l;
  }
}

data {
  int nobs;
  int<lower = -1, upper = 1> x[nobs]; // Whether there was a change and its direction
  int<lower = 1> nsites;
  int<lower = 1> id[nobs];
  real<lower = 0>vp; // Variance of beta distributions
  real<lower = 0>vq; // Variance of beta distributions
}

transformed data {
  // In case we need to transform
}

parameters {
  // Probably needed for the multisite version
  real<lower = 0, upper = 1>P; // Prob of change = P(x != 0)
  real<lower = 0, upper = 1>Q; // P(change | x != 0)
  
  real<lower = 0, upper = 1>p[nsites]; // Prob of change = P(x != 0)
  real<lower = 0, upper = 1>q[nsites]; // P(change | x != 0)
}

transformed parameters {
}

model {
  // Priors
  P ~ uniform(0,1);
  Q ~ uniform(0,1);
  
  p ~ beta(P * vp, (1-P) * vp);
  q ~ beta(Q * vq, (1-Q) * vq);
  
  for(i in 1:nobs){
    x[i] ~ bern_mix(p[id[i]], q[id[i]]);
  }
}
