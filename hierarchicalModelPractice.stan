data{
    int<lower=1> N; //number of observations
    int<lower=1> N_tax; //number of taxa
    //int<lower=1> N_comm; //number of communities
    int<lower=1,upper=N_tax> Taxon[N];//vector of species indices
    //int<lower=1,upper=N_comm> Community[N];//vector of community indices
    vector[N] Value;  //efficiency estimate for each taxon in each community
}
parameters{
    real<lower=0.00001> shape1; 
    real<lower=0.00001> shape2;
    vector<lower=0.00001, upper=0.9999>[N_tax] a_hat;  //each taxon will have a single mean efficiency with a single variance...
    real<lower=0> a_sigma;  //...shared among all taxa
}
model{
    a_hat ~ beta( shape1 , shape2 );

    shape1 ~ lognormal( 0 , 0.3 );
    shape2 ~ lognormal( 0 , 0.3 );

  for (i in 1:N) {
      Value[i] ~ normal(a_hat[Taxon[i]], a_sigma);
    }
    
}
generated quantities{
    real dev;
    dev = 0;
    dev = dev + (-2)*beta_lpdf( a_hat | shape1 , shape2 );
}
