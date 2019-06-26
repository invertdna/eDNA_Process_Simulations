data{
    int<lower=1> N; //number of observations
    int<lower=1> N_tax; //number of taxon-dataset combinations
    int<lower=1> N_datasets; //number of datasets
    //int<lower=1> N_comm; //number of communities
    int<lower=1,upper=N_tax> Taxon[N];//vector of species indices
    int<lower=1,upper=N_tax> Dataset[N];//vector of dataset indices
    //int<lower=1> N_a_hat;//vector of mean amp efficiencies  nrow(unique(dataForStan[, c("Taxon", "Dataset")]))
    //int<lower=1,upper=N_comm> Community[N];//vector of community indices
    vector[N] Value;  //efficiency estimate for each taxon in each community
}
parameters{
    vector<lower=0.00001>[N_datasets] shape1; 
    vector<lower=0.00001>[N_datasets] shape2;
    vector<lower=0.00001, upper=0.9999>[N_tax] a_hat;  //each taxon will have a single mean efficiency with a single variance...
    vector<lower=0>[N_datasets] a_sigma;  //...shared among all taxa within a dataset
}
model{
  
    shape1 ~ lognormal( 0 , 0.3 );
    shape2 ~ lognormal( 0 , 0.3 );

  for (i in 1:N) {
      Value[i] ~ normal(a_hat[Taxon[i]], a_sigma[Dataset[i]]);
      a_hat[Taxon[i]] ~ beta( shape1[Dataset[i]] , shape2[Dataset[i]] );
  }
    
}
// generated quantities{
//     real dev;
//     dev = 0;
//     dev = dev + (-2)*beta_lpdf( a_hat | shape1 , shape2 );
// }
