data{
    int<lower=1> N;
    real A[N];
}
parameters{
    real shape1;
    real shape2;
}
model{
    shape2 ~ lognormal( 0 , 0.3 );
    shape1 ~ lognormal( 0 , 0.3 );
    A ~ beta( shape1 , shape2 );
}
generated quantities{
    real dev;
    dev = 0;
    dev = dev + (-2)*beta_lpdf( A | shape1 , shape2 );
}
