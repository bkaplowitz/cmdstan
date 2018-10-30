data {
    int<lower=1> T_j; // number of time periods
    int<lower=1> J_max; //maximum maturity J over all periods
    int<lower=1> M_max; // maximum number of assets in each period
    matrix[T_j,J_max] S[M_max]; // of assets sold each period //declare for each time period
    matrix<lower=0>[T_j,M_max] P; // value of P_observed in period t for maturity j // change to P
    }

parameters {
    matrix [3,T_j] a; //parameter values for base litterman-scheinkman model, constraints neeeded here for init. Highly constrained.... causes NA warnings, but makes sampler work.
    vector<lower=0>[M_max] sigma; //sd, fix covariance matrix later
    vector<lower=0>[3] sds; //sd, fix covariance matrix later
    }

transformed parameters {
    matrix[T_j,J_max] r; // fitted yield to maturity
    matrix[3,T_j] a_tilde; // # of params by t
    for (t in 1:T_j) {
        for (j in 1:J_max) {
            r[t,j]=-1*(a[1,t]+a[2,t]*j+a[3,t]*j^2)*j;
            }
        for (k in 1:3) {
            a_tilde[k,t] = mean(a[k,1:t]); // can be changed to a weighted mean
            }  
        }
    }
model {
    sigma ~ inv_gamma(0.001, 0.001);
    sds ~inv_gamma(0.001, 0.001);
    for (m in 1:M_max) {
        for (t in 1:T_j) {
            if (P[t,m]>0)  {
                a[1:3,t]~normal(a_tilde[1:3,t],sds[1:3]);  //can add correlations later if  we want, including between different T's and different a values.
                P[t,m]~normal((S[m,t,1:J_max]*exp(r[t,1:J_max])'),sigma[m]); //slicing for efficiency.
                }
            }            
        }
    }
