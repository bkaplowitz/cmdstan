data {
    int<lower=1> T_j; // number of time periods
    int<lower=1> J_max; //maximum maturity J over all periods
    int<lower=1> M_max; // maximum number of assets in each period
    matrix<lower=0>[M_max,J_max] S[T_j]; // of assets sold each period //declare for each time period
    matrix<lower=0>[T_j,M_max] P; // value of P_observed in period t for maturity j // change to P
    }

parameters {
    matrix<lower=-1,upper=1> [3,T_j] a; //parameter values for base litterman-scheinkman model, constraints neeeded here for init. Highly constrained.... causes NA warnings, but makes sampler work.
    vector<lower=0>[M_max] sigma; //sd, fix covariance matrix later
    vector<lower=0>[T_j] sds; //sd, fix covariance matrix later
     //matrix<lower=0>[T_j,M_max] P_tilde;   

    }

transformed parameters {
    matrix[T_j,J_max] r_tilde; 
    matrix[T_j,J_max] r;
    matrix[3,T_j] a_tilde; //t by # of params
    for (t in 1:T_j) {
        for (k in 1:3) {
            a_tilde[k,t] = mean(a[k,1:t]); // can be changed to a weighted mean
            }  
        for (j in 1:J_max) {
            r[t,j]=-1*(a[1,t]+a[2,t]*j+a[3,t]*j^2)*j;
            }
        }
    }
model {
    sigma ~ gamma(0.001, 0.001);
    sds ~gamma(0.001, 0.001);
   //L ~ lkj_corr_cholesky(2);
    for (t in 1:T_j) {
        for (m in 1:M_max) {
            if (P[t,m]>0)  {
                P[t,m]~normal((S[t,m,1:J_max]*exp(r[t,1:J_max])'),sigma[m]); //slicing for efficiency.

                //print((S[t,m,1:J_max]*exp(r[t,1:J_max])'))
                //print((S[t,m,1:J_max]))
                //print(exp(r[t,1:J_max])')
                //P[t,m]~normal(P_tilde[t,m],sigma[m]); //slicing for efficiency.    
                }
            }
            for (k in 1:3) {
                a[k,t]~normal(a_tilde[k,t],sds[t]);  //can add correlations later if  we want, including between different T's and different a values.
            }
        }
    }
generated quantities {
 // matrix[M_max,M_max] Omega;
 // matrix[M_max,M_max] Sigma;
 // Omega = multiply_lower_tri_self_transpose(L); //correlation matrix
 // Sigma = quad_form_diag(Omega, sigma); //variance matrix
 // do difference of P and P_tilde being normal (0,sigma)
 //do alphas ~n(0,sigma), where sigma has covariance 1/t_diff
    }
