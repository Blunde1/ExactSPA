// Copyright (C) 2018 Berent Lunde

#include "includes/ExactSPA.hpp"


// [[Rcpp::export]]
Rcpp::List nll_mjd(NumericVector X, double dt,
                      double r, double lsigma, double llambda, double mu, double lnu,
                      double length, double n)
{
    int nobs = X.size();

    // Warm up phase - solve inner problems and store solutions
    double sigma = exp(lsigma), lambda = exp(llambda), nu = exp(lnu);
    cgf_mjd<double> cgf(r, sigma, lambda, mu, nu, X[0], dt);
    NumericVector sp_warmup(nobs-1);
    for(int i=1; i<nobs; i++){
        cgf.set_x0(X[i-1]);
        sp_warmup[i-1] = sp_newton_opt(0.0, X[i], cgf);
    }

    // Build nll - start taping
    adept::Stack stack;
    adtype ad_r=r, ad_lsigma=lsigma, ad_llambda=llambda, ad_mu=mu, ad_lnu=lnu;

    adtype ad_x0=X[0], ad_dt=dt, ad_x=X[1], ad_lfx=0.0, ad_nll=0.0;

    stack.new_recording();

    adtype ad_sigma = exp(ad_lsigma), ad_lambda = exp(ad_llambda), ad_nu = exp(ad_lnu);

    cgf_mjd<adtype> ad_cgf(ad_r, ad_sigma, ad_lambda, ad_mu, ad_nu, ad_x0, ad_dt);

    adtype ad_s=0.0;

    for(int j=1; j<nobs; j++){
        // Update cgf
        ad_x0 = X[j-1];
        ad_x = X[j];
        ad_cgf.set_x0(ad_x0);

        // Solve inner problem
        ad_s = sp_warmup[j-1];
        ad_s = sp_newton_opt(ad_s, ad_x, ad_cgf, false, 2);

        // Create standardized CF
        lcf_stand_mjd<adtype> ad_lcf_stand( ad_cgf, ad_x, ad_s );

        // Calculate log exact spa
        ad_lfx = log_exact_spa(ad_x, ad_cgf, ad_s, ad_lcf_stand, length, n);

        ad_nll -= ad_lfx;

    }

    // Calculate gradient
    adtype res = ad_nll/1.0;
    res.set_gradient(1.0);
    stack.compute_adjoint();

    // Return nll and derivatives
    double r_grad = ad_r.get_gradient();
    double lsigma_grad = ad_lsigma.get_gradient();
    double llambda_grad = ad_llambda.get_gradient();
    double mu_grad = ad_mu.get_gradient();
    double lnu_grad = ad_lnu.get_gradient();

    return Rcpp::List::create(
        Named("r_grad")  = r_grad,
        Named("sigma_grad")  = lsigma_grad,
        Named("lambda_grad") = llambda_grad,
        Named("mu_grad") = mu_grad,
        Named("nu_grad") = lnu_grad,
        Named("nll") = res.value());
}


/*** R
x = seq(-0.1, .1, length.out=100)
y <- rep(0, length(x))
for(i in 1:length(x)){
    y[i] <- exp(-nll_mjd(c(0,x[i]), 1/250, .08, log(.1), log(100), -.001, log(.015), 12, 64)$nll)
}
plot(x,y)
    */

// [[Rcpp::export]]
Rcpp::List nll_nchisq(NumericVector X,
                      double ldf, double lncp,
                      double length, double n)
{
    int nobs = X.size();

    // Warm up phase - solve inner problems and store solutions
    double df = exp(ldf), ncp = exp(lncp);
    cgf_nchisq<double> cgf(df, ncp);
    NumericVector sp_warmup(nobs);
    for(int i=0; i<nobs; i++){
        sp_warmup[i] = sp_newton_opt(0.0, X[i], cgf);
        //cout << "saddlepoints iter" << i<< " value: "<< sp_warmup[i] << endl;
    }

    // Build nll - start taping
    adept::Stack stack;
    adtype ad_ldf = ldf, ad_lncp = lncp;

    adtype ad_lfx=0.0, ad_nll=0.0, ad_x = X[0], ad_s=0.0;

    stack.new_recording();

    adtype ad_df = exp(ad_ldf), ad_ncp = exp(ad_lncp);

    cgf_nchisq<adtype> ad_cgf(ad_df, ad_ncp);

    for(int j=0; j<nobs; j++){

        // solve inner problem
        ad_x = X[j];
        ad_s = sp_warmup[j];
        ad_s = sp_newton_opt(ad_s, ad_x, ad_cgf, false, 2);

        // Create standardized CF
        lcf_stand_nchisq<adtype> ad_lcf_stand( ad_cgf, ad_x, ad_s );

        // Calculate log exact spa
        ad_lfx = log_exact_spa(ad_x, ad_cgf, ad_s, ad_lcf_stand, length, n);

        // Update likelihood
        ad_nll -= ad_lfx;
    }

    // Calculate gradient
    adtype res = ad_nll/1.0;
    res.set_gradient(1.0);
    stack.compute_adjoint();

    // Return nll and derivatives
    double ldf_grad = ad_ldf.get_gradient();
    double lncp_grad = ad_lncp.get_gradient();

    return Rcpp::List::create(
        Named("ldf_grad") = ldf_grad,
        Named("lncp_grad") = lncp_grad,
        Named("nll") = res.value()
    );
}

/*** R
df = 100
ncp = 40
n = 1000
x <- sort(rchisq(n, df, ncp))
    hist(x, freq=F)
    curve(dchisq(x, df, ncp), add=TRUE)

    -sum(dchisq(x, df, ncp, log=TRUE))
    nll_nchisq(x, log(df), log(ncp), 20.0, 128)

    y <- rep(0, length(x))
    for(i in 1:length(x)){
        y[i] <- exp(-nll_nchisq(as.vector(x[i]),log(df),log(ncp), 12.0, 64)$nll)
    }
    lines(x,y,col="blue")
        */
