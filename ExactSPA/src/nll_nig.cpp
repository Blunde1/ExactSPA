// Copyright (C) 2018 Berent Lunde

#include "includes/ExactSPA.hpp"
//#include <conio.h>
//#include <unistd.h>
#include <iostream>

template<typename T>
struct cgf_nig{

    // Members
    T chi, psi, mu, gamma;

    // Constructors
    cgf_nig(T chi_, T psi_, T mu_, T gamma_) : chi(chi_), psi(psi_), mu(mu_), gamma(gamma_) {}
    cgf_nig(const cgf_nig& cgf_) : chi(cgf_.chi), psi(cgf_.psi), mu(cgf_.mu), gamma(cgf_.gamma) {}

    // Operators
    T operator()(T s){
        return s*mu + sqrt(chi*psi) - sqrt(chi*psi - chi*s*s - 2.0*s*chi*gamma);
    }
    cType<T> operator()(cType<T> s){
        cType<T> res;
        res = s*mu + cType<T>(sqrt(chi*psi)) - sqrt(cType<T>(chi*psi) - chi*s*s - T(2.0)*s*chi*gamma);
        return res;
    }


    // Member functions
    T deriv(T s){
        return mu - (sqrt(chi)*(-2.0*(gamma + s))) / (2.0*sqrt(psi-2.0*gamma*s-s*s));
    }
    T dderiv(T s){
        return sqrt(chi)*(psi+gamma*gamma)/( pow((psi-s*(2.0*gamma+s)),1.5) );
    }

};

// Line search Newton method
template<typename T>
T sp_newton_opt(T s_start, T x, cgf_nig<T>& cgf, bool check_values=true, int maxniter=1000){
    int i;
    double alpha = 1.0;
    T EPS = 1.0e-12, f_eps = 1.0e50;
    T g, h, s_old=s_start, s_new=s_start;
    for(i=0; i<maxniter; i++){
        g = cgf.deriv(s_old) - x;
        h = cgf.dderiv(s_old);
        s_new = s_old -  alpha * g/h;
        //cout << i << " " << alpha << endl;
        //cout << s_new << " " << g << " " << h << cgf(s_new) << endl;
        //sleep(1);

        if(check_values){
            if(isnan(cgf(s_new))){
                alpha = 0.5*alpha;
                s_old = s_start;
                i=0;
            }
            else if(abs(s_new - s_old) <= EPS){
                break;
            } else if(cgf(s_new) <= (cgf(s_old) + f_eps) ){
                s_old = s_new;
            } else{
                alpha = 0.5*alpha;
                s_old = s_start;
                i=0;
            }
        } else{
            s_old=s_new;
        }
    }
    //cout << "Iterations: " << i << endl;
    return s_new;
}


/* standardized log characteristic function */
template<typename T>
struct lcf_stand_nig{
    // Members
    cgf_nig<T>& cgf;
    T x, sp;

    // Constructor
    lcf_stand_nig(cgf_nig<T> &cgf_, T x_, T sp_) : cgf(cgf_), x(x_), sp(sp_) {}

    // Operator
    cType<T> operator()(cType<T> s){
        cType<T> res, i(0,1);
        T h = cgf.dderiv(sp);
        res = cType<T>(-cgf(sp)) - i*s*x / cType<T>(sqrt(h)) +
            cgf(  sp + s*i / cType<T>(sqrt(h)) );

        return res;
    }
};

/* Log exact spa */
template<typename T>
T log_exact_spa(T x, cgf_nig<T> cgf, T sp, lcf_stand_nig<T> lcf_stand, double length, int n){

    T lfx, fx_standardized;
    fx_standardized = ift_simpson0(lcf_stand, length, n);
    lfx = cgf(sp) - sp*x + log(fx_standardized) -
        0.5*log(cgf.dderiv(sp));

    return lfx;
}

/** Simpson's rule for IFT
 */
template<typename T>
T ift_simpson0(lcf_stand_nig<T> lcf_stand, double length, int n){
    T simpson_integral = 0.0;
    T a = 0.0, b = length;
    T h = (b-a) / n;

    simpson_integral = real(exp(lcf_stand(cType<T>(a)))) + real(exp(lcf_stand(cType<T>(b))));

    for(int j=1; j<=(n/2); j++){
        simpson_integral += 4.0 * real(exp(lcf_stand(cType<T>(h*(2*j-1)))));
    }
    for(int j=1; j<=(n/2-1); j++){
        simpson_integral += 2.0 * real(exp(lcf_stand(cType<T>(h*(2*j)))));
    }
    simpson_integral = simpson_integral*h/(3.0*M_PI);

    return simpson_integral;
}

// [[Rcpp::export]]
Rcpp::List nll_nig(NumericVector X,
                   double lchi, double lpsi, double mu, double gamma,
                   double length, double n)
{
    int nobs = X.size();

    // Warm up phase - solve inner problems and store solutions
    double chi = exp(lchi), psi = exp(lpsi);
    cgf_nig<double> cgf(chi, psi, mu, gamma);
    NumericVector sp_warmup(nobs);
    for(int i=0; i<nobs; i++){
        sp_warmup[i] = sp_newton_opt(0.0, X[i], cgf);
        //cout << "saddlepoints iter" << i<< " value: "<< sp_warmup[i] << endl;
    }

    // Build nll - start taping
    adept::Stack stack;
    adtype ad_lchi=lchi, ad_lpsi=lpsi, ad_mu=mu, ad_gamma=gamma;
    adtype ad_lfx=0.0, ad_nll=0.0, ad_x=X[0], ad_s=0.0;
    stack.new_recording();

    adtype ad_chi=exp(ad_lchi), ad_psi=exp(ad_lpsi);
    cgf_nig<adtype> ad_cgf(ad_chi, ad_psi, ad_mu, ad_gamma);

    for(int j=0; j<nobs; j++){

        // solve inner problem
        ad_x = X[j];
        ad_s = sp_warmup[j];
        ad_s = sp_newton_opt(ad_s, ad_x, ad_cgf, false, 2);

        // Create standardized CF
        lcf_stand_nig<adtype> ad_lcf_stand( ad_cgf, ad_x, ad_s );

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
    double lchi_grad = ad_lchi.get_gradient();
    double lpsi_grad = ad_lpsi.get_gradient();
    double mu_grad = ad_mu.get_gradient();
    double gamma_grad = ad_gamma.get_gradient();

    return Rcpp::List::create(
        Named("lchi_grad") = lchi_grad,
        Named("lpsi_grad") = lpsi_grad,
        Named("mu_grad") = mu_grad,
        Named("gamma_grad") = gamma_grad,
        Named("nll") = res.value()
    );
}

/*** R
chi = 3.0e-4; psi=1e3; mu=-3e-4; gamma=2
n = 1000
require(SuppDists)
rNIG=function(n,param){
    # param[1] : chi
    # param[2] : psi
    # param[3] : mu
    # param[4] : gamma
    #Alternative parametrization used in SuppDists
    nu<-sqrt(param[1]/param[2])
    #Simulate n inverse gaussian observations
    IGsim<-rinvGauss(n,nu,param[1])
    #Constructing n NIG observations based upon IGsim
    NIGsim<-param[3]+param[4]*IGsim+sqrt(IGsim)*rnorm(n)
    #returning the simulated NIG observations
    return(NIGsim)
}
x <- sort(rNIG(n, c(chi, psi, mu, gamma)))

loglikNIG=function(param,data){
    #This routine requires the observations to be stored in
    # a global variable "data"
    param[1]<-exp(param[1]) #Reparametrization
    param[2]<-exp(param[2])
    # param[1] : chi
    # param[2] : psi
    # param[3] : mu
    # param[4] : gamma
    # The square root expression
    y = sqrt((param[1] + (data - param[3])^2)*(param[2] + param[4]^2))
    # The log-likelihood
    loglik =
        sum(
            0.5*log(param[1]) +
                log(param[2] + param[4]^2) +
                sqrt(param[1]*param[2]) -
                log(pi) +
                log(besselK(y, -1, expon.scaled = TRUE)) - y +
                (data - param[3])*param[4] -
                log(y)
        )
    # Return the functional value
    return(loglik)
}

#The log-likelihood for a single observation=1 and all parameters=1 should read -0.9558895
nll_nig(1, 0, 0, 1, 1, 20, 256)
#nll_nig(x[50], log(chi), log(psi), mu, gamma, 20, 128)
exp(-nll_nig(mean(x), log(chi), log(psi), mu, gamma, 20, 128)$nll)
exp(-nll_nig(-0.1, log(chi), log(psi), mu, gamma, 20, 128)$nll)

nll_nig(x, log(chi), log(psi), mu, gamma, 100, 512)
loglikNIG(c(log(chi), log(psi), mu, gamma), x)
y <- y2 <- y3 <- rep(0, length(x))
for(i in 1:length(x)){
    y[i] <- exp(-nll_nig(x[i], log(chi), log(psi), mu, gamma, 50, 128)$nll)
    y2[i] <- exp(-nll_nig(x[i], log(chi), log(psi), mu, gamma, 8, 128)$nll)
    y3[i] <- exp(loglikNIG(c(log(chi), log(psi), mu, gamma), x[i]))
}
hist(x, freq=F, ylim=c(0,max(y)))
lines(x,y,col="red")
lines(x,y2,col="blue")
lines(x,y3, col="green")
    */
