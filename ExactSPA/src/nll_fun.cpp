// Copyright (C) 2018 Berent Lunde

#include "includes/ExactSPA.hpp"
//#include "C:/Users/Berent/Projects/it-ift/implementation v5/ExactSPA/ExactSPA/src/includes/ExactSPA.hpp"


// [[Rcpp::export]]
Rcpp::List nll_mjd(NumericVector X, double dt,
                      double r, double lsigma, double llambda, double mu, double lnu,
                      double length, double n,
                      int ift_type)
{
    int nobs = X.size();

    // Warm up phase - solve inner problems and store solutions
    double sigma = exp(lsigma), lambda = exp(llambda), nu = exp(lnu);
    cgf_mjd<double> cgf(r, sigma, lambda, mu, nu, X[0], dt);
    NumericVector sp_warmup(nobs-1);

    switch(ift_type)
    {
    case 1: case 2:
        for(int i=1; i<nobs; i++){
            cgf.set_x0(X[i-1]);
            sp_warmup[i-1] = sp_newton_opt(0.0, X[i], cgf);
        }
        break;
    }

    // Build nll - start taping
    adept::Stack stack;
    adtype ad_r=r, ad_lsigma=lsigma, ad_llambda=llambda, ad_mu=mu, ad_lnu=lnu;

    adtype ad_x0=X[0], ad_dt=dt, ad_x=X[1], ad_lfx=0.0, ad_nll=0.0, ad_fx=0.0, ad_w=0.0;

    stack.new_recording();

    adtype ad_sigma = exp(ad_lsigma), ad_lambda = exp(ad_llambda), ad_nu = exp(ad_lnu);
    adtype ad_k = exp(ad_mu + 0.5*ad_nu*ad_nu) - 1.0;
    adtype ad_intensity = 0.0, ad_mean = 0.0, ad_SD = 0.0;

    cgf_mjd<adtype> ad_cgf(ad_r, ad_sigma, ad_lambda, ad_mu, ad_nu, ad_x0, ad_dt);

    adtype ad_s=0.0;

    for(int j=1; j<nobs; j++){
        // Update cgf
        ad_x0 = X[j-1];
        ad_x = X[j];
        ad_cgf.set_x0(ad_x0);

        switch(ift_type)
        {
        case 1: case 2: // SPA versions
            // Solve inner problem
            ad_s = sp_warmup[j-1];
            ad_s = sp_newton_opt(ad_s, ad_x, ad_cgf, false, 2);

            switch(ift_type)
            {
            case 1: /* Ordinary spa */
                ad_lfx = log_standard_spa(ad_x, ad_cgf, ad_s);
                break;

            case 2: /* Saddlepoint adjusted ift */
                // Create standardized CF
                lcf_stand_mjd<adtype> ad_lcf_stand( ad_cgf, ad_x, ad_s );
                // Calculate log exact spa
                ad_lfx = log_exact_spa(ad_x, ad_cgf, ad_s, ad_lcf_stand, length, n);
                break;
            }
            break;
        case 3: // Direct IFT - simpson
            ad_fx = ift_simpson(ad_cgf, ad_x, length, n);
            if(ad_fx > 1.0e-100){
                ad_lfx = log(ad_fx);
            }else{
                ad_lfx = log(1.0e-14);
            }
            break;
        case 4: // Truncated
            int maxjumps = 20;
            ad_fx = 0.0;
            for(int i=0; i<maxjumps; i++){
                ad_intensity = ad_lambda * ad_dt;
                ad_w = pois_dens(i, ad_intensity); // maybe will not work
                if(ad_w < 1e-14){
                    break;
                    }
                ad_mean = ad_x0 + (ad_r - 0.5*ad_sigma*ad_sigma - ad_lambda*ad_k)*ad_dt + i*mu;
                ad_SD = sqrt(ad_dt*ad_sigma*ad_sigma + i*ad_nu*ad_nu);
                ad_fx += ad_w * norm_dens(ad_x, ad_mean, ad_SD);
                }
            ad_lfx = log(ad_fx);
            break;
        }

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
if(FALSE){

    x = seq(-0.1, .1, length.out=100)
    y <- y2 <- y3 <- y4 <- rep(0, length(x))
    for(i in 1:length(x)){
        y[i] <- exp(-nll_mjd(c(0,x[i]), 1/250, .08, log(.1), log(100), -.001, log(.015), 12, 64, 2)$nll)
        y2[i] <- exp(-nll_mjd(c(0,x[i]), 1/250, .08, log(.1), log(100), -.001, log(.015), 12, 64, 1)$nll)
        y3[i] <- exp(-nll_mjd(c(0,x[i]), 1/250, .08, log(.1), log(100), -.001, log(.015), 700, 64, 3)$nll)
        y4[i] <- exp(-nll_mjd(c(0,x[i]), 1/250, .08, log(.1), log(100), -.001, log(.015), 700, 64, 4)$nll)
    }
    plot(x,y)
    lines(x,y2, col="blue")
    lines(x,y3, col="green")
    lines(x,y4, col="yellow")
    plot(x,log(y), type="b")
    lines(x,log(y3),col="red")

}
*/

// [[Rcpp::export]]
Rcpp::List nll_nchisq(NumericVector X,
                      double ldf, double lncp,
                      double length, double n,
                      int ift_type)
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

        switch(ift_type){
        case 1: // Ordinary spa
            ad_lfx = log_standard_spa(ad_x, ad_cgf, ad_s);
            break;

        case 2: // Saddlepoint adjusted ift
            // Create standardized CF
            lcf_stand_nchisq<adtype> ad_lcf_stand( ad_cgf, ad_x, ad_s );

            // Calculate log exact spa
            ad_lfx = log_exact_spa(ad_x, ad_cgf, ad_s, ad_lcf_stand, length, n);
            break;
        }

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
if(FALSE){
    df = 100
    ncp = 40
    n = 1000
    x <- sort(rchisq(n, df, ncp))
    hist(x, freq=F)
    curve(dchisq(x, df, ncp), add=TRUE)

    -sum(dchisq(x, df, ncp, log=TRUE))
    nll_nchisq(x, log(df), log(ncp), 20.0, 128, 2)
    nll_nchisq(x, log(df), log(ncp), 20.0, 128, 1)

    y <- y2 <- rep(0, length(x))
    for(i in 1:length(x)){
        y[i] <- exp(-nll_nchisq(as.vector(x[i]),log(df),log(ncp), 12.0, 64, 2)$nll)
        y2[i] <- exp(-nll_nchisq(as.vector(x[i]),log(df),log(ncp), 12.0, 64, 1)$nll)
    }
    lines(x,y,col="blue")
    lines(x,y2, col="green")
}
*/



// [[Rcpp::export]]
Rcpp::List nll_nig(NumericVector X,
                   double lchi, double lpsi, double mu, double gamma,
                   double length, double n,
                   int ift_type)
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
    adtype ad_fx = 0.0;
    stack.new_recording();

    adtype ad_chi=exp(ad_lchi), ad_psi=exp(ad_lpsi);
    cgf_nig<adtype> ad_cgf(ad_chi, ad_psi, ad_mu, ad_gamma);

    for(int j=0; j<nobs; j++){

        ad_x = X[j];
        // solve inner problem
        switch(ift_type){
        case 1: case 2: // SPA versions
            // Solve inner problem
            ad_s = sp_warmup[j];
            ad_s = sp_newton_opt(ad_s, ad_x, ad_cgf, false, 2);

            switch(ift_type){
            case 1: // Ordinary spa
                ad_lfx = log_standard_spa(ad_x, ad_cgf, ad_s);
                break;

            case 2: // Saddlepoint adjusted ift
                // Create standardized CF
                lcf_stand_nig<adtype> ad_lcf_stand( ad_cgf, ad_x, ad_s );

                // Calculate log exact spa
                ad_lfx = log_exact_spa(ad_x, ad_cgf, ad_s, ad_lcf_stand, length, n);
                break;
            }
            break;
        case 3: // Direct IFT - simpson
            ad_fx = ift_simpson(ad_cgf, ad_x, length, n);
            if(ad_fx > 1.0e-100){
                ad_lfx = log(ad_fx);
            }else{
                ad_lfx = log(1.0e-14);
            }
            // ad_lfx = max()
            // ad_lfx = log(ift_simpson(ad_cgf, ad_x, length, n));
            // if(!isfinite(ad_lfx)){
            //     ad_lfx = log(1.0e-100);
            // }
        }

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
if(FALSE){

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
        nll_nig(1, 0, 0, 1, 1, 20, 256, 2)
        nll_nig(1, 0, 0, 1, 1, 20, 256, 1)
        nll_nig(1, 0, 0, 1, 1, 20, 256, 3)
    #nll_nig(x[50], log(chi), log(psi), mu, gamma, 20, 128)
            exp(-nll_nig(mean(x), log(chi), log(psi), mu, gamma, 20, 128, 2)$nll)
            exp(-nll_nig(-0.1, log(chi), log(psi), mu, gamma, 20, 128, 2)$nll)

            nll_nig(x, log(chi), log(psi), mu, gamma, 100, 512, 2)
            nll_nig(x, log(chi), log(psi), mu, gamma, 1500, 252, 3)
            loglikNIG(c(log(chi), log(psi), mu, gamma), x)

            y <- y2 <- y3 <- y4 <- y5 <- rep(0, length(x))
            for(i in 1:length(x)){
                y[i] <- exp(-nll_nig(x[i], log(chi), log(psi), mu, gamma, 50, 128, 2)$nll)
                y2[i] <- exp(-nll_nig(x[i], log(chi), log(psi), mu, gamma, 8, 128, 2)$nll)
                y3[i] <- exp(loglikNIG(c(log(chi), log(psi), mu, gamma), x[i]))
                y4[i] <- exp(-nll_nig(x[i], log(chi), log(psi), mu, gamma, 8, 128, 1)$nll)
                y5[i] <- exp(-nll_nig(x[i], log(chi), log(psi), mu, gamma, 500, 1048, 3)$nll)
            }
            hist(x, freq=F, ylim=c(0,max(y)))
                lines(x,y,col="red")
                lines(x,y2,col="blue")
                lines(x,y3, col="green")
                lines(x,y4, col="black", type="l", lty=4)
                lines(x, y5, lwd=2, type="l", col="red", lty=3)

    # testing direct simpson ift
    y5 <- rep(0, length(x))
    for(i in 1:length(x)){
        y5[i] <- exp(-nll_nig(x[i], log(chi), log(psi), mu, gamma, 500, 1048, 3)$nll)
    }
lines(x, y5, lwd=3, type="b", col="red")
}
*/


    // [[Rcpp::export]]
    Rcpp::List ldens_tweedie(double x,
                             double t_xi, double lphi, double lmu,
                             double length, double n,
                             int ift_type)
    {

        // Warm up phase - solve inner problems and store solutions
        double xi = exp(t_xi) + 2, phi=exp(lphi), mu=exp(lmu);

        // Build nll - start taping
        adept::Stack stack;
        adtype ad_t_xi = t_xi, ad_lphi=lphi, ad_lmu = lmu;

        adtype ad_lfx=0.0, ad_x = x, ad_s=0;

        stack.new_recording();

        adtype ad_xi = exp(ad_t_xi) + 2, ad_phi = exp(ad_lphi), ad_mu = exp(ad_lmu);

        cgf_tweedie<adtype> ad_cgf(ad_xi, ad_mu, ad_phi);

        switch(ift_type){

        case 1: case 2: // SPA versions
        {
            // Solve inner problem
            // not taped
            cgf_tweedie<double> cgf(xi, mu, phi);
            ad_s = sp_newton_opt(0.0, x, cgf); //sp;
            // taped
            ad_s = sp_newton_opt(ad_s, ad_x, ad_cgf, false, 2);

            // spa or sp adjusted ift
            switch(ift_type){
            case 1: // ordinary spa
                ad_lfx = log_standard_spa(ad_x, ad_cgf, ad_s);
                break;

            case 2: // Saddlepoint adjusted ift

                lcf_stand_tweedie<adtype> ad_lcf_stand(ad_cgf, ad_x, ad_s);
                ad_lfx = log_exact_spa(ad_x, ad_cgf, ad_s, ad_lcf_stand, length, n);
                break;
            }
        }
            break;
        case 3: // Direct IFT - simpson
            ad_lfx = log(ift_simpson(ad_cgf, ad_x, length, n));

        }

        // Evaluate stack
        adtype res = ad_lfx/1.0;
        res.set_gradient(1.0);
        stack.compute_adjoint();

        // Calculate gradient
        double t_xi_grad = ad_t_xi.get_gradient();
        double lphi_grad = ad_lphi.get_gradient();
        double lmu_grad = ad_lmu.get_gradient();

        // Return lfx and corresponding derivatives w.r.t. param
        return Rcpp::List::create(
            Named("t_xi_grad") = t_xi_grad,
            Named("lphi_grad") = lphi_grad,
            Named("lmu_grad") = lmu_grad,
            Named("lfx") = res.value()
        );
    }

/*** R
if(FALSE){

    # Parameters and comparison library
    library(tweedie)
    p <- 4
    mu <- 2
    phi <- 2
    y <- seq(0, 30, length=100)
    fy <- log(dtweedie( y=y, power=p, mu=mu, phi=phi))

    # compare to my spa
    f.spa <- sapply(y, function(x){ldens_tweedie(x,t_xi = log(p-2),
                                                     lphi=log(phi),
                                                     lmu=log(mu),
                                                     length=40,
                                                     n=100,
                                                     ift_type=1)$lfx})


    # Compare to the saddlepoint density
    f.saddle <- log(dtweedie.saddle( y=y, power=p, mu=mu, phi=phi))

    # compare to spaift
    f.spaift <- sapply(y, function(x){ldens_tweedie(x,t_xi = log(p-2),
                                                        lphi=log(phi),
                                                        lmu=log(mu),
                                                        length=1000,
                                                        n=500,
                                                        ift_type=2)$lfx})

    # Plot log density and compare
    plot(y, fy, type="l")#, ylim=range(fy, f.spa, f.saddle, f.spaift, na.rm=T))
    points(y, f.spa, col=2, type = "p", pch=2)
    lines( y, f.saddle, col=3 )
    lines(y, f.spaift, col=4, type = "p", pch=3)

    legend("topright",
           c("R tweedie","R tweedie.saddle", "Saddlepoint approximation", "Saddlepoint adjusted IFT"),
           lty=c(1,1,NA,NA), col=c(1,3,2,4), pch=c(NA,NA,2,3))


    # Plot density and compare
    plot(y, exp(fy), type="l")#, ylim=range(fy, f.spa, f.saddle, f.spaift, na.rm=T))
    points(y, exp(f.spa), col=2, type = "p", pch=2)
    lines( y, exp(f.saddle), col=3 )
    lines(y, exp(f.spaift), col=4, type = "p", pch=3)

    legend("topright",
           c("R tweedie","R tweedie.saddle", "Saddlepoint approximation", "Saddlepoint adjusted IFT"),
           lty=c(1,1,NA,NA), col=c(1,3,2,4), pch=c(NA,NA,2,3))




    # Check tails
    y <- seq(0,50, length=100)
    lfy <- log(dtweedie( y=y, power=p, mu=mu, phi=phi))
    lfy.spaift <- sapply(y, function(x){ldens_tweedie(x,t_xi = log(p-2),
                                                          lphi=log(phi),
                                                          lmu=log(mu),
                                                          length=1000,
                                                          n=500,
                                                          ift_type=2)$lfx})
    # compare to my spa
    lf.spa <- sapply(y, function(x){ldens_tweedie(x,t_xi = log(p-2),
                                                     lphi=log(phi),
                                                     lmu=log(mu),
                                                     length=40,
                                                     n=100,
                                                     ift_type=1)$lfx})
    plot(y, lfy, type="l" )
    points(y, lfy.spaift, type="p", pch=2, col=2)
    points(y, lf.spa, type = "p", pch=3,col=3)

    legend("topright",
           c("R tweedie","R tweedie.saddle", "Saddlepoint approximation", "Saddlepoint adjusted IFT"),
           lty=c(1,1,NA,NA), col=c(1,3,2,4), pch=c(NA,NA,2,3))

}
*/
