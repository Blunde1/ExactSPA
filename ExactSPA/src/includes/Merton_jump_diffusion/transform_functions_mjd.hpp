// Copyright (C) 2016-2018 Berent Lunde

#ifndef __transform_functions_mjd_hpp_included__
#define __transform_functions_mjd_hpp_included__

#include "../complex.hpp"

template<typename T>
struct cgf_mjd{

    // Members
    T r, sigma, lambda, mu, nu, x0, dt;

    // Constructors
    cgf_mjd(T r_, T sigma_, T lambda_, T mu_, T nu_, T x0_, T dt_) :
        r(r_), sigma(sigma_), lambda(lambda_), mu(mu_), nu(nu_), x0(x0_), dt(dt_) {}
    cgf_mjd(const cgf_mjd& cgf_) :
        r(cgf_.r), sigma(cgf_.sigma), lambda(cgf_.lambda),
        mu(cgf_.mu), nu(cgf_.nu), x0(cgf_.x0), dt(cgf_.dt){}


    // Set methods
    void set_x0(T x0) { this -> x0 = x0; }

    // Operator
    T operator()(T s){

        T k = exp(T(mu)+0.5*T(nu)*T(nu)) - 1.0;
        T cgf_val = s*x0 + s*dt*(r-lambda*k- 0.5*sigma*sigma ) +
            0.5*s*s*sigma*sigma*dt +
            lambda*dt*(exp(s*mu+0.5*s*s*nu*nu) - 1.0);
        return cgf_val;
    }
    cType<T> operator()(cType<T> s){
        T k = exp(T(mu)+0.5*T(nu)*T(nu)) - 1.0;
        cType<T> res, i(0,1);
        res = s*x0 + s*dt * cType<T>(r-lambda*k- T(0.5)*sigma*sigma) +
            T(0.5)*s*s*sigma*sigma*dt +
            cType<T>(lambda*dt)*(exp(s*mu + T(0.5)*s*s*nu*nu) - T(1.0));
        return res;
    }

    // Member functions
    T deriv(T s){
        T k = exp(T(mu)+0.5*T(nu)*T(nu)) - 1.0;
        T cgf_deriv = x0 + dt*(r-lambda*k-sigma*sigma*0.5) +
            s*sigma*sigma*dt +
            (mu+nu*nu*s) * lambda*dt*exp(s*mu+s*s*nu*nu*0.5);
        return cgf_deriv;
    }
    T dderiv(T s){
        T cgf_dderiv = sigma*sigma*dt + lambda*dt*exp(s*mu+s*s*nu*nu*0.5)*
            (nu*nu + (mu+nu*nu*s)*(mu+nu*nu*s));
        return cgf_dderiv;
    }

};

/* standardized log characteristic function */
template<typename T>
struct lcf_stand_mjd{
    // Members
    cgf_mjd<T>& cgf;
    T x, sp;

    // Constructor
    lcf_stand_mjd(cgf_mjd<T> &cgf_, T x_, T sp_) : cgf(cgf_), x(x_), sp(sp_) {}

    // Operator
    cType<T> operator()(cType<T> s){
        cType<T> res, i(0,1);
        T h = cgf.dderiv(sp);
        res = cType<T>(-cgf(sp)) - i*s*x / cType<T>(sqrt(h)) +
            cgf(  sp + s*i / cType<T>(sqrt(h)) );

        return res;
    }
};


#endif //__transform_functions_mjd_hpp_included__
