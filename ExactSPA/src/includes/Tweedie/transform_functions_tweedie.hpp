// Copyright (C) 2016-2018 Berent Lunde

#ifndef __transform_functions_tweedie_hpp_included__
#define __transform_functions_tweedie_hpp_included__

#include "../complex.hpp"

template<typename T>
struct cgf_tweedie{

    // Members
    T xi, mu, phi;

    // Constructors
    cgf_tweedie(T xi_, T mu_, T phi_) : xi(xi_), mu(mu_), phi(phi_) {}
    cgf_tweedie(const cgf_tweedie& cgf_) : xi(cgf_.xi), mu(cgf_.mu), phi(cgf_.phi) {}

    // Operators
    T operator()(T s){

        //mu^(2-xi)/(phi*(2-xi)) * ((1-s*phi*(xi-1)*mu^(xi-1))^((2-xi)/(1-xi)) - 1)

        return pow(mu, 2-xi)/(phi*(2-xi)) * ( pow(1-s*phi*(xi-1)*pow(mu,xi-1), ((2-xi)/(1-xi)))  -1);
    }
    cType<T> operator()(cType<T> s){
        cType<T> res;
        res = cType<T>(pow(mu, 2-xi)/(phi*(2-xi))) *
            ( pow(T(1.0)-s*phi*cType<T>(xi-T(1))*cType<T>(pow(mu,xi-1)), cType<T>((T(2)-xi)/(T(1)-xi)))  -T(1) );
        //res = ncp*s/(T(1.0)-T(2.0)*s) - cType<T>(df/T(2.0))*log(T(1.0)-T(2.0)*s);
        return res;
    }


    // Member functions
    T deriv(T s){
        //     mu * (1-s*phi*(xi-1)*mu^(xi-1))^(1/(1-xi))

        return mu * pow(1-s*phi*(xi-1)*pow(mu, xi-1) , 1/(1-xi) );
    }
    T dderiv(T s){
        //     phi*mu^xi * (1-s*phi*(xi-1)*mu^(xi-1))^(xi/(1-xi))

        return phi*pow(mu,xi) * pow(1-s*phi*(xi-1)*pow(mu, xi-1) , xi/(1-xi) );
    }

};


/* standardized log characteristic function */
template<typename T>
struct lcf_stand_tweedie{
    // Members
    cgf_tweedie<T>& cgf;
    T x, sp;

    // Constructor
    lcf_stand_tweedie(cgf_tweedie<T> &cgf_, T x_, T sp_) : cgf(cgf_), x(x_), sp(sp_) {}

    // Operator
    cType<T> operator()(cType<T> s){
        cType<T> res, i(0,1);
        T h = cgf.dderiv(sp);
        res = cType<T>(-cgf(sp)) - i*s*x / cType<T>(sqrt(h)) +
            cgf(  sp + s*i / cType<T>(sqrt(h)) );

        return res;
    }
};


#endif //__transform_functions_tweedie_hpp_included__
