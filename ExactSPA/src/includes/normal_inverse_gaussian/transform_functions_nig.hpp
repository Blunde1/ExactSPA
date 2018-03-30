// Copyright (C) 2016-2018 Berent Lunde

#ifndef __transform_functions_nig_hpp_included__
#define __transform_functions_nig_hpp_included__

#include "../complex.hpp"

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


#endif //__transform_functions_nig_hpp_included__
