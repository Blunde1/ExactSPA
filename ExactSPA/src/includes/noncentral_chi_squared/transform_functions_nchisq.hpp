// Copyright (C) 2016-2018 Berent Lunde

#ifndef __transform_functions_nchisq_hpp_included__
#define __transform_functions_nchisq_hpp_included__

#include "../complex.hpp"

template<typename T>
struct cgf_nchisq{

    // Members
    T df, ncp;

    // Constructors
    cgf_nchisq(T df_, T ncp_) : df(df_), ncp(ncp_) {}
    cgf_nchisq(const cgf_nchisq& cgf_) : df(cgf_.df), ncp(cgf_.ncp) {}

    // Operators
    T operator()(T s){
        return ncp*s/(1.0-2.0*s) - df/2.0*log(1.0-2.0*s);
    }
    cType<T> operator()(cType<T> s){
        cType<T> res;
        res = ncp*s/(T(1.0)-T(2.0)*s) - cType<T>(df/T(2.0))*log(T(1.0)-T(2.0)*s);
        return res;
    }


    // Member functions
    T deriv(T s){
        return (ncp - 2.0*df*s+df)/pow(1.0-2.0*s,2);
    }
    T dderiv(T s){
        return 8.0*ncp*s/pow(1.0-2.0*s,3)+
            4.0*ncp/pow(1.0-2.0*s,2)+
            2.0*df/pow(1.0-2.0*s,2);
    }

};


/* standardized log characteristic function */
template<typename T>
struct lcf_stand_nchisq{
    // Members
    cgf_nchisq<T>& cgf;
    T x, sp;

    // Constructor
    lcf_stand_nchisq(cgf_nchisq<T> &cgf_, T x_, T sp_) : cgf(cgf_), x(x_), sp(sp_) {}

    // Operator
    cType<T> operator()(cType<T> s){
        cType<T> res, i(0,1);
        T h = cgf.dderiv(sp);
        res = cType<T>(-cgf(sp)) - i*s*x / cType<T>(sqrt(h)) +
            cgf(  sp + s*i / cType<T>(sqrt(h)) );

        return res;
    }
};


#endif //__transform_functions_nchisq_hpp_included__
