// Copyright (C) 2016-2018 Berent Lunde

#ifndef __OPTIMIZATION_TWEEDIE_HPP_INCLUDED__
#define __OPTIMIZATION_TWEEDIE_HPP_INCLUDED__

#include "transform_functions_tweedie.hpp"

// Line search Newton method
template<typename T>
T sp_newton_opt(T s_start, T x, cgf_tweedie<T>& cgf, bool check_values=true, int maxniter=1000){
    int i;
    double alpha = 1.0;
    T EPS = 3.0e-12, f_eps = 1.0e50;
    T g, h, s_old=s_start, s_new=s_start;
    for(i=0; i<maxniter; i++){
        g = cgf.deriv(s_old) - x;
        h = cgf.dderiv(s_old);
        s_new = s_old - alpha * g/h;
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


#endif //__OPTIMIZATION_TWEEDIE_HPP_INCLUDED__
