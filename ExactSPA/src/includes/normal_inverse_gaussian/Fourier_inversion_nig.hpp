// Copyright (C) 2016-2018 Berent Lunde

#ifndef __FOURIER_INVERSION_NIG_HPP_INCLUDED__
#define __FOURIER_INVERSION_NIG_HPP_INCLUDED__

#include "../complex.hpp"


/** Simpson's rule for IFT
 */
template<typename T>
T ift_simpson0(lcf_stand_nig<T> lcf_stand, double length, int n){
    T simpson_integral = 0.0, tmp_val = 0.0;
    T a = 0.0, b = length;
    T h = (b-a) / n;

    // a
    tmp_val = real(exp(lcf_stand(cType<T>(a))));
    if(isfinite(tmp_val)){
        simpson_integral += tmp_val;
    }
    // b
    tmp_val = real(exp(lcf_stand(cType<T>(b))));
    if(isfinite(tmp_val)){
        simpson_integral += tmp_val;
    }

    //simpson_integral = real(exp(lcf_stand(cType<T>(a)))) + real(exp(lcf_stand(cType<T>(b))));

    for(int j=1; j<=(n/2); j++){
        tmp_val = 4.0 * real(exp(lcf_stand(cType<T>(h*(2*j-1)))));
        if(isfinite(tmp_val)){
            simpson_integral += tmp_val;
        }
    }
    for(int j=1; j<=(n/2-1); j++){
        tmp_val = 2.0 * real(exp(lcf_stand(cType<T>(h*(2*j)))));
        if(isfinite(tmp_val)){
            simpson_integral += tmp_val;
        }
    }

    simpson_integral = simpson_integral*h/(3.0*M_PI);

    return simpson_integral;
}

/* Simpson's rule for IFT - non-standardised cf */
template<typename T>
T ift_simpson(cgf_nig<T> cgf, T x, double length, int n){
    T simpson_integral = 0.0;
    T a = 0.0, b = length;
    T h = (b-a) / n;
    cType<T> i(0,1);

    simpson_integral = real(exp(cgf(cType<T>(a*i)) - i*a*x )) + real(exp(cgf(cType<T>(b*i)) - i*b*x));

    for(int j=1; j<=(n/2); j++){
        simpson_integral += 4.0 * real(exp(cgf(cType<T>(i*cType<T>(h*(2*j-1)))) - i*cType<T>(h*(2.0*j-1)*x)));
    }
    for(int j=1; j<=(n/2-1); j++){
        simpson_integral += 2.0 * real(exp(cgf(cType<T>(i*cType<T>(h*(2*j)))) - i*cType<T>(h*(2.0*j)*x) ));
    }
    simpson_integral = simpson_integral*h/(3.0*M_PI);

    return simpson_integral;
}


#endif //__FOURIER_INVERSION_NIG_HPP_INCLUDED__
