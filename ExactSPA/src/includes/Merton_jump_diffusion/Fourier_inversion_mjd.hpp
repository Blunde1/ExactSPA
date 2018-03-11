// Copyright (C) 2016-2018 Berent Lunde

#ifndef __FOURIER_INVERSION_HPP_INCLUDED__
#define __FOURIER_INVERSION_HPP_INCLUDED__

#include "../complex.hpp"


/** Simpson's rule for IFT
 */
template<typename T>
T ift_simpson0(lcf_stand_mjd<T> lcf_stand, double length, int n){
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


#endif //__FOURIER_INVERSION_HPP_INCLUDED__
