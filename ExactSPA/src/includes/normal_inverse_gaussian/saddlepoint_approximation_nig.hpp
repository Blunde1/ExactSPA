// Copyright (C) 2016-2018 Berent Lunde

#ifndef __SADDLEPOINT_APPROXIMATION_NIG_HPP_INCLUDED__
#define __SADDLEPOINT_APPROXIMATION_NIG_HPP_INCLUDED__

#include "transform_functions_nig.hpp"

/* Log exact spa */
template<typename T>
T log_exact_spa(T x, cgf_nig<T> cgf, T sp, lcf_stand_nig<T> lcf_stand, double length, int n){

    T lfx, fx_standardized;
    fx_standardized = ift_simpson0(lcf_stand, length, n);
    lfx = cgf(sp) - sp*x + log(fx_standardized) -
        0.5*log(cgf.dderiv(sp));

    return lfx;
}


#endif //__SADDLEPOINT_APPROXIMATION_NIG_HPP_INCLUDED__
