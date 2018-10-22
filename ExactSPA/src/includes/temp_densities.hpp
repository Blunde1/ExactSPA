// Copyright (C) 2018 Berent Lunde

#ifndef __TEMP_DENSITIES_HPP_INCLUDED__
#define __TEMP_DENSITIES_HPP_INCLUDED__


template<typename T>
T pois_dens(int x, T lambda){
    T lfactorialX = 0.0;
    for(int i=x; i > 0; i--){
        lfactorialX = lfactorialX + log(i);
    }
    T lfx = x*log(lambda) - lambda - lfactorialX;
    return exp(lfx);
}

template<typename T>
T norm_dens(T x, T mu, T sigma){
    T lfx = - 0.5 * (x-mu)*(x-mu)/(sigma*sigma) - 0.5*log(2.0*M_PI) - log(sigma);
    return exp(lfx);
}

#endif // __TEMP_DENSITIES_HPP_INCLUDED__
