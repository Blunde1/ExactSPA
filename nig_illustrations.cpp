// Berent Lunde
// 29.05.2018

// Load exactspa lib

#include <math.h>

#include "ExactSPA/src/includes/complex.hpp"
// Normal inverse Gaussian
#include "ExactSPA/src/includes/normal_inverse_gaussian/optimization_nig.hpp"
#include "ExactSPA/src/includes/normal_inverse_gaussian/transform_functions_nig.hpp"
#include "ExactSPA/src/includes/normal_inverse_gaussian/Fourier_inversion_nig.hpp"
#include "ExactSPA/src/includes/normal_inverse_gaussian/saddlepoint_approximation_nig.hpp"


// [[Rcpp::export]]
double reStandCFNIG(
        double s, // stand lcf evaluation
        double x, // nig realisation
        double lchi, double lpsi, double mu, double gamma // param
                       )
{
    // Give interface to Re stand cf
    
    // Find saddepoint
    double chi = exp(lchi), psi = exp(lpsi);
    cgf_nig<double> cgf(chi, psi, mu, gamma);
    double sp = sp_newton_opt(0.0, x, cgf);
    lcf_stand_nig<double> lcf_stand( cgf, x, sp );
    
    // Cast to complex
    cType<double> s_complex(s);
    
    return real(exp(lcf_stand(s_complex)));
}



// [[Rcpp::export]]
double reCFNIG(
        double s, // stand lcf evaluation
        double x, // nig realisation
        double lchi, double lpsi, double mu, double gamma // param
)
{
    // Give interface to Re cf
    
    double chi = exp(lchi), psi = exp(lpsi);
    cgf_nig<double> cgf(chi, psi, mu, gamma);
    // Cast to complex
    cType<double> s_complex(s);
    cType<double> i(0,1);
    return real(exp(cgf(-i*s_complex)));
}