// Copyright (C) 2018 Berent Lunde

#ifndef __EXACTSPA_HPP_INCLUDED__
#define __EXACTSPA_HPP_INCLUDED__

// External
#include <RcppEigen.h>
//#include "../external/Eigen/Dense"
#include "../external/adept.h"
#include "../external/adept_source.h"
#include <math.h>
//#include <string>

//using namespace Eigen;
using namespace std;
using namespace Rcpp;
using adept::adouble;
typedef adept::adouble adtype;

// Internal
#include "complex.hpp"

// Merton Jump Diffusion
#include "Merton_jump_diffusion/optimization_mjd.hpp"
#include "Merton_jump_diffusion/transform_functions_mjd.hpp"
#include "Merton_jump_diffusion/Fourier_inversion_mjd.hpp"
#include "Merton_jump_diffusion/saddlepoint_approximation_mjd.hpp"

// Noncentral Chi Square
#include "noncentral_chi_squared/optimization_nchisq.hpp"
#include "noncentral_chi_squared/transform_functions_nchisq.hpp"
#include "noncentral_chi_squared/Fourier_inversion_nchisq.hpp"
#include "noncentral_chi_squared/saddlepoint_approximation_nchisq.hpp"

// Normal inverse Gaussian
#include "normal_inverse_gaussian/optimization_nig.hpp"
#include "normal_inverse_gaussian/transform_functions_nig.hpp"
#include "normal_inverse_gaussian/Fourier_inversion_nig.hpp"
#include "normal_inverse_gaussian/saddlepoint_approximation_nig.hpp"

// Tweedie


#endif // __EXACTSPA_HPP_INCLUDED__
