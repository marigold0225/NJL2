#ifndef GAUSS_H
#define GAUSS_H

#include <vector>

namespace Gauss
{
    extern const double Pi;
    extern const std::vector<double> weights;
    extern const std::vector<double> nodes;

    extern const int int_phi_i, int_rho_i, int_E_i, int_phi_0, int_Z_i, int_E;

    double fermi_dirac(double E_i, double mu_star, double T_i);
    double gauss_integral(double down, double up, double (*integrand)(int case_id, double x, double mu_star, double T_i, double M_i), int case_id, double mu_star, double T_i, double M_i);
    double integrand_function(int case_id, double x, double mu_star, double T_i, double M_i);

}

#endif
