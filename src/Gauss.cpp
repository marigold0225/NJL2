#include "Gauss.h"
#include <cmath>

namespace Gauss
{
    const double Pi = 3.1415926535;
    const int int_phi_i = 1, int_rho_i = 2, int_E_i = 3, int_phi_0 = 4, int_Z_i = 5, int_E = 6;
    const std::vector<double> weights = {0.3626837833783620, 0.3626837833783620,
                                         0.3137066458778873, 0.3137066458778873,
                                         0.2223810344533745, 0.2223810344533745,
                                         0.1012285362903763, 0.1012285362903763};
    const std::vector<double> nodes = {-0.1834346424956498, 0.1834346424956498,
                                       -0.5255324099163290, 0.5255324099163290,
                                       -0.7966664774136267, 0.7966664774136267,
                                       -0.9602898564975363, 0.9602898564975363};

    double fermi_dirac(double E_i, double mu_star, double T_i)
    {
        return 1.0 / (1.0 + exp((E_i - mu_star) / T_i));
    }

    double gauss_integral(double down, double up, double (*integrand)(int case_id, double x, double mu_star, double T_i, double M_i), int case_id, double mu_star, double T_i, double M_i)
    {
        double result = 0.0;
        double t, f;

        for (int i = 0; i < 8; i++)
        {
            t = (up - down) * nodes[i] / 2 + (up + down) / 2;
            f = integrand(case_id, t, mu_star, T_i, M_i);
            result += weights[i] * f;
        }

        result *= (up - down) / 2;
        return result;
    }

    double integrand_function(int case_id, double t, double mu_star, double T_i, double M_i)
    {
        double E = sqrt(t * t + M_i * M_i);
        double ni_plus, ni_minus, Z_plus, Z_minus;
        switch (case_id)
        {
        case int_phi_i:
            ni_plus = fermi_dirac(E, mu_star, T_i);
            ni_minus = fermi_dirac(E, -mu_star, T_i);
            return 4 * Pi * t * t * M_i * (ni_plus + ni_minus - 1) / E / (2 * Pi) / (2 * Pi) / (2 * Pi);
        case int_rho_i:
            ni_plus = fermi_dirac(E, mu_star, T_i);
            ni_minus = fermi_dirac(E, -mu_star, T_i);
            return 4 * Pi * t * t * (ni_plus - ni_minus) / (2 * Pi) / (2 * Pi) / (2 * Pi);
            break;
        case int_E_i:
            return 4 * Pi * t * t * E / (2 * Pi) / (2 * Pi) / (2 * Pi);
            break;
        case int_phi_0:
            return 4 * Pi * t * t * M_i / E / (2 * Pi) / (2 * Pi) / (2 * Pi);
            break;
        case int_Z_i:
            ni_plus = fermi_dirac(E, mu_star, -T_i);
            ni_minus = fermi_dirac(E, -mu_star, -T_i);
            Z_plus = T_i * log(1 / ni_plus);
            Z_minus = T_i * log(1 / ni_minus);
            return 4 * Pi * t * t * (Z_plus + Z_minus) / (2 * Pi) / (2 * Pi) / (2 * Pi);
            break;
        case int_E:
            ni_plus = fermi_dirac(E, mu_star, T_i);
            ni_minus = fermi_dirac(E, -mu_star, T_i);
            return 4 * Pi * t * t * E * (-ni_plus - ni_minus + 1) / (2 * Pi) / (2 * Pi) / (2 * Pi);
            break;
        default:
            return 0.0;
        }
    }
}
