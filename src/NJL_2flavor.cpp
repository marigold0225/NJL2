#include "Gauss.h"
#include "NJL_2flavor.h"
#include <cmath>
#include <iostream>
#include <string>
#include <fstream>

using namespace std;

namespace NJL_2flavors
{
    Constants::Constants()
    {
        n_steps = 1500;
        t_steps = 120;
        pressure.resize(n_steps);
        energy.resize(n_steps);
        omega.resize(n_steps);
        omega_zero.resize(n_steps);
    }
    double Constants::compute_effect_mass(double m, double rho1, double rho2, double phi1, double phi2)
    {
        double M_i = m - 2 * G_s * (phi1 + phi2) - 2 * G_sv * (phi1 + phi2) * pow((rho1 + rho2), 2) * pow(planck, 6);
        return M_i;
    }
    double Constants::compute_che_star(double mu_i, double rho1, double rho2, double phi1, double phi2)
    {
        double mu_star = mu_i + 2 * G_sv * (rho1 + rho2) * pow((phi1 + phi2), 2) * pow(planck, 3);
        return mu_star;
    }
    double Constants::compute_phi(double cut_k, double mu_star, double T_i, double M_i)
    {
        double phi_i = 2 * Nc * Gauss::gauss_integral(0.0, cut_k, Gauss::integrand_function, 1, mu_star, T_i, M_i);
        return phi_i;
    }
    double Constants::compute_rho(double cut_k, double mu_star, double T_i, double M_i)
    {
        double rho_i = 2 * Nc * Gauss::gauss_integral(0.0, cut_k, Gauss::integrand_function, 2, mu_star, T_i, M_i) / planck / planck / planck;
        return rho_i;
    }
    void Constants::compute_omega(int index, double T_input)
    {
        if (index >= 0 && index < n_steps)
        {
            omega[index] = G_s * pow((phi_u + phi_d), 2) + 3 * G_sv * pow((phi_u + phi_d), 2) * pow((rho_u + rho_d), 2) * pow(planck, 6) - 2 * Nc * (Gauss::gauss_integral(0.0, Lambda, Gauss::integrand_function, 3, che_u_star, T_input, effect_M_u) + Gauss::gauss_integral(0.0, Lambda, Gauss::integrand_function, 3, che_d_star, T_input, effect_M_d)) - 2 * Nc * (Gauss::gauss_integral(0.0, Lambda, Gauss::integrand_function, 5, che_u_star, T_input, effect_M_u) + Gauss::gauss_integral(0.0, Lambda, Gauss::integrand_function, 5, che_d_star, T_input, effect_M_d));
        }

        else
        {
            cout << "out of range" << endl;
        }
    }
    void Constants::compute_omega_zero(int index, double T_input)
    {
        if (index >= 0 && index < n_steps)
        {
            double phi_u0 = -2 * Nc * Gauss::gauss_integral(0.0, Lambda, Gauss::integrand_function, 4, che_u_star, T_input, 325.1);
            double phi_d0 = -2 * Nc * Gauss::gauss_integral(0.0, Lambda, Gauss::integrand_function, 4, che_d_star, T_input, 325.1);
            omega_zero[index] = -2 * Nc * (Gauss::gauss_integral(0.0, Lambda, Gauss::integrand_function, 3, che_u_star, T_input, 325.1) + Gauss::gauss_integral(0.0, Lambda, Gauss::integrand_function, 3, che_d_star, T_input, 325.1)) + G_s * pow((phi_u0 + phi_d0), 2);
        }
        else
        {
            cout << "out of range" << endl;
        }
    }
    void Constants::compute_self_coulping(double rho_input, double T_input, int case_id)
    {
        double rho_u_new, rho_d_new;
        double phi_u_new, phi_d_new;
        double up_phi, down_phi;
        double up_che, down_che;
        double diff_phi, diff_rho, phi_total, che_total;
        const double tol_rho = 1.0e-5, tol_phi = 2.0;
        const int max_outer_iterations = 1000, max_inner_iterations = 1000;
        int outer_iteration, inner_iteration;

        rho_u = rho_input / 2.;
        rho_d = rho_input / 2.;
        down_phi = -2.0 * pow(252, 3);
        up_phi = 0;

        bool converged = false;

        for (outer_iteration = 1; outer_iteration <= max_outer_iterations; ++outer_iteration)
        {
            phi_total = 0.5 * (up_phi + down_phi);
            up_che = 2000;
            down_che = 0;
            phi_u = phi_total / 2.0;
            phi_d = phi_total / 2.0;
            effect_M_u = compute_effect_mass(M_u, rho_u, rho_d, phi_u, phi_d);
            effect_M_d = compute_effect_mass(M_d, rho_u, rho_d, phi_u, phi_d);

            for (inner_iteration = 1; inner_iteration <= max_inner_iterations; ++inner_iteration)
            {
                che_total = 0.5 * (up_che + down_che);
                che_u = che_total / 2.0;
                che_d = che_total / 2.0;

                che_u_star = compute_che_star(che_u, rho_u, rho_d, phi_u, phi_d);
                che_d_star = compute_che_star(che_d, rho_u, rho_d, phi_u, phi_d);

                phi_u_new = compute_phi(Lambda, che_u_star, T_input, effect_M_u);
                phi_d_new = compute_phi(Lambda, che_d_star, T_input, effect_M_d);

                diff_phi = phi_u_new - phi_total / 2.0;
                if (diff_phi < 0.0)
                {
                    down_che = che_total;
                }
                else
                {
                    up_che = che_total;
                }

                if (std::abs(diff_phi) < tol_phi)
                    break;
            }

            rho_u_new = compute_rho(Lambda, che_u_star, T_input, effect_M_u);
            rho_d_new = compute_rho(Lambda, che_d_star, T_input, effect_M_d);

            diff_rho = rho_u_new - rho_input / 2.0;
            if (diff_rho < 0.0)
            {
                down_phi = phi_total;
            }
            else
            {
                up_phi = phi_total;
            }

            if (std::abs(diff_rho) < tol_rho)
            {
                converged = true;
                break;
            }
        }
        switch (case_id)
        {
        case 1:
            if (converged)
            {
                std::cout << "T_input: " << T_input << "\n"
                          << "Converged in " << outer_iteration << " outer iterations. "
                          << "Converged in " << inner_iteration << " inner iterations.\n"
                          << "phi_u: " << phi_u / pow(planck, 3) << " phi_d: " << phi_d / pow(planck, 3) << "\n"
                          << "che_u: " << che_u << " che_d: " << che_d << " che_total: " << che_total << "\n"
                          << "che_u_star: " << che_u_star << " che_d_star: " << che_d_star << "\n"
                          << "effect_M_u: " << effect_M_u << " effect_M_d: " << effect_M_d << "\n"
                          << "rho_u: " << rho_u << " rho_d: " << rho_d << "\n"
                          << "diff_phi: " << diff_phi << " diff_rho: " << diff_rho << "\n"
                          << "-----------------------------\n";
            }
            else
            {
                std::cout << "Didn't converge in " << outer_iteration << " iterations.\n";
            }
            break;
        case 2:
            // Handle this case
            break;
        }
    }
    void Constants::compute_NJL2(int case_id)
    {

        double rho_step, T_steps;
        std::vector<double> G_sv_values(3);
        std::vector<double> rho_values(n_steps);
        std::vector<double> che_values(n_steps);
        std::vector<double> dp_drho(n_steps);
        switch (case_id)
        {
        case 1:
            G_sv_values = {0.0, 100.0 / std::pow(Lambda, 8), -100.0 / std::pow(Lambda, 8)};
            for (int j = 0; j < 3; j++)
            {
                G_sv = G_sv_values[j];
                std::string filename = "P_rho_G_sv_" + std::to_string(j + 1) + ".txt";
                std::ofstream outFile(filename);

                T = 36;
                rho_step = 1.5 / n_steps;
                rho = 0.0;

                for (int i = 0; i < n_steps; i++)
                {
                    rho += rho_step;
                    rho_values[i] = rho;
                    compute_self_coulping(rho, T, case_id);
                    compute_omega(i, T);
                    compute_omega_zero(i, T);
                    pressure[i] = (-omega[i] + omega_zero[i]) / std::pow(planck, 3);
                    che_values[i] = che_u;
                    outFile << rho_values[i] << " " << pressure[i] << " " << che_values[i] << std::endl;
                }

                outFile.close();
                std::cout << "Calculation for G_sv = " << G_sv << " completed. Results saved to " << filename << std::endl;
            }
            break;

        case 2:
            G_sv_values = {0.0, 100.0 / std::pow(Lambda, 8), -200.0 / std::pow(Lambda, 8)};
            for (int j = 0; j < 3; j++)
            {
                G_sv = G_sv_values[j];
                std::string filename = "spinodal_G_sv_" + std::to_string(j + 1) + ".txt";
                std::ofstream outFile(filename);
                T_steps = 120.0 / t_steps;
                rho_step = 1.5 / n_steps;
                T = 2.0;
                for (int k = 0; k < t_steps; k++)
                {
                    T += T_steps;
                    rho = 0.0;
                    for (int i = 0; i < n_steps; i++)
                    {
                        rho += rho_step;
                        rho_values[i] = rho;
                        compute_self_coulping(rho, T, case_id);
                        compute_omega(i, T);
                        compute_omega_zero(i, T);
                        pressure[i] = (-omega[i] + omega_zero[i]) / std::pow(planck, 3);
                        che_values[i] = che_u;
                    }
                    for (int i = 0; i < n_steps; i++)
                    {
                        if (i == 0)
                        {
                            dp_drho[i] = (pressure[1] - pressure[0]) / (rho_values[1] - rho_values[0]);
                        }
                        else if (i == n_steps - 1)
                        {
                            dp_drho[i] = (pressure[n_steps - 1] - pressure[n_steps - 2]) / (rho_values[n_steps - 1] - rho_values[n_steps - 2]);
                        }
                        else
                        {
                            dp_drho[i] = (pressure[i + 1] - pressure[i - 1]) / (rho_values[i + 1] - rho_values[i - 1]);
                        }
                    }
                    bool positive = true;
                    for (int i = 0; i < n_steps; i++)
                    {
                        if (dp_drho[i] <= 0.0)
                        {
                            positive = false;
                            break;
                        }
                    }
                    if (positive)
                    {
                        std::cout << "the derivatives > 0 for T=:" << T << std::endl;
                        break;
                    }
                    else
                    {
                        for (int i = 0; i < n_steps - 1; i++)
                        {
                            if (dp_drho[i] * dp_drho[i + 1] <= 0.0)
                            {
                                double derivations_point = (rho_values[i] + rho_values[i + 1]) / 2;
                                std::cout << "rho=:" << derivations_point << "T=:" << T << std::endl;
                                outFile << derivations_point << " " << T << std::endl;
                            }
                        }
                    }
                }
                outFile.close();
                std::cout << "Calculation for G_sv = " << G_sv << " completed. Results saved to " << filename << std::endl;
            }
            break;

        default:
            std::cout << "Invalid case_id" << std::endl;
            break;
        }
        std::cout << "All calculations completed." << std::endl;
    }
}
