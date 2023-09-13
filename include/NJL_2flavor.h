#ifndef NJL_2flavor
#define NJL_2flavor
#include "Gauss.h"

namespace NJL_2flavors
{

    class Constants
    {
    public:
        // Constructor
        Constants();
        double compute_effect_mass(double m, double rho_1, double rho_2, double phi_1, double phi_2);
        double compute_rho(double cut_k, double mu_star, double T_i, double M_i);
        double compute_phi(double cut_k, double mu_star, double T_i, double M_i);
        double compute_che_star(double mu_i, double rho_1, double rho_2, double phi_1, double phi_2);
        void compute_self_coulping(double rho_input, double T_input, int case_id);
        void compute_omega(int x, double T_input);
        void compute_omega_zero(int x, double T_input);
        void compute_NJL2(int case_id);

    private:
        // Constants
        static constexpr double M_u = 5.5;
        static constexpr double M_d = 5.5;
        static constexpr double Lambda = 651;
        static constexpr double G_s = 2.135 / (Lambda * Lambda);
        static constexpr double Nc = 3;
        static constexpr double planck = 197.33;
        int n_steps;
        int t_steps;
        double G_sv;

        // State variables
        double rho_u;
        double rho_d;
        double rho;
        double effect_M_u;
        double effect_M_d;
        double phi_u;
        double phi_d;
        double phi_total;
        double che_u;
        double che_d;
        double che_u_star;
        double che_d_star;
        double che_total;
        double T;
        std::vector<double> pressure;
        std::vector<double> energy;
        std::vector<double> omega;
        std::vector<double> omega_zero;
    };
} // namespace NJL

#endif