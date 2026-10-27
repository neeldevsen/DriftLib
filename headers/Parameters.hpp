#pragma once
#include <type_traits>

template <typename S_t, typename K_t, typename r_t, typename sig_t, typename T_t, typename q_t = r_t>
class GBMParams
{
    public:
        using U = std::common_type_t<long double, S_t, K_t, r_t, sig_t, T_t, q_t>;
        U S {}; // spot
        U K {}; // strike
        U r {}; // interest rate
        U sig {}; // volatility
        U T {};  // time to maturity
        U q {}; // contiuous divided yield
};

template <typename ka_t, typename th_t, typename xi_t, typename rho_t>
class HestonParams
{
    public:
        using U = std::common_type_t<long double, ka_t, th_t, xi_t, rho_t>;
        U ka_v {}; // mean reversion speed
        U th_v {}; // long run average volatility
        U xi_v {}; // volatility of volatility
        U rho_v {}; // correlation factor
};

template <typename la_t, typename muj_t, typename sij_t>
class JumpDiffusionsParams
{
    public:
        using U = std::common_type_t<long double, la_t, muj_t, sij_t>;
        U la {}; // jump intensity 
        U mu_j {}; // mean of log jump size
        U si_j {}; // std of log jump size
};

template <typename r0_t, typename a_t, typename sir_t, typename thr_t>
class HullWhiteParams
{
    public:
        using U = std::common_type_t<long double, r0_t, a_t, sir_t, thr_t>;
        U r0 {};
        U a {};
        U sir {};
        U th {};
};
