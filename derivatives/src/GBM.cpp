#include <iostream>
#include <thread>
#include <random>
#include <cmath>
#include <typeinfo>
#include <type_traits>
#include <optional>
#include <iomanip>
#include <vector>
#include "../headers/Parameters.hpp"
#include "../headers/Algorithms.hpp"
#include "../headers/Functions.hpp"
#include "../headers/Structures.hpp"

template <typename U, typename V>
class Options
{
    public:
        U call {};
        V put {};
};

template <typename S_t, typename K_t, typename r_t, typename sig_t, typename T_t, typename q_t = r_t>
class GBM: private GBMParams<S_t, K_t, r_t, sig_t, T_t, q_t>
{
    private:
        template <typename U>
        class FDMCoefficients
        {
            private:
                std::vector<U> alphas {};
                std::vector<U> betas {};
                std::vector<U> gammas {};
            
            public:
                explicit FDMCoefficients(const int M)
                : alphas(M-1), betas(M-1), gammas(M-1) {}

                friend class GBM;

                void generateCoefficients(const GBM& params, const U dS, const int M)
                {
                    for (int i{0}; i < M-1; ++i)
                    {
                        U S_i {(i+1)*dS};
                        (this->alphas)[i] = static_cast<U>(0.5) * (std::pow(params.sig*S_i, static_cast<U>(2))) / (std::pow(dS, static_cast<U>(2))) - ((params.r - params.q) * S_i) / (static_cast<U>(2) * dS);
                        (this->betas)[i] = -(std::pow(params.sig * S_i, static_cast<U>(2))) / (std::pow(dS, static_cast<U>(2))) - params.r;
                        (this->gammas)[i] = static_cast<U>(0.5) * (std::pow(params.sig*S_i, static_cast<U>(2))) / (std::pow(dS, static_cast<U>(2))) + ((params.r - params.q) * S_i) / (static_cast<U>(2) *dS);
                    }
                }

                Structures::TridiagonalMatrix<U> generateDiagonals(const U dT, const int M)
                {
                    Structures::TridiagonalMatrix<U> matrix {M};

                    for (int i {0}; i < M-2; ++i)
                    {
                        matrix.middle[i] = static_cast<U>(1) - static_cast<U>(0.5) * dT * (this->betas)[i];
                        matrix.lower[i] = - static_cast<U>(0.5) * dT * (this->alphas)[i+1];
                        matrix.upper[i] = - static_cast<U>(0.5) * dT * (this->gammas)[i];
                    }
                    matrix.middle[M-2] = static_cast<U>(1) - static_cast<U>(0.5) * dT * (this->betas)[M-2];
                    return matrix;
                }

                std::vector<U> generateRHSCallCN(const GBM& params, std::vector<U> V, const int n, const int M, const U dT, const U S_max) 
                {
                    V[M] = S_max * std::exp(-params.q * (static_cast<U>(n)) * dT) - params.K * std::exp(-params.r*static_cast<U>(n)*dT); 
                    V[0] = 0;
                    std::vector<U> RHS (M-1);
                    for (int i{}; i < M-1; ++i)
                    {
                        RHS[i] = static_cast<U>(0.5) * dT * (this->alphas[i]) * V[i] + (static_cast<U>(1) + static_cast<U>(0.5)* dT * this->betas[i])*V[i+1] + static_cast<U>(0.5)*  dT * (this->gammas)[i] * V[i+2];
                    }
                    RHS[M - 2] += (static_cast<U>(0.5) * dT * (this->gammas)[M - 2]) * (S_max * std::exp(-params.q * (static_cast<U>(n + 1)) * dT) - params.K * std::exp(-params.r * (static_cast<U>(n + 1)) * dT));
                    return RHS; 
                }

                std::vector<U> generateRHSPutCN(const GBM& params, std::vector<U> V, const int n, const int M, const U dT, const U S_max) 
                {
                    V[0] = params.K * std::exp(-params.r*static_cast<U>(n)*dT); 
                    V[M] = 0;
                    std::vector<U> RHS (M-1);
                    for (int i{}; i < M-1; ++i)
                    {
                        RHS[i] = static_cast<U>(0.5) * dT * (this->alphas[i]) * V[i] + (static_cast<U>(1) + static_cast<U>(0.5)* dT * this->betas[i])*V[i+1] + static_cast<U>(0.5)*  dT * (this->gammas)[i] * V[i+2];
                    }
                    RHS[0] += (static_cast<U>(0.5) * dT * (this->alphas)[0]) * (params.K * std::exp(-params.r * (static_cast<U>(n + 1)) * dT));
                    return RHS; 
                }

                std::vector<U> generateRHSCallBE(const GBM& params, std::vector<U> V, const int n, const int M, const U dT, const U S_max) 
                {
                    V[M] = S_max * std::exp(-params.q * (static_cast<U>(n)) * dT) - params.K * std::exp(-params.r*static_cast<U>(n)*dT); 
                    std::vector<U> RHS (M-1);
                    for (int i{}; i < M-1; ++i)
                    {
                        RHS[i] = V[i+1];
                    }
                    RHS[M - 2] += (dT * (this->gammas)[M - 2]) * (S_max * std::exp(-params.q * (static_cast<U>(n + 1)) * dT) - params.K * std::exp(-params.r * (static_cast<U>(n + 1)) * dT));
                    return RHS; 
                }

                std::vector<U> generateRHSPutBE(const GBM& params, std::vector<U> V, const int n, const int M, const U dT, const U S_max) 
                {
                    V[0] = params.K * std::exp(-params.r*static_cast<U>(n)*dT); 
                    std::vector<U> RHS (M-1);
                    for (int i{}; i < M-1; ++i)
                    {
                        RHS[i] = V[i+1];                    
                    }
                    RHS[0] += (dT * (this->alphas)[0]) * (params.K * std::exp(-params.r * (static_cast<U>(n + 1)) * dT));
                    return RHS; 
                }

                // American now

                std::vector<U> generateRHSCallCNAmerican(const GBM& params, std::vector<U> V, const int n, const int M, const U dT, const U S_max) 
                {
                    V[M] = std::max(S_max - params.K , S_max * std::exp(-params.q * (static_cast<U>(n)) * dT) - params.K * std::exp(-params.r*static_cast<U>(n)*dT)); 
                    V[0] = 0;
                    std::vector<U> RHS (M-1);
                    for (int i{}; i < M-1; ++i)
                    {
                        RHS[i] = static_cast<U>(0.5) * dT * (this->alphas[i]) * V[i] + (static_cast<U>(1) + static_cast<U>(0.5)* dT * this->betas[i])*V[i+1] + static_cast<U>(0.5)*  dT * (this->gammas)[i] * V[i+2];
                    }
                    RHS[M - 2] += (static_cast<U>(0.5) * dT * (this->gammas)[M - 2]) * std::max(S_max - params.K , S_max * std::exp(-params.q * (static_cast<U>(n + 1)) * dT) - params.K * std::exp(-params.r*static_cast<U>(n + 1)*dT)); 
                    return RHS; 
                }

                std::vector<U> generateRHSPutCNAmerican(const GBM& params, std::vector<U> V, const int n, const int M, const U dT, const U S_max) 
                {
                    V[0] = params.K;
                    V[M] = 0;
                    std::vector<U> RHS (M-1);
                    for (int i{}; i < M-1; ++i)
                    {
                        RHS[i] = static_cast<U>(0.5) * dT * (this->alphas[i]) * V[i] + (static_cast<U>(1) + static_cast<U>(0.5)* dT * this->betas[i])*V[i+1] + static_cast<U>(0.5)*  dT * (this->gammas)[i] * V[i+2];
                    }
                    RHS[0] += (static_cast<U>(0.5) * dT * (this->alphas)[0]) * (params.K);
                    return RHS; 
                }

                std::vector<U> generateRHSCallBEAmerican(const GBM& params, std::vector<U> V, const int n, const int M, const U dT, const U S_max) 
                {
                    V[M] = std::max(S_max - params.K , S_max * std::exp(-params.q * (static_cast<U>(n)) * dT) - params.K * std::exp(-params.r*static_cast<U>(n)*dT)); 
                    std::vector<U> RHS (M-1);
                    for (int i{}; i < M-1; ++i)
                    {
                        RHS[i] = V[i+1];
                    }
                    RHS[M - 2] += (dT * (this->gammas)[M - 2]) * std::max(S_max - params.K , S_max * std::exp(-params.q * (static_cast<U>(n + 1)) * dT) - params.K * std::exp(-params.r*static_cast<U>(n + 1)*dT)); 
                    return RHS; 
                }

                std::vector<U> generateRHSPutBEAmerican(const GBM& params, std::vector<U> V, const int n, const int M, const U dT,  const U S_max) 
                {
                    V[0] = params.K;
                    V[M] = 0;
                    std::vector<U> RHS (M-1);
                    for (int i{}; i < M-1; ++i)
                    {
                        RHS[i] = V[i+1];                    
                    }
                    RHS[0] += (dT * (this->alphas)[0]) * (params.K);
                    return RHS; 
                }
        };


    public:
        GBM(S_t S, K_t K, r_t r, sig_t sig, T_t T, q_t q=0.0) 
        : GBMParams<S_t, K_t, r_t, sig_t, T_t, q_t>{S, K, r, sig, T, q} 
        {
            if (sig <= 0.0 || T <= 0.0 || S <= 0.0 || K <= 0.0)
            {
                throw std::invalid_argument("Invalid argument given\n");
            }
        }

        using U = std::common_type_t<long double, S_t, K_t, r_t, sig_t, T_t, q_t>;

        U putCallParity(const Options<U,U>& option) const 
        {
            return this->S * std::exp(-this->q * this->T)  + option.put - this->K * std::exp(-this->r * this->T) - option.call;
        }

        Options<U,U> closedForm() const noexcept
        {
            Options<U,U> option {};
            U sqrt_maturity {std::sqrt(this->T)};

            auto d1 {(std::log(this->S / this->K) + (this->r - this->q + (static_cast<U>(0.5) * this->sig * this->sig)) * (this->T))/ (this->sig * sqrt_maturity)};
            auto d2 {d1 - this->sig * sqrt_maturity};

            option.call = this->S  * std::exp(-this->q * this->T) * Functions::normalCDF(d1) - this->K * std::exp(-this->r * this->T) * Functions::normalCDF(d2);
            option.put = this->K * std::exp(-this->r * this->T) * Functions::normalCDF(-d2) - this->S * std::exp(-this->q * this->T) * Functions::normalCDF(-d1);
            return option;
        }

        Options<U,U> monteCarlo(const int N, bool av=false, bool cv=false, const std::optional<unsigned int> seed=std::nullopt) const
        {
            if (N <= 0)
            {
                throw std::invalid_argument("Invalid argument given\n");
            }
            else if (av && N % 2 == 1)
            {
                throw std::invalid_argument("Even number of simulations must be given for antithetic variates \n");
            }

            Options<U,U> option {};
            U exp_rt {std::exp(-this->r * this->T)};

            U sumY {};
            U sumY2 {};

            U sumXcall {};
            U sumXYcall {};

            U sumXput {};
            U sumXYput {};

            std::mt19937 mt{seed ? *seed : std::random_device{}()};
            
            std::normal_distribution W {static_cast<U>(0), std::sqrt(this->T)};

            for (int i {}; i < N; i += (av ? 2 : 1))
            {
                U Z {W(mt)};
                U S_final {this->S * std::exp((this->r -  this->q - static_cast<U>(0.5) * this->sig * this->sig) * this->T + (this->sig * Z))};
                
                
     
                sumXcall += exp_rt * std::max(S_final - this->K, static_cast<U>(0));
                sumXput += exp_rt * std::max(this->K - S_final, static_cast<U>(0));

                if (cv)
                {
                    U X_call {exp_rt * std::max(S_final - this->K, static_cast<U>(0))};
                    U X_put {exp_rt * std::max(this->K - S_final, static_cast<U>(0))};
                    U Y {exp_rt * S_final};

                    sumY += Y;
                    sumY2 += std::pow(Y, static_cast<U>(2));

                    sumXYcall += X_call * Y;
                    sumXYput += X_put * Y;
                }

                if (av)
                {
                    S_final = this->S * std::exp((this->r -  this->q - static_cast<U>(0.5) * this->sig * this->sig) * this->T - (this->sig * Z));
                    U X_call {exp_rt * std::max(S_final - this->K, static_cast<U>(0))};
                    sumXcall += X_call;
                    U X_put {exp_rt * std::max(this->K - S_final, static_cast<U>(0))};
                    sumXput += X_put;

                    if (cv)
                    {
                        U Y {exp_rt * S_final};
                        sumY += Y;
                        sumY2 += std::pow(Y, static_cast<U>(2));
                        sumXYcall += X_call * Y;
                        sumXYput += X_put * Y;
                    }
                }
            }

            sumXcall /= N;
            sumXput /= N;
            sumXYcall /= N;
            sumXYput /= N;
            sumY /= N;
            sumY2 /= N;
            
            U covXYcall {sumXYcall - sumXcall * sumY};
            U covXYput {sumXYput - sumXput * sumY};
            U varY {sumY2 - std::pow(sumY, static_cast<U>(2))};

            option.call = cv ? sumXcall - (covXYcall / varY) * (sumY - this->S*std::exp(-this->q*this->T)) : sumXcall;
            option.put = cv ? sumXput - (covXYput / varY) * (sumY - this->S*std::exp(-this->q*this->T)) : sumXput;
            return option;
        }

        Options<U,U> binomialTree(const int N, bool american=false) const
        {
            if (N <= 0)
            {
                throw std::invalid_argument("Invalid argument given\n");
            }

            Options<U,U> option {};

            U sqrt_deltat {std::sqrt(this->T / static_cast<U>(N))};

            U u {std::exp(this->sig * sqrt_deltat)};
            U d {std::exp(-this->sig * sqrt_deltat)};
            U p {(std::exp((this->r - this->q) * this->T / static_cast<U>(N)) - d) / (u - d)};

            if (p < 0 || p > 1)
            {
                throw std::invalid_argument("Probability is outside normal probability values.\n");
            }

            std::vector<U> V_call (N + 1);
            std::vector<U> V_put (N + 1);
            for (int j {}; j < N + 1; ++j)
                {
                    U S_ij {this->S * std::pow(static_cast<U>(u), static_cast<U>(j)) * std::pow(static_cast<U>(d), static_cast<U>(N-j))};
                    V_call[j] = std::max(S_ij - this->K, static_cast<U>(0));
                    V_put[j] = std::max(this->K - S_ij, static_cast<U>(0));
                }

            if (!american)
            {
            for (int i {}; i < N ; ++i)
            {
                for (int k {}; k < N-i; ++k)
                {
                    V_call[k] = p * V_call[k+1] + (1-p) * V_call[k];
                    V_put[k] = p * V_put[k+1] + (1-p) * V_put[k];
                }
                V_call.pop_back();
                V_put.pop_back();
            }   
      
                option.call = std::exp(-this->r * this->T) * V_call[0];
                option.put = std::exp(-this->r * this->T) * V_put[0];
                return option;
            }

            else
            {
                U exp_rdt {std::exp(-this->r * this->T / static_cast<U>(N))};
                U S_ij {0};
                U continuation_call {0};
                U intrinsic_call {0};

                U continuation_put {0};
                U intrinsic_put {0};

                for (int i {}; i < N ; ++i)
                {
                    for (int k {}; k < N-i; ++k)
                    {
                        S_ij = this->S * std::pow(static_cast<U>(u) , static_cast<U>(k)) * std::pow(static_cast<U>(d), static_cast<U>(N-i-k-1));

                        continuation_call = exp_rdt * (p * V_call[k+1] + (1-p) * V_call[k]);
                        intrinsic_call = std::max(S_ij - this->K , static_cast<U>(0));
                        V_call[k] = std::max(continuation_call, intrinsic_call);
                        
                        continuation_put = exp_rdt * (p * V_put[k+1] + (1-p) * V_put[k]);
                        intrinsic_put = std::max(this->K - S_ij, static_cast<U>(0));
                        V_put[k] = std::max(continuation_put, intrinsic_put);
                    }

                V_call.pop_back();
                V_put.pop_back();
            }

                option.call = V_call[0];
                option.put = V_put[0];
                return option;
            }    
        }

        Options<U,U> trinomialTree (const int N, bool american=false) const
        {
            if (N <= 0)
            {
                throw std::invalid_argument("Invalid argument\n");
            }
            
            Options<U,U> option {};
            U sqrt_2deltat {std::sqrt(static_cast<U>(2) * this->T / static_cast<U>(N))};

            auto u {std::exp(this->sig * sqrt_2deltat)};
            auto d {std::exp(-this->sig * sqrt_2deltat)};
            auto p_u {std::pow((std::exp((this->r - this->q) * this->T / (2 * N)) - std::exp(-this->sig * sqrt_2deltat / static_cast<U>(2) )) / ((std::exp(this->sig * sqrt_2deltat / static_cast<U>(2))) - std::exp(-this->sig * sqrt_2deltat / static_cast<U>(2))), static_cast<U>(2))};
            auto p_d {std::pow((std::exp(this->sig * sqrt_2deltat / static_cast<U>(2)) - std::exp((this->r - this->q) * this->T / (2 * N))) / ((std::exp(this->sig * sqrt_2deltat / static_cast<U>(2))) - std::exp(-this->sig * sqrt_2deltat / static_cast<U>(2))), static_cast<U>(2))};
            auto p_m {static_cast<U>(1) - p_u - p_d};
            std::vector<U> V_call (2 * N + 1);
            std::vector<U> V_put (2 * N + 1);
            for (int j {}; j < 2 * N + 1; ++j)
                {
                    U S_ij {this->S * std::pow(static_cast<U>(u) , static_cast<U>(j- N))};
                    V_call[j] = std::max(S_ij - this->K, static_cast<U>(0));
                    V_put[j] = std::max(this->K - S_ij, static_cast<U>(0));
                }

            if (!american)
            {
                for (int i {}; i < N ; ++i)
                {
                
                    for (int k {}; k < 2* N-2*i -1; ++k)
                    {
                        V_call[k] = p_u * V_call[k+2] + p_m * V_call[k+1] + p_d * V_call[k];
                        V_put[k] = p_u * V_put[k+2] + p_m * V_put[k+1] + p_d * V_put[k];
                    }
                    V_call.pop_back();
                    V_call.pop_back();

                    V_put.pop_back();
                    V_put.pop_back();
                }   
            
                option.call = std::exp(-this->r * this->T) * V_call[0];
                option.put = std::exp(-this->r * this->T) * V_put[0];
                return option;
            }
           else
           {
                U exp_rdt {std::exp(-this->r * this->T / static_cast<U>(N))};
                U S_k {0};
                U continuation_call {0};
                U intrinsic_call {0};

                U continuation_put {0};
                U intrinsic_put {0};

                for (int i {}; i < N ; ++i)
                {
                
                    for (int k {}; k < 2* N-2*i -1; ++k)
                    {
                        S_k = this->S * std::pow(static_cast<U>(u) , static_cast<U>(k - N + i + 1));
                        continuation_call = exp_rdt * (p_u * V_call[k+2] + p_m * V_call[k+1] + p_d * V_call[k]);
                        intrinsic_call = std::max(S_k - this->K, static_cast<U>(0));
                        V_call[k] = std::max(continuation_call, intrinsic_call);

                        continuation_put = exp_rdt * (p_u * V_put[k+2] + p_m * V_put[k+1] + p_d * V_put[k]);
                        intrinsic_put = std::max(this->K - S_k, static_cast<U>(0));
                        V_put[k] = std::max(continuation_put, intrinsic_put);
                    }

                    V_call.pop_back();
                    V_call.pop_back();

                    V_put.pop_back();
                    V_put.pop_back();
                }   
            
                option.call = V_call[0];
                option.put = V_put[0];
                return option;
            } 
        }

        Options<U,U> finiteDifferencesBackwardEuler(const int M, const int N, const U multiplier) const
        {
            if (N <= 0 || M <= 0 || multiplier <= 0)
            {
                throw std::invalid_argument("Invalid arguments\n");
            }

            Options<U,U> option {};

            U S_max {multiplier * std::max(this->S, this->K)};
            U dS {S_max / static_cast<U>(M)};
            U dT {this->T / static_cast<U>(N)};
            std::vector<U> S (M+1);

            for (int i{}; i < M + 1; ++i)
            {
                S[i] = static_cast<U>(i) * dS;
            }

            std::vector<U> V_call (M + 1);
            for (int i{}; i < M + 1; ++i)
            {
                V_call[i] = std::max(S[i] - this->K, static_cast<U>(0));
            }

            std::vector<U> V_put (M + 1);
            for (int i{}; i < M + 1; ++i)
            {
                V_put[i] = std::max(this->K - S[i], static_cast<U>(0));
            }

            FDMCoefficients<U> coeff {M};
            coeff.generateCoefficients(*this, dS, M);

            Structures::TridiagonalMatrix matrix {coeff.generateDiagonals(static_cast<U>(2) * dT, M)};
    
            for (int n {0}; n < N; ++n)
            {
                std::vector<U> RHS_call {coeff.generateRHSCallBE(*this, V_call, n, M, dT, S_max)};
                std::vector<U> *V_temp_call = new std::vector{Algorithms::thomasAlgorithm(matrix, RHS_call, M)};
                std::copy((*V_temp_call).begin(), (*V_temp_call).end(), V_call.begin()+1);
                delete V_temp_call;
                V_call[0] = static_cast<U>(0);
                V_call[M] = S_max * std::exp(-this->q * (static_cast<U>(n + 1)) * dT) - this->K * std::exp(-this->r* dT *(static_cast<U>(n+1)));

                std::vector<U> RHS_put {coeff.generateRHSPutBE(*this, V_put, n, M, dT, S_max)};
                std::vector<U> *V_temp_put = new std::vector{Algorithms::thomasAlgorithm(matrix, RHS_put, M)};
                std::copy((*V_temp_put).begin(), (*V_temp_put).end(), V_put.begin()+1);
                delete V_temp_put;
                V_put[0] = this->K * std::exp(-this->r * dT * (static_cast<U>(n+1)));
                V_put[M] = static_cast<U>(0);
            }

            int j {static_cast<int>(this->S / dS)};
            j = std::max(static_cast<int>(0), static_cast<int>(std::min(M-1,static_cast<int>(j))));
            option.call = V_call[j] + (this->S - S[j]) * (V_call[j+1] - V_call[j]) / (S[j+1] - S[j]);
            option.put = V_put[j] + (this->S - S[j]) * (V_put[j+1] - V_put[j]) / (S[j+1] - S[j]);
            return option;
        }


        Options<U,U> finiteDifferencesCrankNicolson(const int M, const int N, const U multiplier, bool american=false, int rannacher_smoothing=2) const
        {
            if (N <= 2 || M <= 2 || multiplier <= 0 || rannacher_smoothing > N)
            {
                throw std::invalid_argument("Invalid arguments\n");
            }

            Options<U,U> option {};

            U S_max {multiplier * std::max(this->S, this->K)};
            U dS {S_max / static_cast<U>(M)};
            U dT {this->T / static_cast<U>(N)};
            std::vector<U> S (M+1);

            for (int i{}; i < M + 1; ++i)
            {
                S[i] = static_cast<U>(i) * dS;
            }

            std::vector<U> V_call (M + 1);
            for (int i{}; i < M + 1; ++i)
            {
                V_call[i] = std::max(S[i] - this->K, static_cast<U>(0));
            }

            std::vector<U> V_put (M + 1);
            for (int i{}; i < M + 1; ++i)
            {
                V_put[i] = std::max(this->K - S[i], static_cast<U>(0));
            }

            FDMCoefficients<U> coeff {M};
            coeff.generateCoefficients(*this, dS, M);

            Structures::TridiagonalMatrix matrix {coeff.generateDiagonals(dT, M)};
    
            for (int n {0}; n < 2 * rannacher_smoothing; ++n)
            {
                std::vector<U> RHS_call {american ? coeff.generateRHSCallBEAmerican(*this, V_call, n, M, dT / static_cast<U>(2) , S_max) : coeff.generateRHSCallBE(*this, V_call, n, M, dT / static_cast<U>(2) , S_max)};
                std::vector<U> *V_temp_call = new std::vector{american ? Algorithms::thomasAlgorithmAmericanCall(matrix, this->K, RHS_call, M, dS) : Algorithms::thomasAlgorithm(matrix, RHS_call, M)};
                std::copy((*V_temp_call).begin(), (*V_temp_call).end(), V_call.begin()+1);
                delete V_temp_call;
                V_call[0] = static_cast<U>(0);
                V_call[M] = american ? std::max(S_max - this->K, S_max * std::exp(-this->q * (static_cast<U>(n + 1) / static_cast<U>(2)) * dT) - this->K * std::exp(-this->r* dT * (static_cast<U>(n+1) / static_cast<U>(2)))) : S_max * std::exp(-this->q * (static_cast<U>(n + 1) / static_cast<U>(2)) * dT) - this->K * std::exp(-this->r* dT * (static_cast<U>(n+1) / static_cast<U>(2)));

                std::vector<U> RHS_put {american ? coeff.generateRHSPutBEAmerican(*this, V_put, n, M, dT / static_cast<U>(2), S_max) : coeff.generateRHSPutBE(*this, V_put, n, M, dT / static_cast<U>(2) , S_max)};
                std::vector<U> *V_temp_put = new std::vector{american ? Algorithms::thomasAlgorithmAmericanPut(matrix, this->K, RHS_put, M, dS) : Algorithms::thomasAlgorithm(matrix, RHS_put, M)};
                std::copy((*V_temp_put).begin(), (*V_temp_put).end(), V_put.begin()+1);
                delete V_temp_put;
                V_put[0] = american ? this->K : this->K * std::exp(-this->r * dT * (static_cast<U>(n + 1) / static_cast<U>(2)));
                V_put[M] = std::max(this->K - S_max, static_cast<U>(0));
            }

            for (int n {rannacher_smoothing}; n < N; ++n)
            {
                std::vector<U> RHS_call {american ? coeff.generateRHSCallCNAmerican(*this, V_call, n, M, dT , S_max) : coeff.generateRHSCallCN(*this, V_call, n, M, dT, S_max)};
                std::vector<U> *V_temp_call = new std::vector{american ? Algorithms::thomasAlgorithmAmericanCall(matrix, this->K, RHS_call, M, dS) : Algorithms::thomasAlgorithm(matrix, RHS_call, M)};
                std::copy((*V_temp_call).begin(), (*V_temp_call).end(), V_call.begin()+1);
                delete V_temp_call;
                V_call[0] = static_cast<U>(0);
                V_call[M] = american ? std::max(S_max - this->K, S_max * std::exp(-this->q * (static_cast<U>(n + 1) * dT) - this->K * std::exp(-this->r* dT * (static_cast<U>(n+1))))) : S_max * std::exp(-this->q * (static_cast<U>(n + 1)) * dT) - this->K * std::exp(-this->r* dT * (static_cast<U>(n+1) / static_cast<U>(2)));

                std::vector<U> RHS_put {american ? coeff.generateRHSPutCNAmerican(*this, V_put, n, M, dT, S_max) : coeff.generateRHSPutCN(*this, V_put, n, M, dT, S_max)};
                std::vector<U> *V_temp_put = new std::vector{american ? Algorithms::thomasAlgorithmAmericanPut(matrix, this->K, RHS_put, M, dS) : Algorithms::thomasAlgorithm(matrix, RHS_put, M)};
                std::copy((*V_temp_put).begin(), (*V_temp_put).end(), V_put.begin()+1);
                delete V_temp_put;
                V_put[0] = american ? this->K : this->K * std::exp(-this->r * dT * (static_cast<U>(n+1)));
                V_put[M] = std::max(this->K - S_max, static_cast<U>(0));
            }

            int j {static_cast<int>(this->S / dS)};
            j = std::max(static_cast<int>(0), static_cast<int>(std::min(M-1,static_cast<int>(j))));
            option.call = V_call[j] + (this->S - S[j]) * (V_call[j+1] - V_call[j]) / (S[j+1] - S[j]);
            option.put = V_put[j] + (this->S - S[j]) * (V_put[j+1] - V_put[j]) / (S[j+1] - S[j]);
            return option;

        }
    };
        

template <typename S_t, typename K_t, typename r_t, typename sig_t, typename T_t, typename q_t>
GBM(S_t, K_t, r_t, sig_t, T_t, q_t) -> GBM<S_t, K_t, r_t, sig_t, T_t, q_t>;

template <typename S_t, typename K_t, typename r_t, typename sig_t, typename T_t>
GBM(S_t, K_t, r_t, sig_t, T_t) -> GBM<S_t, K_t, r_t, sig_t, T_t>;


int main()
{
    GBM optionData {100.0, 100.0, 0.01, 0.3, 1.0, 0.1};

    Options BS {optionData.closedForm()};
    Options ABT {optionData.monteCarlo(1000000, true, true, 42)};
    Options ATT {optionData.monteCarlo(1000000, true, false, 42)};
    Options ANT {optionData.monteCarlo(1000000, false, false, 42)};


     // Options AFDM {optionData.finiteDifferencesCrankNicolson(10000, 10000, 4, true)};

    std::cout << std::setprecision(15);

    std::cout << "Black Scholes Call: " << BS.call << '\n';
    std::cout << "Black Scholes Put: "  << BS.put << '\n';
    std::cout << "AV + CV Call: " << ABT.call << '\n';
    std::cout << "AV + CV Put: "  << ABT.put << '\n';
    std::cout << "AV Call: " << ATT.call << '\n';
    std::cout << "AV Put: " << ATT.put << '\n';
    std::cout << "None Call: " << ANT.call << '\n';
    std::cout << "None Put: " << ANT.put << '\n';

 //   std::cout << AFDM.call << '\n';
  //  std::cout << AFDM.put << '\n';

        

    
    return 0;
}







