#pragma once
#include <vector>
#include "Structures.hpp"

class Algorithms
{
    public:
        template <typename U>
        static std::vector<U> thomasAlgorithm(Structures::TridiagonalMatrix<U> matrix, std::vector<U> RHS, const int M)
        {
            // following algorithm from cfd-online.com
            std::vector<U> x (M-1);
            // forward elimination
            for (int i{1}; i < M-1; ++i)
            {
                U m {matrix.lower[i-1] / matrix.middle[i-1]};
                matrix.middle[i] -= m * matrix.upper[i-1];
                RHS[i] -= m * RHS[i-1];
            }

            x[M-2] = RHS[M-2] / matrix.middle[M-2];
            // backward substitution
            for (int i {M-3}; i >=0; --i)
            {
                x[i] = (RHS[i] - matrix.upper[i] * x[i+1]) / matrix.middle[i];
            }

            return x;
        }

        template <typename U>
        static std::vector<U> thomasAlgorithmAmericanCall(Structures::TridiagonalMatrix<U> matrix, const U& K, std::vector<U> RHS, const int M, const U dS)
        {
            // following algorithm from cfd-online.com
            std::vector<U> x (M-1);
            // backward elimination
            for (int i{M-2}; i >= 1; --i)
            {
                U m {matrix.upper[i-1] / matrix.middle[i]};
                matrix.middle[i-1] -= m * matrix.lower[i-1];
                RHS[i-1] -= m * RHS[i];
            }

            x[0] = std::max(dS - K, RHS[0] / matrix.middle[0]);
            // forward substitution
            for (int i {1}; i < M-1; ++i)
            {
                U S_i {static_cast<U>(i+1) * dS};
                x[i] = std::max(S_i - K, (RHS[i] - matrix.lower[i-1] * x[i-1]) / matrix.middle[i]);
            }
            return x;
        }

        template <typename U>
        static std::vector<U> thomasAlgorithmAmericanPut(Structures::TridiagonalMatrix<U> matrix, const U& K, std::vector<U> RHS, const int M, const U dS)
        {
            // following algorithm from cfd-online.com
            std::vector<U> x (M-1);
            // forward elimination
            for (int i{1}; i < M-1; ++i)
            {
                U m {matrix.lower[i-1] / matrix.middle[i-1]};
                matrix.middle[i] -= m * matrix.upper[i-1];
                RHS[i] -= m * RHS[i-1];
            }

            x[M-2] = std::max(K - static_cast<U>(M-1) * dS, RHS[M-2] / matrix.middle[M-2]);
            // backward substitution
            for (int i {M-3}; i >=0; --i)
            {
                U S_i {static_cast<U>(i+1) * dS};
                x[i] = std::max(K - S_i, (RHS[i] - matrix.upper[i] * x[i+1]) / matrix.middle[i]);
            }
            return x;
        }

};