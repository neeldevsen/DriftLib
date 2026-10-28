#pragma once


class Structures
{
    public:
        template <typename U>
        class TridiagonalMatrix
        {
            private:
                std::vector<U> lower {};
                std::vector<U> middle {};
                std::vector<U> upper {};

            public:
                explicit TridiagonalMatrix(const int M)
                : lower(M-1), middle(M), upper(M-1) {}

                friend class Algorithms;

                template <typename S_t, typename K_t, typename r_t, typename sig_t, typename T_t, typename q_t>
                friend class GBM;
        };
};
