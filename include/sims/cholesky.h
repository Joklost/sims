#ifndef MANETSIMS_CHOLESKY_H
#define MANETSIMS_CHOLESKY_H

#include <vector>
#include <future>
#include <cmath>
#include <algorithm>

#include "math.h"

#include <common/equality.h>

namespace sims {
    namespace cholesky {

        template<typename T>
        sims::linalg::vecvec<T> slow_cholesky(const sims::linalg::vecvec<T> matrix) {
            auto size = matrix.size();
            sims::linalg::vecvec<T> vec{};
            vec.resize(size, ::std::vector<T>(size));

            for (auto i = 0; i < size; ++i) {
                for (auto j = 0; j <= i; ++j) {
                    T sum = (T) 0;
                    for (auto k = 0; k < j; ++k) {
                        sum += vec[i][k] * vec[j][k];
                    }

                    if (i == j) {
                        auto val = ::std::sqrt(matrix[i][i] - sum);
                        vec[i][j] = (::std::isnan(val) ? 0.0 : val);
                    } else {
                        vec[i][j] = ((T) 1) / vec[j][j] * (matrix[i][j] - sum);
                    }

                }
            }

            return vec;
        }

        template<typename T>
        sims::linalg::vecvec<T> cholesky(const sims::linalg::vecvec<T> matrix) {
            auto size = matrix.size();
            sims::linalg::vecvec<T> vec{};
            vec.resize(size, ::std::vector<T>(size));

            for (auto row = 0; row < size; ++row) {
                for (auto column = 0; column <= row; ++column) {
                    T sum = (T) 0;
                    for (auto i = 0; i < column; ++i) {
                        sum += vec[row][i] * vec[column][i];
                    }

                    vec[row][column] = row == column ?
                                       ::std::sqrt(::std::abs(matrix[row][row] - sum)) :
                                       ((T) 1) / vec[column][column] * (matrix[row][column] - sum);
                }
            }

            return vec;
        }

        template<typename T>
        bool is_positive_definite(const sims::linalg::vecvec<T> &a) {
            auto size = a.size();

            for (auto r = 0; r < size; ++r) {
                for (auto c = 0; c <= r; ++c) {

                    if (c == r) {
                        T sum = (T) 0;
                        for (auto j = 0; j < c; j++) {
                            sum += a[r][j] * a[r][j];
                        }
                        if ((a[c][c] - sum) < 0) {
                            return false;
                        }
                    }
                }
            }

            return true;
        }


        template<typename T>
        sims::linalg::vecvec<T> iterative_spectral(sims::linalg::vecvec<T> &c, unsigned long max_it, double val) {
            auto size = c.size();
            auto cold{c};
            auto j = 0;
            for (j = 0; j < max_it; ++j) {
                for (auto i = 0; i < size; ++i) {
                    cold[i][i] = 1.0;
                }

                auto eigen = eig(cold, max_it);
                auto s = eigen.vectors;
                auto l = eigen.values;

                if (diag(l) > 0.0) {
                    break;
                }

                /* Replace negative eigenvalues with a small positive value. */
                ::std::replace_if(l.begin(), l.end(), [](T x) {
                    return x < 0.0;
                }, val);

                sims::linalg::vecvec<T> ll{l};
                cold = s * ll * transpose(s);
            }

            return cold;
        }


        template<typename T>
        sims::linalg::vecvec<T> adjusted_gradient(sims::linalg::vecvec<T> &c) {
            auto sigma = 0.5; // sigma, balance between speed and accuracy
            auto tol = 1e-2;  // stopping criteria
            auto steps = 100; // number of iterations

            auto n = c.size();

            auto b{c};

            for (auto i = 0; i < n; ++i) {
                b[i][i] = ::std::sqrt(n) / 2;
            }

            sims::linalg::vecvec<T> row_sums{};
            row_sums.resize(n, ::std::vector<T>(1));
            for (auto i = 0; i < n; ++i) {
                auto sum = ::std::accumulate(b[i].begin(), b[i].end(), (T) 0.0, [](T sum, T element) -> T {
                    return sum + ::std::pow(element, 2);
                });
                row_sums[i][0] = ((T) 1) / ::std::sqrt(sum);
            }

            b = row_sums * b;
            auto bb = b * b;
            sims::linalg::vecvec<T> bbc = bb - c;
            sims::linalg::vecvec<T> r{};

            return sims::linalg::vecvec<T>{};
        }


        template<typename T>
        sims::linalg::vecvec<T> cholesky_v2(const sims::linalg::vecvec<T> a) {
            auto size = a.size();
            sims::linalg::vecvec<T> ret{};
            ret.resize(size, ::std::vector<T>(size));

            for (auto r = 0; r < size; ++r) {
                for (auto c = 0; c <= r; ++c) {

                    if (c == r) {
                        T sum = (T) 0;
                        for (auto j = 0; j < c; j++) {
                            sum += ret[r][j] * ret[r][j];
                        }
                        if (sum > a[c][c])
                            sum = (T) 0;

                        auto val = a[c][c] - sum;
                        ret[c][c] = ::std::sqrt(val);
                    } else {
                        T sum = (T) 0;
                        for (auto j = 0; j < c; j++) {
                            sum += ret[r][j] * ret[c][j];
                        }

                        assert(!common::is_equal(ret[c][c], (T) 0));
                        ret[r][c] = ((T) 1) / ret[c][c] * (a[r][c] - sum);
                    }
                }
            }

            return ret;
        }

    }
}

#endif /* MANETSIMS_CHOLESKY_H */
