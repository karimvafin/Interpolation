//
// Created by Карим Вафин on 27.09.2022.
//

#ifndef INTERPOLATION_CUBICSPLINE_HPP
#define INTERPOLATION_CUBICSPLINE_HPP

#include "interpolation/newton/NewtonInterpolation.hpp"
#include "interpolation/solvers/ThreeDiagonalSolver.hpp"

namespace Interpolation {

    template<std::size_t N>
    class CubicSpline {
        std::array<scalar, N + 1> xs_;
        std::array<std::array<scalar, 4>, N> C_;

        static constexpr std::array<std::array<scalar, 4>, N>
        build_spline_coeffs(const std::array<scalar, N + 1> &xs, const std::array<scalar, N + 1> &ys);

    public:
        constexpr CubicSpline(const std::array<scalar, N + 1> &xs, const std::array<scalar, N + 1> &ys) :
                xs_(xs),
                C_(build_spline_coeffs(xs, ys)) {}

        [[nodiscard]] constexpr scalar interpolate(scalar x) const;

    };

    template<std::size_t N>
    constexpr std::array<std::array<scalar, 4>, N>
    CubicSpline<N>::build_spline_coeffs(const std::array<scalar, N + 1> &xs,
                                        const std::array<scalar, N + 1> &ys) {
        std::array<std::array<scalar, 4>, N> res;
        std::array<scalar, N - 1> col;
        Utils::ThreeDiagonalMatrix<N - 1> matrix;

        for (int i = 0; i < N - 1; i++) {
            col[i] = 6 * ((ys[i] - ys[i + 1]) / (xs[i] - xs[i + 1]) - (ys[i + 1] - ys[i + 2]) / (xs[i + 1] - xs[i + 2])) /
                     (xs[i] - xs[i + 2]);
        }

        matrix.set_element(0, 0, static_cast<scalar>(0));
        matrix.set_element(0, 1, static_cast<scalar>(2));
        matrix.set_element(0, 2, (xs[2] - xs[1]) / (xs[2] - xs[0]));
        for (int i = 1; i < N - 2; i++) {
            matrix.set_element(i, 0, (xs[i] - xs[i - 1]) / (xs[i + 1] - xs[i - 1]));
            matrix.set_element(i, 1, static_cast<scalar>(2));
            matrix.set_element(i, 2, (xs[i + 1] - xs[i]) / (xs[i + 2] - xs[i]));
        }
        matrix.set_element(N - 2, 0, (xs[N - 1] - xs[N - 2]) / (xs[N] - xs[N - 2]));
        matrix.set_element(N - 2, 1, static_cast<scalar>(2));
        matrix.set_element(N - 2, 2, static_cast<scalar>(0));

        for (int i = 0; i < matrix.get_data().size() / 3; i++)
        {
            std::cout << matrix.get_data()[3 * i];
            std::cout << " ";
            std::cout << matrix.get_data()[3 * i + 1];
            std::cout << " ";
            std::cout << matrix.get_data()[3 * i + 2];
            std::cout << std::endl;
        }

        std::array<scalar, N - 1> c = solve_three_diagonal_matrix(matrix, col);

        res[0][0] = ys[1];
        res[0][1] = (ys[0] - ys[1]) / (xs[0] - xs[1]) + c[0] * (xs[1] - xs[0]) / 3;
        res[0][2] = c[0];
        res[0][3] = c[0] / (xs[1] - xs[0]);
        for (int i = 1; i < N - 1; i++) {
            res[i][0] = ys[i + 1];
            res[i][1] = (ys[i] - ys[i + 1]) / (xs[i] - xs[i + 1]) + c[i] * (xs[i + 1] - xs[i]) / 3 + c[i - 1] * (xs[i + 1] - xs[i]) / 6;
            res[i][2] = c[i];
            res[i][3] = (c[i] - c[i - 1]) / (xs[i + 1] - xs[i]);
        }

        res[N - 1][0] = ys[N];
        res[N - 1][1] = (ys[N - 1] - ys[N]) / (xs[N - 1] - xs[N]) + c[N - 2] * (xs[N] - xs[N - 1]) / 6;
        res[N - 1][2] = static_cast<scalar>(0);
        res[N - 1][3] = - c[N - 2] / (xs[N] - xs[N - 1]);

        return res;
    }

    template <std::size_t N>
    constexpr scalar CubicSpline<N>::interpolate(scalar x) const {
        scalar tmp = xs_[1];
        int i = 0;
        while (x > tmp) {
            i++;
            tmp = xs_[i + 1];
        }
        std::cout << i << std::endl;
        scalar dif = (x - xs_[i + 1]);
        return C_[i][0] + C_[i][1] * dif + C_[i][2] / 2 * dif * dif + C_[i][3] / 6 * dif * dif * dif;
    }

} // Interpolation
#endif //INTERPOLATION_CUBICSPLINE_HPP
