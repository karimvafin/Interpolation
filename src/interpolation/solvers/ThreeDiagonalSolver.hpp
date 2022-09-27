//
// Created by Карим Вафин on 27.09.2022.
//

#include "interpolation/utils/ThreeDiagonalMatrix.hpp"

#ifndef INTERPOLATION_THREEDIAGONALSOLVER_HPP
#define INTERPOLATION_THREEDIAGONALSOLVER_HPP

namespace Interpolation {

    template <std::size_t N>
    std::array<scalar, N> solve_three_diagonal_matrix(const Utils::ThreeDiagonalMatrix<N>& matrix, const std::array<scalar, N>& col) {
        std::array<scalar, N> result;

        std::array<std::array<scalar, 2>, N - 1> params;

        params[0][0] = - matrix(0, 2) / matrix(0, 1);
        params[0][1] = col[0] / matrix(0, 1);

        for (int i = 1; i < col.size() - 1; i++)
        {
            params[i][0] = - matrix(i, 2) / (matrix(i, 0) * params[i - 1][0] + matrix(i, 1));
            params[i][1] = (col[i] - matrix(i, 0) * params[i - 1][1]) /
                           (matrix(i, 0) * params[i - 1][0] + matrix(i, 1));
        }

        result[N - 1] = (col[N - 1] - matrix(N - 1, 0) * params[N - 2][1]) /
                        (matrix(N - 1, 0) * params[N - 2][0] + matrix(N - 1, 1));

        for (int i = N - 2; i >= 0; i--)
        {
            result[i] = params[i][0] * result[i + 1] + params[i][1];
        }

        return result;
    }

} // Interpolation

#endif //INTERPOLATION_THREEDIAGONALSOLVER_HPP
