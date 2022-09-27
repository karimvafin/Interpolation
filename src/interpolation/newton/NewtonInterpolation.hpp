//
// Created by Карим Вафин on 27.09.2022.
//

#ifndef INTERPOLATION_NEWTONINTERPOLATION_HPP
#define INTERPOLATION_NEWTONINTERPOLATION_HPP

#include <array>
#include "interpolation/types/Types.hpp"

namespace Interpolation {

    /** Количество узлов интерполяции N + 1 **/
    class NewtonInterpolation {
    public:
        template<std::size_t N>
        static constexpr std::array<scalar, N + 1>
        newton_differencies(const std::array<scalar, N + 1> &xs, const std::array<scalar, N + 1> &fx);

        template<std::size_t N>
        static constexpr scalar
        newton_interpolation(scalar x, const std::array<scalar, N + 1> &c, const std::array<scalar, N + 1> xs);
    };

    template<std::size_t N>
    constexpr std::array<scalar, N + 1> NewtonInterpolation::newton_differencies(const std::array<scalar, N + 1> &xs,
                                                                                 const std::array<scalar, N + 1> &fx) {
        std::array<scalar, N + 1> diffs(fx);
        std::array<scalar, N + 1> res;

        for (int diff_num = 0; diff_num < N; diff_num++) {
            res[diff_num] = diffs[0];
            for (int diff = 0; diff < N - diff_num; diff++)
                diffs[diff] = (diffs[diff] - diffs[diff + 1]) / (xs[diff] - xs[diff + diff_num + 1]);
        }
        res[N] = diffs[0];
        return res;
    }

    template<std::size_t N>
    constexpr scalar NewtonInterpolation::newton_interpolation(scalar x,
                                                               const std::array<scalar, N + 1> &c,
                                                               const std::array<scalar, N + 1> xs) {
        scalar res = c[0];
        for (int i = 0; i < N; i++) {
            scalar tmp = c[i + 1];
            for (int j = 0; j < i + 1; j++)
                tmp *= (x - xs[j]);
            res += tmp;
        }
        return res;
    }
}

#endif //INTERPOLATION_NEWTONINTERPOLATION_HPP
