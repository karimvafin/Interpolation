//
// Created by Карим Вафин on 27.09.2022.
//

#ifndef INTERPOLATION_THREEDIAGONALMATRIX_HPP
#define INTERPOLATION_THREEDIAGONALMATRIX_HPP

#include <array>
#include "interpolation/types/Types.hpp"

namespace Interpolation::Utils {

    template <std::size_t N>
    class ThreeDiagonalMatrix {
        std::array<scalar, 3 * N> data_;

    public:
        constexpr explicit ThreeDiagonalMatrix() = default;
        constexpr explicit ThreeDiagonalMatrix(const std::array<scalar, 3 * N>& data) : data_(data) {}

        [[nodiscard]] constexpr std::array<scalar, 3 * N> get_data() const;
        [[nodiscard]] constexpr scalar operator()(unsigned int i, unsigned int j) const;
        void set_element(unsigned int i, unsigned int j, scalar elem);
    };

    template <std::size_t N>
    constexpr std::array<scalar, 3 * N> ThreeDiagonalMatrix<N>::get_data() const {
        return data_;
    }

    /** Нумерация i и j с нуля до N - 1 **/
    template <std::size_t N>
    constexpr scalar ThreeDiagonalMatrix<N>::operator()(unsigned int i, unsigned int j) const {
        return data_[3 * i + j];
    }

    template <std::size_t N>
    void ThreeDiagonalMatrix<N>::set_element(unsigned int i, unsigned int j, scalar elem) {
        data_[3 * i + j] = elem;
    }
} // Interpolation::Utils

#endif //INTERPOLATION_THREEDIAGONALMATRIX_HPP
