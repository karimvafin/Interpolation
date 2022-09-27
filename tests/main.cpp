#include <iostream>
#include "src/interpolation/newton/NewtonInterpolation.hpp"
#include "src/interpolation/cubic_spline/CubicSpline.hpp"
#include <cmath>

template <unsigned int N>
constexpr unsigned int factorial() {
    unsigned int res = 1;
    for (int i = 1; i <= N; i++)
        res *= i;
    return res;
}

using namespace Interpolation;
int main() {

    // интерполяция
    constexpr std::size_t N = 10;
    std::array<scalar, N + 1> fx = {0, 0.033, 0.067, 0.100, 0.134, 0.168, 0.203, 0.238, 0.273, 0.309, 0.346};
    std::array<scalar, N + 1> xs;
    for (int i = 0; i < xs.size(); i++) xs[i] = static_cast<scalar>(i) / 10;

    auto c = NewtonInterpolation::newton_differencies<N>(xs, fx);
    scalar x = 0.95;
    scalar res = NewtonInterpolation::newton_interpolation<N>(x, c, xs);

    std::cout << "Результат интерполяции полиномами: " << res << std::endl;

    CubicSpline<N> I(xs, fx);

    std::cout << "Результат интерполяции сплайном: " << I.interpolate(x) << std::endl;

    // ошибка интерполяции

    // интерполяции полиномами
    scalar f_11 = 15;
    scalar err_polynom = std::pow((xs[1] - xs[0]), N + 1) / (N + 1) * f_11;
    std::cout << "Погрешность интерполяции полинонами <= " << err_polynom << std::endl;

    // интерполяция сплайном
    scalar f_4 = 0.1;
    scalar err_spline = f_4 * std::pow((xs[1] - xs[0]), 4) / 24 / 16;
    std::cout << "Погрешность интерполяции cплайном <= " << err_spline << std::endl;
}
