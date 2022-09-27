#include <iostream>
#include "src/interpolation/newton/NewtonInterpolation.hpp"
#include "src/interpolation/cubic_spline/CubicSpline.hpp"

using namespace Interpolation;
int main() {
    constexpr std::size_t N = 10;
    std::array<scalar, N + 1> fx = {0, 0.033, 0.067, 0.100, 0.134, 0.168, 0.203, 0.238, 0.273, 0.309, 0.346};
    std::array<scalar, N + 1> xs;
    for (int i = 0; i < xs.size(); i++) xs[i] = static_cast<scalar>(i) / 10;

    auto c = NewtonInterpolation::newton_differencies<N>(xs, fx);
    scalar x = 0.95;
    scalar res = NewtonInterpolation::newton_interpolation<N>(x, c, xs);

    CubicSpline<N> I(xs, fx);
    int a = 1;

    std::cout << I.interpolate(x) << std::endl;
}
