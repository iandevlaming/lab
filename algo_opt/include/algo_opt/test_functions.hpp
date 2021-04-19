#pragma once

#include <algo_opt/alias.hpp>

#include <gnuplot-iostream/gnuplot-iostream.h>

#include <array>
#include <cmath>
#include <optional>
#include <type_traits>

namespace algo_opt
{
binary_f_t make_rosenbrock(double a = 1, double b = 100)
{
  // f = (a - x)^2 + b * (y - x^2)^2
  return [a, b](const Vectord<2> &pt) {
    const auto &x = pt(0);
    const auto &y = pt(1);
    return std::pow((a - x), 2) + b * std::pow((y - x * x), 2);
  };
}

// auto diff (or even getting funclib working) would be a better move _but_ here
// we are
binary_grad_t make_grad_rosenbrock(double a = 1, double b = 100)
{
  // df/dx = -2 * (a - x) - 4 * x * b * (y - x^2)
  // df/dy = 2 * b * (y - x^2)
  return [a, b](const Vectord<2> &pt) -> Vectord<2> {
    const auto &x = pt(0);
    const auto &y = pt(1);
    const auto dfdx = -2.0 * (a - x) - 4.0 * x * b * (y - x * x);
    const auto dfdy = 2.0 * b * (y - x * x);
    Vectord<2> grad(dfdx, dfdy);
    return grad;
  };
}

binary_hess_t make_hess_rosenbrock(double a = 1, double b = 100)
{
  // d2f/dx2 = 2 - 4 * b * y + 12 * b * x^2
  // d2f/dxdy = -4 * b * x
  // d2f/dy2 = -4 * b * x
  // d2f/dydx = 2b
  return [a, b](const Vectord<2> &pt) -> Matrixd<2> {
    const auto &x = pt(0);
    const auto &y = pt(1);
    const auto d2fdx2 = 2.0 - 4.0 * b * y + 12.0 * b * x * x;
    const auto d2fdxdy = -4.0 * b * x;
    const auto d2fdydx = -4.0 * b * x;
    const auto d2fdy2 = 2.0 * b;
    // clang-format off
    Matrixd<2> hess;
    hess << d2fdx2,  d2fdxdy,
            d2fdydx, d2fdy2;
    // clang-format on
    return hess;
  };
}

binary_level_t make_level_rosenbrock(double a = 1, double b = 100)
{
  return [a, b](double x, double f) -> std::optional<std::array<double, 2>> {
    // y = sqrt((f - (a - x)^2) / b) + x^2
    auto g = (f - pow((a - x), 2)) / b;
    if (g < 0.0)
      return {};
    auto h = std::sqrt(g);
    auto x2 = x * x;
    return std::array<double, 2>({-h + x2, h + x2});
  };
}

void rosenbrock_viz(gnuplotio::PlotGroup *plots,
                    ao::bracket_t<2> x_lim,
                    ao::bracket_t<2> y_lim,
                    const std::vector<double> &levels,
                    double d = 0.1,
                    double a = 1,
                    double b = 100)
{
  std::ranges::sort(x_lim);
  std::ranges::sort(y_lim);
  auto x = std::vector<double>({x_lim[0]});
  while (x.back() < x_lim[1])
    x.push_back(x.back() + d);

  auto f_level = ao::make_level_rosenbrock(a, b);
  auto get_level_points = [&](auto &lower, auto &upper, const auto &li) {
    return [&](const auto &xi) {
      auto yi = f_level(xi, li);
      if (yi.has_value())
      {
        lower.emplace_back(xi, yi.value()[0]);
        upper.emplace_back(xi, yi.value()[1]);
      }
    };
  };

  for (const auto &level : levels)
  {
    auto xy_upper = point_vec_2d_t();
    auto xy_lower = point_vec_2d_t();
    std::for_each(
        x.cbegin(), x.cend(), get_level_points(xy_lower, xy_upper, level));

    auto points = xy_lower;
    std::copy(xy_upper.rbegin(), xy_upper.rend(), std::back_inserter(points));

    auto ss = std::stringstream();
    ss << "with lines notitle";
    plots->add_plot1d(points, ss.str());
  }
}

}