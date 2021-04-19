#pragma once

#include <algo_opt/alias.hpp>

#include <gnuplot-iostream/gnuplot-iostream.h>

#include <algorithm>
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
                    bracket_t<2> x_lim,
                    bracket_t<2> y_lim,
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

  auto f_level = make_level_rosenbrock(a, b);
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

binary_f_t make_flower(double a = 1.0, double b = 1.0, double c = 4.0)
{
  // f = a * sqrt(x^2 + y^2) + b * sin(c * atan2(y, x))
  return [a, b, c](const Vectord<2> &pt) {
    return a * pt.norm() + b * sin(c * atan2(pt(1), pt(0)));
  };
}

binary_grad_t make_grad_flower(double a = 1.0, double b = 1.0, double c = 4.0)
{
  // df/dx = a * x / sqrt(x^2 + y^2) + b * c * cos(c * atan2(y, x)) * [-y / (x^2 + y^2)]
  // df/dy = a * y / sqrt(x^2 + y^2) + b * c * cos(c * atan2(y, x)) * [x / (x^2 + y^2)]
  return [a, b, c](const Vectord<2> &pt) -> Vectord<2> {
    const auto &x = pt(0);
    const auto &y = pt(1);

    const auto norm_sq = pt.squaredNorm();
    const auto norm = std::sqrt(norm_sq);

    const auto d = a / norm;
    const auto e = b * c * cos(c * atan2(y, x)) / norm_sq;

    const auto dfdx = d * x - e * y;
    const auto dfdy = d * y + e * x;
    Vectord<2> grad(dfdx, dfdy);
    return grad;
  };
}

binary_hess_t make_hess_flower(double a = 1.0, double b = 1.0, double c = 4.0)
{
  // d2f/dx2 = a * [1 - x^2 / (x^2 + y^2)] / sqrt(x^2 + y^2) + b * c * y * [2 * x * cos(c * atan2(y, x)) - c * y * sin(c * atan2(y, x))] / (x^2 + y^2)^2
  // d2f/dxdy = -a * x * y * (x^2 + y^2)^(-3/2) + b * c / (x^2 + y^2) * [c * x * y * sin(c * atan2(y, x)) / (x^2 + y^2) + cos(c * atan2(y, x)) * (1 - 2 * x^2 / (x^2 + y^2))]
  // d2f/dy2 = a * [1 - y^2 / (x^2 + y^2)] / sqrt(x^2 + y^2) - b * c * x * [2 * y * cos(c * atan2(y, x)) + c * x * sin(c * atan2(y, x))] / (x^2 + y^2)^2
  // d2f/dydx = -a * x * y * (x^2 + y^2)^(-3/2) + b * c / (x^2 + y^2) * [c * x * y * sin(c * atan2(y, x)) / (x^2 + y^2) + cos(c * atan2(y, x)) * (2 * y^2 / (x^2 + y^2) - 1)]
  return [a, b, c](const Vectord<2> &pt) -> Matrixd<2> {
    const auto &x = pt(0);
    const auto &y = pt(1);

    const auto norm_sq = pt.squaredNorm();
    const auto norm = std::sqrt(norm_sq);

    const auto A = a / norm;
    const auto C = c * atan2(y, x);
    const auto bc = b * c;

    const auto y_hat = y / norm_sq;
    const auto y_y_hat = y * y_hat;
    const auto x_hat = x / norm_sq;
    const auto x_x_hat = x * x_hat;
    const auto x_y_hat = x * y_hat;

    const auto d2fdx2 = A * (1.0 - x_x_hat) + bc * y_hat * (2.0 * x * cos(C) - c * y * sin(C));
    const auto d2fdxdy = -A * x_y_hat + (bc / norm_sq) * (c * x_y_hat * sin(C) + cos(C) * (1.0 - 2.0 * x_x_hat));
    const auto d2fdydx = -A * x_y_hat + (bc / norm_sq) * (c * x_y_hat * sin(C) + cos(C) * (2.0 * y_y_hat - 1.0));
    const auto d2fdy2 = A * (1.0 - y_y_hat) - bc * x_hat * (2.0 * y * cos(C) + c * x * sin(C));
    // clang-format off
    Matrixd<2> hess;
    hess << d2fdx2,  d2fdxdy,
            d2fdydx, d2fdy2;
    // clang-format on
    return hess;
  };
}

point_vec_2d_t level_flower(double f, size_t N, double a = 1.0, double b = 1.0, double c = 4.0)
{
    // f = a * r + b * sin(c * theta)
    // r = (f - b * sin(c * theta)) / a

    auto pts = point_vec_2d_t();
    auto inc = 2.0 * M_PI / static_cast<double>(N);
    for (size_t i = 0; i <= N; ++i)
    {
      auto theta = -M_PI + inc * static_cast<double>(i);
      auto r = (f - b * sin(c * theta)) / a;

      auto x = r * cos(theta);
      auto y = r * sin(theta);

      pts.emplace_back(x, y);
    }

  return pts;
}

void flower_viz(gnuplotio::PlotGroup *plots,
                    bracket_t<2> x_lim,
                    bracket_t<2> y_lim,
                    const std::vector<double> &levels,
                    size_t N = 360,
                    double a = 1.0,
                    double b = 1.0,
                    double c = 4.0)
{
  std::ranges::sort(x_lim);
  std::ranges::sort(y_lim);

  auto isOutOfBounds = [x_lim, y_lim](auto const& pt){return std::get<0>(pt) < x_lim[0] || std::get<0>(pt) > x_lim[1] || std::get<1>(pt) < y_lim[0] || std::get<1>(pt) > y_lim[1];};
  for (const auto& level : levels)
  {
    auto level_pts = level_flower(level, N, a, b, c);
    level_pts.erase(std::remove_if(level_pts.begin(), level_pts.end(), isOutOfBounds), level_pts.end());

    auto ss = std::stringstream();
    ss << "with lines notitle";
    plots->add_plot1d(level_pts, ss.str());
  }
}
}