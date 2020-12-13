#include <algo_opt/local_descent.hpp>
#include <algo_opt/stochastic_methods.hpp>

#include <gnuplot-iostream/gnuplot-iostream.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <iostream>
#include <sstream>
#include <string>
#include <utility>
#include <vector>

namespace ao = algo_opt;
namespace gio = gnuplotio;
using point_vec_2d_t = std::vector<std::tuple<double, double>>;

void rosenbrock_viz(gio::PlotGroup *plots,
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

  auto f_level = ao::makeLevelRosenbrock(a, b);
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

void point_log_viz(gio::PlotGroup *plots,
                   const std::vector<ao::Vectord<2>> &log,
                   bool print = false)
{
  auto points = point_vec_2d_t();
  for (unsigned int i = 0; i < log.size(); ++i)
  {
    ao::Vectord<2> xi = log[i];
    points.emplace_back(xi(0), xi(1));
    if (print)
      std::cout << "[" << xi(0) << ", " << xi(1) << "]" << std::endl;
  }

  plots->add_plot1d(points, "with lines lt rgb 'black' notitle");
}

int main(int, char **)
{
  auto gp = Gnuplot();
  auto plots = gp.plotGroup();

  const auto x_lim = ao::bracket_t<2>({-2.0, 2.0});
  const auto y_lim = x_lim;
  const auto levels =
      std::vector<double>({0.1, 0.5, 1, 2, 4, 8, 16, 32, 64, 128, 256});
  const auto d = 0.001;
  const auto a = 1.0;
  const auto b = 10.0;
  rosenbrock_viz(&plots, x_lim, y_lim, levels, d, a, b);

  auto f = ao::makeRosenbrock(a, b);
  ao::Vectord<2> x0(-1.5, 1.5);

  auto log = std::vector<ao::Vectord<2>>();
  ao::descent::mesh_adaptive_direct_search(f, x0, 1e-4, &log);

  point_log_viz(&plots, log);

  gp << "set xrange [" << x_lim[0] << ":" << x_lim[1] << "]\n";
  gp << "set yrange [" << y_lim[0] << ":" << y_lim[1] << "]\n";
  gp << plots;

  return 0;
}