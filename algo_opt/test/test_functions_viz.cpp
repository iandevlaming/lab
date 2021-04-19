#include <algo_opt/alias.hpp>
#include <algo_opt/test_functions.hpp>

#include <gnuplot-iostream/gnuplot-iostream.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <iostream>
#include <sstream>
#include <utility>
#include <vector>

namespace gio = gnuplotio;

int main(int, char **)
{
  auto gp = Gnuplot();
  auto plots = gp.plotGroup();

  const auto x_lim = ao::bracket_t<2>({-4.0, 4.0});
  const auto y_lim = x_lim;
  // const auto levels =
  //     std::vector<double>({0.1, 0.5, 1, 2, 4, 8, 16, 32, 64, 128, 256});
  // const auto d = 0.01;
  // const auto a = 1.0;
  // const auto b = 10.0;
  // ao::rosenbrock_viz(&plots, x_lim, y_lim, levels, d, a, b);

  const auto levels =
      std::vector<double>({1.0, 2.0, 3.0, 4.0});
  ao::flower_viz(&plots, x_lim, y_lim, levels, 360, 1.0, 1.0, 4.0);

  gp << "set xrange [" << x_lim[0] << ":" << x_lim[1] << "]\n";
  gp << "set yrange [" << y_lim[0] << ":" << y_lim[1] << "]\n";
  gp << plots;

  return 0;
}