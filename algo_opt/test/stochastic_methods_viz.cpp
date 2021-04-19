#include <algo_opt/alias.hpp>
#include <algo_opt/local_descent.hpp>
#include <algo_opt/stochastic_methods.hpp>
#include <algo_opt/test_functions.hpp>

#include <gnuplot-iostream/gnuplot-iostream.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <iostream>
#include <sstream>
#include <string>
#include <utility>
#include <vector>

namespace gio = gnuplotio;

void point_log_viz(gio::PlotGroup *plots,
                   const std::vector<ao::Vectord<2>> &log,
                   bool print = false)
{
  auto points = ao::point_vec_2d_t();
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
  ao::rosenbrock_viz(&plots, x_lim, y_lim, levels, d, a, b);

  auto f = ao::make_rosenbrock(a, b);
  ao::Vectord<2> x0(-1.5, 1.5);

  auto log = std::vector<ao::Vectord<2>>();
  ao::descent::mesh_adaptive_direct_search(f, x0, 1e-4, &log);

  point_log_viz(&plots, log);

  gp << "set xrange [" << x_lim[0] << ":" << x_lim[1] << "]\n";
  gp << "set yrange [" << y_lim[0] << ":" << y_lim[1] << "]\n";
  gp << plots;

  return 0;
}