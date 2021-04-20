#include <algo_opt/alias.hpp>
#include <algo_opt/local_descent.hpp>
#include <algo_opt/test_functions.hpp>

#include <boost/tuple/tuple.hpp>
#include <gnuplot-iostream/gnuplot-iostream.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <iostream>
#include <sstream>
#include <utility>
#include <vector>

namespace gio = gnuplotio;

void backtracking_line_search_viz(gio::PlotGroup *plots,
                                  ao::binary_f_t f,
                                  ao::binary_grad_t df,
                                  const Eigen::Vector2d &x0,
                                  double alpha)
{
  const auto n = 10;
  Eigen::Vector2d xi = x0;
  for (auto i = 0; i < n; ++i)
  {
    Eigen::Vector2d d = -1.0 * df(xi);
    d.normalize();

    Eigen::Vector2d x_full = xi + alpha * d;
    auto full_line = ao::point_vec_2d_t();
    full_line.emplace_back(xi(0), xi(1));
    full_line.emplace_back(x_full(0), x_full(1));
    plots->add_plot1d(full_line, "with lines lt rgb 'red' notitle");

    Eigen::Vector2d x_alpha = ao::backtracking_line_search(f, df, xi, d, alpha);
    auto alpha_line = ao::point_vec_2d_t();
    alpha_line.emplace_back(xi(0), xi(1));
    alpha_line.emplace_back(x_alpha(0), x_alpha(1));
    plots->add_plot1d(alpha_line, "with lines lt rgb 'black' notitle");

    xi = x_alpha;
  }
}

void strong_backtracking_line_search_viz(gio::PlotGroup *plots,
                                         ao::binary_f_t f,
                                         ao::binary_grad_t df,
                                         const Eigen::Vector2d &x0,
                                         double alpha)
{
  const auto n = 10;
  Eigen::Vector2d xi = x0;
  for (auto i = 0; i < n; ++i)
  {
    Eigen::Vector2d d = -1.0 * df(xi);
    d.normalize();

    Eigen::Vector2d x_full = xi + alpha * d;
    auto full_line = ao::point_vec_2d_t();
    full_line.emplace_back(xi(0), xi(1));
    full_line.emplace_back(x_full(0), x_full(1));
    plots->add_plot1d(full_line, "with lines lt rgb 'red' notitle");

    Eigen::Vector2d x_alpha =
        ao::strong_backtracking_line_search(f, df, xi, d, alpha);
    auto alpha_line = ao::point_vec_2d_t();
    alpha_line.emplace_back(xi(0), xi(1));
    alpha_line.emplace_back(x_alpha(0), x_alpha(1));
    plots->add_plot1d(alpha_line, "with lines lt rgb 'black' notitle");

    xi = x_alpha;
  }
}

int main(int, char **)
{
  auto gp = Gnuplot();
  auto plots = gp.plotGroup();

  const auto x_lim = ao::bracket_t<2>({-2.0, 2.0});
  const auto y_lim = x_lim;
  const auto levels =
      std::vector<double>({0.1, 0.5, 1, 2, 4, 8, 16, 32, 64, 128, 256});
  const auto d = 0.01;
  const auto a = 1.0;
  const auto b = 10.0;
  ao::rosenbrock_viz(&plots, x_lim, y_lim, levels, d, a, b);

  auto f = ao::make_rosenbrock(a, b);
  auto df = ao::make_grad_rosenbrock(a, b);
  Eigen::Vector2d xi(-1.5, -1.5);
  double alpha = 1.0;
  // backtracking_line_search_viz(&plots, f, df, xi, alpha);
  strong_backtracking_line_search_viz(&plots, f, df, xi, alpha);

  gp << "set xrange [" << x_lim[0] << ":" << x_lim[1] << "]\n";
  gp << "set yrange [" << y_lim[0] << ":" << y_lim[1] << "]\n";
  gp << plots;

  return 0;
}