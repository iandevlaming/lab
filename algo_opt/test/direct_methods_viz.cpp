#include <algo_opt/direct_methods.hpp>
#include <algo_opt/local_descent.hpp>
#include <algo_opt/test_functions.hpp>

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

void pattern_log_viz(gio::PlotGroup *plots,
                     const std::vector<std::vector<ao::Vectord<2>>> &log)
{
  auto line_points = ao::point_vec_2d_t();
  auto pattern_points = ao::point_vec_2d_t();
  for (unsigned int i = 0; i < log.size(); ++i)
  {
    std::vector<ao::Vectord<2>> points = log[i];
    ao::Vectord<2> line_point = points.front();
    line_points.emplace_back(line_point(0), line_point(1));
    for (const auto &pattern_point : points)
      pattern_points.emplace_back(pattern_point(0), pattern_point(1));
  }

  plots->add_plot1d(line_points, "with lines lt rgb 'red' notitle");
  plots->add_plot1d(pattern_points, "with points lt rgb 'black' notitle");
}

void nelder_mead_log_viz(gio::PlotGroup *plots,
                         const std::vector<std::array<ao::Vectord<2>, 3>> &log)
{
  auto colors = std::vector<std::string>();
  colors.emplace_back("red");
  colors.emplace_back("green");
  colors.emplace_back("blue");
  colors.emplace_back("black");
  unsigned int i = 0;
  for (const auto &pattern : log)
  {
    auto points = ao::point_vec_2d_t();
    for (const auto &pt : pattern)
    {
      points.emplace_back(pt(0), pt(1));
    }
    points.emplace_back(pattern.front()(0), pattern.front()(1));
    auto ss = std::stringstream();
    ss << "with lines lt rgb '" << colors[i] << "' notitle";
    plots->add_plot1d(points, ss.str());
    i = (i + 1) % colors.size();
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
  const auto d = 0.001;
  const auto a = 1.0;
  const auto b = 10.0;
  ao::rosenbrock_viz(&plots, x_lim, y_lim, levels, d, a, b);

  auto f = ao::make_rosenbrock(a, b);
  ao::Vectord<2> x0(-1.5, 1.5);
  // auto point_log = std::vector<ao::Vectord<2>>();

  // ao::descent::cyclic_coordinate_descent(f, x0, 1e-4, &point_log);
  // ao::descent::accelerated_cyclic_coordinate_descent(f, x0, 1e-4,
  // &point_log); ao::descent::powell(f, x0, 1e-4, &point_log);
  // point_log_viz(&plots, point_log);

  // auto pattern_log = std::vector<std::vector<ao::Vectord<2>>>();

  // ao::descent::hooke_jeeves(f, x0, 0.5, 1e-4, 0.5, &pattern_log);
  // auto D = std::vector<ao::Vectord<2>>();
  // D.emplace_back(1.0, 0.0);
  // D.emplace_back(0.0, 1.0);
  // D.emplace_back(-sqrt(2.0) / 2.0, -sqrt(2.0) / 2.0);
  // ao::descent::generalized_pattern_search(
  //     f, x0, 0.5, D, 1e-4, 0.5, &pattern_log);
  // pattern_log_viz(&plots, pattern_log);

  auto nelder_mead_log = std::vector<std::array<ao::Vectord<2>, 3>>();
  ao::Vectord<2> x1(0.2, 0.2);
  ao::Vectord<2> x2(0.2, -0.2);
  ao::Vectord<2> x3(0.0, 0.0);
  std::array<ao::Vectord<2>, 3> S = {x1, x2, x3};
  ao::descent::nelder_mead(f, S, 1e-4, 1.0, 2.0, 0.5, &nelder_mead_log);

  nelder_mead_log_viz(&plots, nelder_mead_log);

  gp << "set xrange [" << x_lim[0] << ":" << x_lim[1] << "]\n";
  gp << "set yrange [" << y_lim[0] << ":" << y_lim[1] << "]\n";
  gp << plots;

  return 0;
}