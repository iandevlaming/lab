#include <algo_opt/alias.hpp>
#include <algo_opt/first_order.hpp>
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

template <typename P>
void step_viz(gio::PlotGroup *plots,
              ao::binary_f_t f,
              ao::binary_grad_t g,
              const ao::Vectord<2> &x0,
              P params,
              unsigned int n = 10,
              bool print = false)
{
  ao::Vectord<2> xi = x0;
  auto points = ao::point_vec_2d_t();
  points.emplace_back(xi(0), xi(1));

  for (unsigned int i = 0; i < n; ++i)
  {
    std::tie(xi, params) = ao::descent::step(params, f, g, xi);
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
  auto df = ao::make_grad_rosenbrock(a, b);
  Eigen::Vector2d xi(-1.5, -1.5);
  // auto params = ao::descent::GradientDescentParams(0.10);
  // auto params = ao::descent::GradientDescentOptParams(0.0001);
  // auto params = ao::descent::ConjugateGradientDescentParams(df, xi);
  // auto params = ao::descent::MomentumParams<2>();
  // auto params = ao::descent::NesterovMomentumParams<2>();
  // params.alpha = 0.01;
  // params.beta = 1.0;
  // auto params = ao::descent::AdagradParams<2>();
  // auto params = ao::descent::RMSPropParams<2>();
  // auto params = ao::descent::AdadeltaParams<2>();
  // auto params = ao::descent::AdamParams<2>();
  // auto params = ao::descent::HyperGradientDescentParams<2>();
  // auto params = ao::descent::DFPParams<2>();
  // auto params = ao::descent::BFGSParams<2>();
  // auto params = ao::descent::LimitedMemoryBFGSParams<2>();
  // params.m_max = 3;
  auto params =
      ao::descent::NoisyDescentParams<ao::unary_f_t,
                                      ao::descent::GradientDescentOptParams>();
  step_viz(&plots, f, df, xi, params, 10, true);

  gp << "set xrange [" << x_lim[0] << ":" << x_lim[1] << "]\n";
  gp << "set yrange [" << y_lim[0] << ":" << y_lim[1] << "]\n";
  gp << plots;

  return 0;
}