#include <algo_opt/bracketing.hpp>

#include <boost/tuple/tuple.hpp>
#include <gnuplot-iostream/gnuplot-iostream.h>

#include <array>
#include <cmath>
#include <utility>
#include <vector>

namespace ao = algo_opt;

template <typename F, typename = std::enable_if<
                          std::is_invocable_r<double, F, double>::value>>
void bracket_viz(F f, const std::vector<std::vector<double>> &f_data,
                 const std::vector<std::array<double, 2>> &brackets) {
  auto i = 0;
  for (const auto &bracket : brackets) {
    const auto &[a, b] = bracket;
    auto a_pt = std::tuple(a, f(a));
    auto b_pt = std::tuple(b, f(b));

    auto gp = Gnuplot();
    gp << "plot '-' with lines, '-' with points title 'plot " << i++ << "'\n";
    gp.send1d(f_data);
    gp.send1d(std::vector({a_pt, b_pt}));
  }
}

template <typename F, typename = std::enable_if<
                          std::is_invocable_r<double, F, double>::value>>
void quadratic_fit_search_viz(
    F f, const std::vector<std::vector<double>> &f_data,
    const std::vector<ao::QuadraticFitSearchLog> &log) {
  auto i = 0;
  for (const auto &log_i : log) {
    const auto &[c0, c1, c2] = log_i.coefficients;
    auto curve_f = [&](double x) { return c0 + c1 * x + c2 * x * x; };
    auto curve = std::vector<std::vector<double>>();

    for (const auto &pt : f_data)
      curve.push_back(std::vector<double>({pt[0], curve_f(pt[0])}));

    auto bracket_i = std::vector<std::vector<double>>();
    for (const auto &pt : log_i.bracket)
      bracket_i.push_back(std::vector<double>({pt, f(pt)}));

    auto gp = Gnuplot();
    gp << "plot '-' with lines, '-' with lines, '-' with points title 'plot "
       << i++ << "'\n";
    gp.send1d(f_data);
    gp.send1d(curve);
    gp.send1d(bracket_i);
  }
}

int main(int, char **) {
  auto w = 2.0 * M_PI * 0.05;
  auto A = 10.0;
  auto f = [&w, &A](double x) { return A * -cos(w * x); };
  auto df = [&w, &A](double x) { return A * w * sin(w * x); };
  auto f_data = std::vector<std::vector<double>>();
  for (auto x = -10.0; x <= 10.0; x += 0.1)
    f_data.push_back(std::vector<double>({x, f(x)}));

  {
    // bracket minimum viz
    const auto x0 = -1.0;
    const auto s = 0.1;
    const auto k = 2.0;
    auto brackets = std::vector<std::array<double, 2>>();

    ao::bracket_minimum(f, x0, s, k, &brackets);
    bracket_viz(f, f_data, brackets);
  }

  {
    // quadratic fit viz
    auto bracket = std::array<double, 3>({-4.0, -3.0, 2.0});
    auto n = 3;
    auto log = std::vector<ao::QuadraticFitSearchLog>();

    ao::quadratic_fit_search(f, bracket, n, &log);
    quadratic_fit_search_viz(f, f_data, log);
  }

  {
    // bisection viz
    auto bracket = std::array<double, 2>({-3.0, 1.0});
    auto eps = 0.01;
    auto brackets = std::vector<std::array<double, 2>>();

    ao::bisection(df, bracket, eps, &brackets);
    bracket_viz(f, f_data, brackets);
  }
  return 1;
}