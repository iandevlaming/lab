
#include <search/dijkstra.hpp>

#include <algorithm>
#include <array>
#include <cmath>
#include <iostream>
#include <sstream>
#include <vector>

int main(int, char **)
{
  // construct a simple 2d grid
  using I = int;
  using N = std::array<I, 2>;
  using C = std::vector<N>;

  I lb = 0;
  I ub = 5;

  C graph;
  for (I x = lb; x < ub; ++x)
    for (I y = lb; y < ub; ++y)
      graph.push_back(N({x, y}));

  N source = {0, 0};

  auto length = [](const N *a, const N *b) {
    return std::abs((*a)[0] - (*b)[0]) + std::abs((*a)[1] - (*b)[1]);
  };

  auto neighbors = [&](const N *u) {
    auto n = std::vector<const N *>();
    I xu = (*u)[0];
    I yu = (*u)[1];
    for (I xi = xu - 1; xi <= xu + 1; ++xi)
    {
      for (I yi = yu - 1; yi <= yu + 1; ++yi)
      {
        if (xi >= lb && xi < ub && yi >= lb && yi < ub &&
            !(xi == xu && yi == yu))
        {
          N v_copy = {xi, yi};
          auto v_itr = std::ranges::find(graph, v_copy);
          n.push_back(&(*v_itr));
        }
      }
    }
    return n;
  };

  // auto [dist, prev] = search::dijsktra(graph, source, neighbors, length);
  auto [dist, prev] = search::dijkstra<I>(graph, source, neighbors, length);

  auto ss = std::ostringstream();
  for (I x = lb; x < ub; ++x)
  {
    for (I y = lb; y < ub; ++y)
    {
      N v_tmp = {x, y};
      N const *v_ptr = &(*std::ranges::find(graph, v_tmp));
      I d_v = dist[v_ptr];
      ss << d_v << "\t";
    }
    ss << "\n";
  }

  std::cout << ss.str() << std::endl;
}