
#include <search/breadth_first_search.hpp>
#include <search/depth_first_search.hpp>
#include <search/dijkstra.hpp>

#include <algorithm>
#include <array>
#include <cmath>
#include <iostream>
#include <numeric>
#include <sstream>
#include <vector>

int main(int, char **)
{
  // construct a simple 2d grid
  using I = int;
  using N = int;
  using C = std::vector<N>;

  I lb = 0;
  I ub = 5;
  I width = ub - lb;

  C graph(width * width);
  std::iota(graph.begin(), graph.end(), 0);

  auto id = [width](I x, I y) { return y * width + x; };
  auto coord = [width](N u) {
    I y = static_cast<int>(floor(u / width));
    I x = u % width;
    return std::array<I, 2>({x, y});
  };

  I x0 = 0;
  I y0 = 0;
  N source = id(x0, y0);

  auto dir = std::array<std::array<int, 2>, 4>({std::array<int, 2>{1, 0},
                                                std::array<int, 2>{0, 1},
                                                std::array<int, 2>{-1, 0},
                                                std::array<int, 2>{0, -1}});
  auto neighbors = [&](N u) -> std::vector<N> {
    auto n = std::vector<N>();
    auto [xu, yu] = coord(u);
    for (const auto &d : dir)
    {
      I xi = xu + d[0];
      I yi = yu + d[1];
      if (xi >= lb && xi < ub && yi >= lb && yi < ub && !(xi == xu && yi == yu))
        n.push_back(id(xi, yi));
    }
    return n;
  };

  // DIJKSTRA START
  // auto length = [](const N *a, const N *b) {
  //   return std::abs((*a)[0] - (*b)[0]) + std::abs((*a)[1] - (*b)[1]);
  // };

  // auto [dist, prev] = search::dijkstra<I>(graph, source, neighbors, length);

  // auto ss = std::ostringstream();
  // for (I x = lb; x < ub; ++x)
  // {
  //   for (I y = lb; y < ub; ++y)
  //   {
  //     N v_tmp = {x, y};
  //     N const *v_ptr = &(*std::ranges::find(graph, v_tmp));
  //     I d_v = dist[v_ptr];
  //     ss << d_v << "\t";
  //   }
  //   ss << "\n";
  // }

  // std::cout << ss.str() << std::endl;
  // DIJKSTRA END

  // BFS START
  // auto is_goal = [&](N const u) {
  //   auto [xi, yi] = coord(u);
  //   return abs(xi - x0) + abs(yi - y0) >= 4;
  // };
  // auto path = search::breadth_first_search(source, is_goal, neighbors);

  // if (path.empty())
  //   std::cout << "Error! path empty" << std::endl;
  // else
  //   for (auto const &point : path)
  //   {
  //     auto [x, y] = coord(point);
  //     std::cout << "[" << x << ", " << y << "]" << std::endl;
  //   }
  // BFS END

  // DFS START
  auto is_goal = [&](N const u) {
    auto [xi, yi] = coord(u);
    return abs(xi - x0) + abs(yi - y0) >= 4;
  };
  auto path = search::depth_first_search(source, is_goal, neighbors);

  if (path.empty())
    std::cout << "Error! path empty" << std::endl;
  else
    for (auto const &point : path)
    {
      auto [x, y] = coord(point);
      std::cout << "[" << x << ", " << y << "]" << std::endl;
    }
  // DFS END
}