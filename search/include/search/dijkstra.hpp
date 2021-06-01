#include <algorithm>
#include <functional>
#include <queue>
#include <tuple>
#include <type_traits>
#include <unordered_map>
#include <unordered_set>

namespace search
{
template <typename I,
          typename N,
          template <typename> typename C,
          typename NeighborF,
          typename LenF>
std::tuple<std::unordered_map<N const *, I>,
           std::unordered_map<N const *, N const *>>
dijkstra(C<N> const &graph, N const &source, NeighborF neighbors, LenF length)
{
  auto Q = std::unordered_set<const N *>();
  auto dist = std::unordered_map<N const *, I>();
  auto prev = std::unordered_map<N const *, N const *>();

  for (const auto &v : graph)
  {
    if (v == source)
      dist.insert({&v, I{}});
    else
      dist.insert({&v, std::numeric_limits<I>::max()});
    prev.insert({&v, nullptr});
    Q.insert(&v);
  }

  while (!Q.empty())
  {
    auto u_itr = std::ranges::min_element(
        Q, [&](const N *a, const N *b) { return dist[a] < dist[b]; });

    const N *u = *u_itr;
    I d_u = dist[u];

    Q.erase(u_itr);

    auto neighbors_u = neighbors(u);
    for (const N *v : neighbors_u)
    {
      I d_v = d_u + length(u, v);
      if (d_v < dist[v])
      {
        dist[v] = d_v;
        prev[v] = u;
      }
    }
  }

  return {dist, prev};
}
} // namespace search