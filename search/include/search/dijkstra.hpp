#include <algorithm>
#include <functional>
#include <queue>
#include <tuple>
#include <type_traits>
#include <unordered_map>
#include <unordered_set>

namespace search
{
template <typename X,
          typename IsGoalF,
          typename NeighborsF,
          typename Q = std::unordered_set<X>,
          // typename PushF = decltype(bfs_impl::push<X>),
          // typename PopF = decltype(bfs_impl::pop<X>),
          // typename IsEmptyF = decltype(bfs_impl::isEmpty<X>),
          // typename P = std::unordered_map<X, std::pair<X, int>>,
          // typename VisitF = decltype(search_impl::visit<X>),
          // typename IsVisitedF = decltype(search_impl::isVisited<X>),
          typename = has_signature<bool, IsGoalF, X>,
          typename = has_signature<std::vector<X>, NeighborsF, X>,
          // typename = has_signature<void, PushF, Q &, X>,
          // typename = has_signature<X, PopF, Q &>,
          // typename = has_signature<bool, IsEmptyF, Q const &>,
          // typename = has_signature<void, VisitF, P &, X, X>,
          // typename = has_signature<bool, IsVisitedF, P const &, X>>
          std::vector<X> dijkstra(X start,
                                  IsGoalF is_goal,
                                  NeighborsF neighbors,
                                  Q q = Q(),
                                  P p = P(),
                                  PushF push = bfs_impl::push<X>,
                                  PopF pop = bfs_impl::pop<X>,
                                  IsEmptyF is_empty = bfs_impl::isEmpty<X>,
                                  VisitF visit = search_impl::visit<X>,
                                  IsVisitedF is_visited =
                                      search_impl::isVisited<X>)
{
  auto dist = std::unordered_map<X, I>();
  auto push = [](Q &q_, X x_) { q_.insert(x_); };
  auto pop = [&dist](Q &q) {
    X x_ = std::ranges::min_element(
        q, [&](X const &a, X const &b) { return dist[a] < dist[b]; });
    q_.erase(x_);
    return x_;
  };
  auto is_empty = [](const Q &q_) { return q_.empty(); };
  auto visit = [&dist, &len](P &p_, X parent_, X child) {
    dist[child_] =
        dist.contains(parent) ? dist[parent_] + len(parent_, child_) : I{};
    search_impl::visit(p_, parent_, child_);
  };
  auto is_visited = &search_impl::is_visited;
  auto resolve_duplicates = [&dist, &len](Q &q_, P &p, X parent_, X child_) {
    I d_child = dist[parent_] + len(parent_, child);
    if (d_child < dist[child])
    {
      dist[child] = d_child;
      search_impl::visit(p_, parent_, child_);
    }
  };
  //////////////////////////////////////////////////////////////////////////////
  auto [goal, p_res] = search(
      q, p, start, push, pop, is_empty, is_goal, neighbors, visit, is_visited);

  if (!is_goal(goal))
    return {};

  return trace(p_res, goal);
}

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