#pragma once

#include <algorithm>
#include <deque>
#include <iostream>
#include <map>
#include <optional>
#include <ranges>
#include <unordered_map>
#include <unordered_set>
#include <vector>

namespace algo_dm
{
template <template <typename, typename, typename...> typename MapT,
          typename Key,
          typename Value,
          typename... Args>
std::vector<Key> getKeys(const MapT<Key, Value, Args...>& map)
{
  auto keys_view = std::ranges::views::common(map | std::views::keys);
  return std::vector(keys_view.begin(), keys_view.end());
}

template <typename C>
bool isPermutation(const C& container_1, const C& container_2)
{
  auto is_in = [](const auto& cont_1) {
    return [&cont_1](const auto& elem_2) {
      return std::ranges::find(cont_1, elem_2) != cont_1.cend();
    };
  };

  return (container_1.size() == container_2.size()) &&
         std::ranges::all_of(container_1, is_in(container_2));
}

template <template <typename, typename...> typename C,
          typename T,
          typename... Ts>
bool contains(const C<T, Ts...>& container, const T& element)
{
  return std::ranges::find(container, element) != container.cend();
}

template <typename C>
C setDiff(const C& super_set, const C& disjoint_set)
{
  auto is_not_in = [](const auto& set_1) {
    return [&set_1](const auto& elem_2) {
      return std::ranges::find(set_1, elem_2) == set_1.cend();
    };
  };

  auto sub_set_view = super_set | std::views::filter(is_not_in(disjoint_set));
  auto sub_set = C(sub_set_view.begin(), sub_set_view.end());

  return sub_set;
}

template <typename C>
C setUnion(const C& a, const C& b)
{
  auto b_diff_a = setDiff(b, a); // preserve a ordering

  auto result = C(a.cbegin(), a.cend());
  std::ranges::move(b_diff_a, std::inserter(result, result.end()));

  return result;
}

template <typename MapT>
void FactorTable::checkInput(const MapT& table)
{
  if (!table.empty())
  {
    variable_names_ = table.cbegin()->first.getVariableNames();
    auto has_bad_keys = [&](const auto& p) {
      return p.first.getVariableNames() != variable_names_;
    };

    if (std::ranges::any_of(table, has_bad_keys))
      throw std::invalid_argument(
          "All Assignments in a FactorTable must have the same keys");
  }
}

template <typename NodeT, bool B>
std::unordered_set<NodeT> AdjacencyList<NodeT, B>::getNodes() const
{
  auto nodes_view = std::views::keys(graph_);
  auto nodes = std::unordered_set(nodes_view.begin(), nodes_view.end());

  auto values_view = std::views::values(graph_);
  for (const auto& value : values_view)
    nodes = setUnion(nodes, value);

  return nodes;
}

template <typename NodeT, bool B>
void AdjacencyList<NodeT, B>::addEdge(NodeT from_node, NodeT to_node)
{
  graph_[from_node].insert(to_node);
  if (!graph_.contains(to_node))
    graph_[to_node] = {};
}

template <typename NodeT, bool B>
const std::unordered_set<NodeT>&
AdjacencyList<NodeT, B>::getEdges(NodeT node) const
{
  return graph_.at(node);
}

template <typename NodeT>
std::unordered_set<NodeT> AdjacencyList<NodeT, false>::getNodes() const
{
  auto nodes_view = std::ranges::views::common(graph_ | std::views::keys);
  auto nodes = std::unordered_set(nodes_view.begin(), nodes_view.end());

  auto values_view = std::views::values(graph_);
  for (const auto& value : values_view)
    nodes = setUnion(nodes, value);

  return nodes;
}

template <typename NodeT>
void AdjacencyList<NodeT, false>::addEdge(const NodeT& from_node,
                                          const NodeT& to_node)
{
  graph_[from_node].insert(to_node);
  if (!graph_.contains(to_node))
    graph_[to_node] = {};
}

template <typename NodeT>
const std::unordered_set<NodeT>&
AdjacencyList<NodeT, false>::getEdges(const NodeT& node) const
{
  return graph_.at(node);
}

template <typename NodeT>
std::vector<NodeT> topoSort(const AdjacencyList<NodeT>& graph)
{
  auto nodes = graph.getNodes();

  auto num_parents = std::unordered_map<NodeT, int>();
  for (const auto& parent : nodes)
    for (const auto& child : graph.getEdges(parent))
      ++num_parents[child];

  auto is_parentless = [&num_parents](const auto& node) {
    return !num_parents.contains(node);
  };
  auto parentless = nodes | std::views::filter(is_parentless);
  auto node_queue = std::deque(parentless.begin(), parentless.end());
  auto topo_sort = std::vector<NodeT>();

  while (!node_queue.empty())
  {
    topo_sort.push_back(node_queue.front());
    node_queue.pop_front();
    const auto& cur_node = topo_sort.back();

    for (const auto& child : graph.getEdges(cur_node))
    {
      --num_parents[child];
      if (num_parents[child] == 0)
        node_queue.push_back(child);
    }
  }

  auto is_nonzero = [](const auto& p) { return p.second != 0; };
  if (std::ranges::any_of(num_parents, is_nonzero))
    throw std::invalid_argument("Topological is only valid DAGs");

  return topo_sort;
}

} // namespace algo_dm