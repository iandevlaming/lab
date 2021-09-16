#pragma once

#include <algorithm>
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
  auto select_keys = [](const auto& p) { return p.first; };
  auto keys_view = std::views::transform(map, select_keys);
  return std::vector<Key>(keys_view.begin(), keys_view.end());
}

template <template <typename...> typename C, typename... Ts>
bool isPermutation(const C<Ts...>& container_1, const C<Ts...>& container_2)
{
  auto is_in = [](const auto& cont_1) {
    return [&cont_1](const auto& elem_2) {
      return std::ranges::find(cont_1, elem_2) != cont_1.cend();
    };
  };

  return (container_1.size() == container_2.size()) &&
         std::ranges::all_of(container_1, is_in(container_2));
}

// TODO: it might be interesting to have the assignment templated on its
// variable names (will have to use a int id probably), so that the factor table
// could detect at compile time whether or not the assignments used in its
// construction are compatible

class Variable
{
public:
  using Name = std::string;

  Variable(const Name& name, int num_vals);
  Variable(Name&& name, int num_vals);
  const Name& getName() const;
  int getNumVals() const;

private:
  void init();

  Name name_;
  int num_vals_;
};

bool operator==(const Variable& lhs, const Variable& rhs);

std::ostream& operator<<(std::ostream& os, const Variable& var);

std::vector<Variable::Name> getNames(const std::vector<Variable>& variables);

class Assignment
{
public:
  using Key = Variable::Name;
  using Value = int;

  Assignment() = default;
  Assignment(const std::map<Key, Value>& assignment);
  void set(const Key& key, Value value);
  std::optional<Value> get(const Key& key) const;
  void erase(const Key& key);
  std::vector<Key> getVariableNames() const;
  bool operator==(const Assignment& other) const;

  struct Hash
  {
    size_t operator()(const Assignment& assignment) const;
  };

private:
  std::map<Key, Value> assignment_;
};

std::ostream& operator<<(std::ostream& os, const Assignment& a);

class FactorTable
{
public:
  using Key = Assignment;
  using Value = double;

  struct AssignmentCompare
  {
    bool operator()(const Assignment& lhs, const Assignment& rhs) const;
  };

  FactorTable() = default;
  FactorTable(const std::map<Assignment, double, AssignmentCompare>& table);
  FactorTable(
      const std::unordered_map<Assignment, double, Assignment::Hash>& table);
  FactorTable(const std::vector<Assignment>& assignments,
              const std::vector<double>& probabilities);
  std::optional<Value> get(const Key& key) const;
  void set(const Key& key, Value value);
  std::vector<Key> getAssignments() const;
  std::vector<Assignment::Key> getVariableNames() const;
  void normalize();

private:
  template <typename MapT>
  void checkInput(const MapT& table);
  std::map<Assignment, Value, AssignmentCompare> table_;
  std::vector<Assignment::Key> variable_names_;
};

bool operator==(const FactorTable& lhs, const FactorTable& rhs);

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

std::ostream& operator<<(std::ostream& os, const FactorTable& table);

class Factor
{
public:
  Factor(const std::vector<Variable>& variables, const FactorTable& table);
  Factor(std::vector<Variable>&& variables, FactorTable&& table);
  const std::vector<Variable>& getVariables() const;
  const FactorTable& getFactorTable() const;
  void normalize();

private:
  std::vector<Variable> variables_;
  FactorTable table_;
};

bool operator==(const Factor& lhs, const Factor& rhs);

std::ostream& operator<<(std::ostream& os, const Factor& f);

class AdjacencyList
{
public:
  AdjacencyList(int num_nodes);
  int getNumNodes() const;
  void addEdge(int from_node, int to_node);
  const std::unordered_set<int>& getEdges(int node) const;

private:
  int num_nodes_;
  std::vector<std::unordered_set<int>> graph_;
};

class BayesianNetwork
{
public:
  BayesianNetwork(const std::vector<Variable>& variables,
                  const std::vector<Factor>& factors,
                  const AdjacencyList& graph);
  const std::vector<Variable>& getVariables() const;
  const std::vector<Factor>& getFactors() const;
  const AdjacencyList& getGraph() const;

private:
  std::vector<Variable> variables_;
  std::vector<Factor> factors_;
  AdjacencyList graph_;
};

Assignment select(const Assignment& assignment,
                  const std::vector<Variable::Name>& variable_names);

std::vector<Variable> erase(const std::vector<Variable>& variables,
                            const Variable::Name& name);

std::vector<std::vector<int>>
cartesianProduct(std::vector<Variable>::const_iterator begin,
                 std::vector<Variable>::const_iterator end);
std::vector<std::vector<int>>
cartesianProduct(const std::vector<Variable>& variables);

std::vector<Assignment> assign(const std::vector<Variable>& variables);

double computeProbability(const std::vector<FactorTable>& tables,
                          const Assignment& assignment);
} // namespace algo_dm