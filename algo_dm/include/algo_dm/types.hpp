#pragma once

#include <algorithm>
#include <iostream>
#include <map>
#include <optional>
#include <ranges>
#include <type_traits>
#include <unordered_map>
#include <unordered_set>
#include <vector>

namespace algo_dm
{
template <template <typename, typename, typename...> typename MapT,
          typename Key,
          typename Value,
          typename... Args>
std::vector<Key> getKeys(const MapT<Key, Value, Args...>& map);

template <typename C>
bool isPermutation(const C& container_1, const C& container_2);

template <template <typename, typename...> typename C,
          typename T,
          typename... Ts>
bool contains(const C<T, Ts...>& container, const T& element);

template <typename C>
C setDiff(const C& super_set, const C& disjoint_set);

template <typename C>
C setUnion(const C& a, const C& b);

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

std::ostream& operator<<(std::ostream& os, const FactorTable& table);

class Factor
{
public:
  Factor() = default;
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

Factor operator*(const Factor& lhs, const Factor& rhs);

bool isInScope(const Factor& factor, const Variable::Name name);

std::vector<Variable::Name> getNames(const std::vector<Factor>& factors);

template <typename NodeT, bool = std::is_integral_v<NodeT>>
class AdjacencyList
{
public:
  AdjacencyList() = default;
  std::unordered_set<NodeT> getNodes() const;
  void addEdge(NodeT from_node, NodeT to_node);
  const std::unordered_set<NodeT>& getEdges(NodeT node) const;

private:
  std::unordered_map<NodeT, std::unordered_set<NodeT>> graph_;
};

template <typename NodeT>
class AdjacencyList<NodeT, false>
{
public:
  AdjacencyList() = default;
  std::unordered_set<NodeT> getNodes() const;
  void addEdge(const NodeT& from_node, const NodeT& to_node);
  const std::unordered_set<NodeT>& getEdges(const NodeT& node) const;

private:
  std::unordered_map<NodeT, std::unordered_set<NodeT>> graph_;
};

template <typename NodeT>
std::vector<NodeT> topoSort(const AdjacencyList<NodeT>& graph);

class BayesianNetwork
{
public:
  BayesianNetwork() = default;
  BayesianNetwork(const std::unordered_map<Variable::Name, Factor>& nodes,
                  const AdjacencyList<Variable::Name>& graph);
  std::vector<Variable> getVariables() const;
  std::vector<Factor> getFactors() const;
  const std::unordered_map<Variable::Name, Factor> getNodes() const;
  const AdjacencyList<Variable::Name>& getGraph() const;
  const Factor getFactor(const Variable::Name& node) const;

private:
  std::unordered_map<Variable::Name, Factor> nodes_;
  AdjacencyList<Variable::Name> graph_;
};

Assignment select(const Assignment& assignment,
                  const std::vector<Variable::Name>& variable_names);

std::vector<Variable> select(const std::vector<Variable>& variables,
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

#include <algo_dm/inl/types.inl>