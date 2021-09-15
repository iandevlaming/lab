#pragma once

#include <algorithm>
#include <map>
#include <optional>
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
  auto keys = std::vector<Key>();
  keys.reserve(map.size());

  auto select_keys = [](const auto& p) { return p.first; };
  std::ranges::transform(map, std::back_inserter(keys), select_keys);

  return keys;
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
  std::vector<Key> getKeys() const;
  bool operator==(const Assignment& other) const;

  struct Hash
  {
    size_t operator()(const Assignment& assignment) const;
  };

private:
  std::map<Key, Value> assignment_;
};

class FactorTable
{
public:
  using Key = Assignment;
  using Value = double;

  FactorTable() = default;
  FactorTable(
      const std::unordered_map<Assignment, double, Assignment::Hash>& table);
  FactorTable(const std::vector<Assignment>& assignments,
              const std::vector<double>& probabilities);
  std::optional<Value> get(const Key& key) const;
  void set(const Key& key, Value value);
  std::vector<Key> getKeys() const;
  std::vector<Assignment::Key> getAssignmentKeys() const;
  bool operator==(const FactorTable& rhs) const;
  void normalize();

private:
  std::unordered_map<Assignment, Value, Assignment::Hash> table_;
  std::vector<Assignment::Key> variable_names_;
};

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

Assignment select(const Assignment& assignment,
                  const std::vector<Variable::Name>& variable_names);

std::vector<std::vector<int>>
cartesianProduct(std::vector<Variable>::const_iterator begin,
                 std::vector<Variable>::const_iterator end);
std::vector<std::vector<int>>
cartesianProduct(const std::vector<Variable>& variables);

std::vector<Assignment> assign(const std::vector<Variable>& variables);

} // namespace algo_dm