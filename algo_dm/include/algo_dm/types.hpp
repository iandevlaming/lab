#pragma once

#include <boost/functional/hash.hpp>

#include <algorithm>
#include <map>
#include <numeric>
#include <optional>
#include <sstream>
#include <stdexcept>
#include <unordered_map>
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

  Variable(const Name& name, int num_vals) : name_(name), num_vals_(num_vals)
  {
    init();
  }
  Variable(Name&& name, int num_vals)
      : name_(std::move(name)), num_vals_(num_vals)
  {
    init();
  }
  const Name& getName() const { return name_; }
  int getNumVals() const { return num_vals_; }

private:
  void init()
  {
    if (num_vals_ <= 0)
    {
      auto error_stream = std::ostringstream();
      error_stream << "Variable num_vals must be positive, but got "
                   << num_vals_;
      throw std::domain_error(error_stream.str());
    }
  }

  Name name_;
  int num_vals_;
};

bool operator==(const Variable& lhs, const Variable& rhs)
{
  return lhs.getName() == rhs.getName() && lhs.getNumVals() == rhs.getNumVals();
}

std::vector<Variable::Name> getNames(const std::vector<Variable>& variables)
{
  auto variable_names = std::vector<Variable::Name>();
  variable_names.reserve(variables.size());

  auto select_names = [](const Variable& v) { return v.getName(); };
  std::ranges::transform(
      variables, std::back_inserter(variable_names), select_names);

  return variable_names;
}

class Assignment
{
public:
  using Key = Variable::Name;
  using Value = int;

  Assignment() = default;
  Assignment(const std::map<Key, Value>& assignment) : assignment_(assignment)
  {
  }
  Value& operator[](const Key& key) { return assignment_[key]; }
  const Value& at(const Key& key) const { return assignment_.at(key); }
  std::vector<Key> getKeys() const { return algo_dm::getKeys(assignment_); }
  bool operator==(const Assignment& other) const
  {
    return assignment_ == other.assignment_;
  }

  struct Hash
  {
    size_t operator()(const Assignment& assignment) const
    {
      auto seed = static_cast<size_t>(0);
      auto pairwise_hash = [&seed](const auto& p) {
        boost::hash_combine(seed, boost::hash_value(p.first));
        boost::hash_combine(seed, boost::hash_value(p.second));
      };

      std::ranges::for_each(assignment.assignment_, pairwise_hash);
      return seed;
    }
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
      const std::unordered_map<Assignment, double, Assignment::Hash>& table)
      : table_(table)
  {
    if (!table_.empty())
    {
      variable_names_ = table.cbegin()->first.getKeys();
      auto has_bad_keys = [&](const auto& p) {
        return p.first.getKeys() != variable_names_;
      };

      if (std::ranges::any_of(table_, has_bad_keys))
        throw std::invalid_argument(
            "All Assignments in a FactorTable must have the same keys");
    }
  }

  std::optional<Value> get(const Key& key) const
  {
    if (table_.contains(key))
      return table_.at(key);
    return {};
  }

  void set(const Key& key, Value value)
  {
    if (!table_.empty() && variable_names_ != key.getKeys())
      throw std::invalid_argument(
          "All Assignments in a FactorTable must have the same keys");

    if (table_.empty())
      variable_names_ = key.getKeys();

    table_[key] = value;
  }

  std::vector<Key> getKeys() const { return algo_dm::getKeys(table_); }
  std::vector<Assignment::Key> getAssignmentKeys() const
  {
    return variable_names_;
  }

  bool operator==(const FactorTable& rhs) const { return table_ == rhs.table_; }
  void normalize()
  {
    auto sum_value = [](auto sum, const auto& p) { return sum + p.second; };
    auto z = std::accumulate(table_.cbegin(), table_.cend(), 0.0, sum_value);
    std::ranges::for_each(table_, [&z](auto& p) { p.second /= z; });
  }

private:
  std::unordered_map<Assignment, Value, Assignment::Hash> table_;
  std::vector<Assignment::Key> variable_names_;
};

class Factor
{
public:
  Factor(const std::vector<Variable>& variables, const FactorTable& table)
      : variables_(variables), table_(table)
  {
  }
  Factor(std::vector<Variable>&& variables, FactorTable&& table)
      : variables_(std::move(variables)), table_(std::move(table))
  {
  }
  const std::vector<Variable>& getVariables() const { return variables_; }
  const FactorTable& getFactorTable() const { return table_; }
  void normalize() { table_.normalize(); }

private:
  std::vector<Variable> variables_;
  FactorTable table_;
};

Assignment select(const Assignment& assignment,
                  const std::vector<Variable::Name>& variable_names)
{
  auto assignment_keys = assignment.getKeys();
  auto is_invalid = [&assignment_keys](const Variable::Name& n) {
    return std::ranges::find(assignment_keys, n) == assignment_keys.cend();
  };
  if (std::ranges::any_of(variable_names, is_invalid))
    throw std::invalid_argument(
        "Cannot select variable names not contained by Assignment");

  auto sub_assignment = Assignment();
  auto sub_assign = [&assignment, &sub_assignment](const Variable::Name& name) {
    sub_assignment[name] = assignment.at(name);
  };
  std::ranges::for_each(variable_names, sub_assign);

  return sub_assignment;
}

// TODO: handle 0 num values
// TODO: clean up
std::vector<std::vector<int>>
cartesianProduct(std::vector<Variable>::const_iterator begin,
                 std::vector<Variable>::const_iterator end)
{
  if (begin == end)
    return {};

  auto product = std::vector<std::vector<int>>();
  auto sub_product = cartesianProduct(begin + 1, end);

  const auto& num_vals_i = begin->getNumVals();
  product.reserve(sub_product.empty() ? num_vals_i
                                      : sub_product.size() * num_vals_i);

  for (int i = 0; i < num_vals_i; ++i)
  {
    if (sub_product.empty())
    {
      auto product_i = std::vector<int>({i});
      product.push_back(product_i);
    }
    else
    {
      for (const auto& sub_product_j : sub_product)
      {
        auto product_i = std::vector<int>({i});
        std::ranges::copy(sub_product_j, std::back_inserter(product_i));
        product.push_back(product_i);
      }
    }
  }

  return product;
}

std::vector<std::vector<int>>
cartesianProduct(const std::vector<Variable>& variables)
{
  return cartesianProduct(variables.cbegin(), variables.cend());
}

std::vector<Assignment> assign(const std::vector<Variable>& variables)
{
  auto names = getNames(variables);
  auto values = cartesianProduct(variables);

  auto assignments = std::vector<Assignment>();
  assignments.reserve(values.size());

  for (const auto& value : values)
  {
    auto assignment = Assignment();
    for (auto i = 0; i < static_cast<int>(value.size()); ++i)
      assignment[names[i]] = value[i];
    assignments.push_back(assignment);
  }

  return assignments;
}
} // namespace algo_dm