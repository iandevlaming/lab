#pragma once

#include <unordered_map>
#include <vector>

namespace algo_dm
{
// TODO: it might be interesting to have the assignment templated on its
// variable names (will have to use a int id probably), so that the factor table
// could detect at compile time whether or not the assignments used in its
// construction are compatible

class Variable
{
public:
  using Name = std::string;

  Variable(const Name &name, int num_vals) : name_(name), num_vals_(num_vals) {}
  Variable(Name &&name, int num_vals)
      : name_(std::move(name)), num_vals_(num_vals)
  {
  }
  const Name &name() const { return name_; }
  int numVals() const { return num_vals_; }

private:
  Name name_;
  int num_vals_;
};

using Assignment = std::unordered_map<Variable::Name, int>;
using FactorTable = std::unordered_map<Assignment, double>;

class Factor
{
public:
  Factor(const std::vector<Variable> &variables, const FactorTable &table)
      : variables_(variable), table_(table)
  {
  }
  Factor(std::vector<Variable> &&variables, FactorTable &&table)
      : variables_(std::move(variable)), table_(std::move(table))
  {
  }
  const std::vector<Variable> &variables() const { return variables_; }
  const std::vector<Variable::Name> &variableNames() const
  {
    auto variable_names = std::vector<Variable::Name>();
    variable_names.reserve(variables_.size());
    std::ranges::transform(variables_,
                           std::back_inserter(variable_names),
                           [](const Variable &var) { return var.name() });
    return variable_names;
  };
  const FactorTable &table() const { return table_; }
  void normalize()
  {
    auto z = std::ranges::accumulate(
        table_, 0.0, [](auto sum, const auto &p) { return sum + p.second; });
    std::ranges::for_each(table_, [](auto &p) { p.second /= z; });
  }

private:
  std::vector<Variable> variables_;
  FactorTable table_;
};

Assignment select(const Assignment &assignment,
                  const std::vector<Variable::Name> &variable_names)
{
  auto sub_assignment = Assignment();
  std::ranges::for_each(variable_names,
                        [&assignment](const Variable::Name &name) {
                          sub_assignment[name] = assignment[name];
                        });
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
  product.reserve(sub_product.empty() ? (*begin)
                                      : sub_product.size() * (*begin));

  auto base_fill = [&product](int i) {
    auto product_i = std::vector<int>({i});
    product.push_back(product_i);
  };
  auto normal_fill =
      [&product, &sub_product](int i) {
        for (const auto &sub_product_j : sub_product)
        {
          auto product_i = std::vector<int>({i});
          std::copy(sub_product, std::back_inserter(product_i));
          product.push_back(product_i);
        }
      }

  for (int i = 0; i < (*begin); ++i)
  {
    if (sub_product.empty())
      base_fill(i);
    else
      normal_fill(i);
  }

  return product;
}

std::vector<int> cartesianProduct(const std::vector<Variable> &variables)
{
  return cartesianProduct(variables.cbegin(), variables.cend());
}

std::vector<Assignment> assign(const std::vector<Variable> &variables)
{
  auto variable_names = std::vector<Variable::Name>();
  variable_names.reserve(variables_.size());
  std::ranges::transform(variables_,
                         std::back_inserter(variable_names),
                         [](const Variable &var) { return var.name() });

  auto values = cartesianProduct(variables);

  auto assignments = std::vector<Assignment>();
  assignments.reserve(values.size());
  for (const auto &value : values)
  {
    auto assignment = Assignment();
    for (auto i = 0; i < static_cast<int>(value.size()); ++i)
      assignment[variable_names[i]] = value[i];
    assignments.push_back(assignment);
  }

  return assignments;
}
} // namespace algo_dm