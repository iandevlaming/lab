#include <algo_dm/types.hpp>

#include <boost/functional/hash.hpp>

#include <numeric>
#include <sstream>
#include <stdexcept>

namespace algo_dm
{
Variable::Variable(const Name& name, int num_vals)
    : name_(name), num_vals_(num_vals)
{
  init();
}

Variable::Variable(Name&& name, int num_vals)
    : name_(std::move(name)), num_vals_(num_vals)
{
  init();
}

const Variable::Name& Variable::getName() const { return name_; }

int Variable::getNumVals() const { return num_vals_; }

void Variable::init()
{
  if (num_vals_ <= 0)
  {
    auto error_stream = std::ostringstream();
    error_stream << "Variable num_vals must be positive, but got " << num_vals_;
    throw std::domain_error(error_stream.str());
  }
}

bool operator==(const Variable& lhs, const Variable& rhs)
{
  return lhs.getName() == rhs.getName() && lhs.getNumVals() == rhs.getNumVals();
}

std::vector<Variable::Name> getNames(const std::vector<Variable>& variables)
{
  auto select_names = [](const Variable& v) { return v.getName(); };
  auto names_view = std::views::transform(variables, select_names);
  return std::vector<Variable::Name>(names_view.begin(), names_view.end());
}

std::ostream& operator<<(std::ostream& os, const Variable& var)
{
  os << "Variable " << var.getName() << " (" << var.getNumVals() << ")";
  return os;
}

Assignment::Assignment(const std::map<Key, Value>& assignment)
    : assignment_(assignment)
{
}

void Assignment::set(const Key& key, Value value) { assignment_[key] = value; }

const Assignment::Value& Assignment::get(const Key& key) const
{
  return assignment_.at(key);
}

bool Assignment::contains(const Key& key) const
{
  return assignment_.contains(key);
}

void Assignment::erase(const Key& key)
{
  if (assignment_.contains(key))
    assignment_.erase(key);
}

std::vector<Assignment::Key> Assignment::getVariableNames() const
{
  return algo_dm::getKeys(assignment_);
}

bool Assignment::operator==(const Assignment& other) const
{
  return assignment_ == other.assignment_;
}

typename Assignment::ItrT Assignment::begin() { return assignment_.begin(); }

typename Assignment::ItrT Assignment::end() { return assignment_.end(); }

typename Assignment::CItrT Assignment::begin() const { return cbegin(); }

typename Assignment::CItrT Assignment::end() const { return cend(); }

typename Assignment::CItrT Assignment::cbegin() const
{
  return assignment_.cbegin();
}

typename Assignment::CItrT Assignment::cend() const
{
  return assignment_.cend();
}

size_t Assignment::Hash::operator()(const Assignment& assignment) const
{
  auto seed = static_cast<size_t>(0);
  auto pairwise_hash = [&seed](const auto& p) {
    boost::hash_combine(seed, boost::hash_value(p.first));
    boost::hash_combine(seed, boost::hash_value(p.second));
  };

  std::ranges::for_each(assignment.assignment_, pairwise_hash);
  return seed;
}

std::ostream& operator<<(std::ostream& os, const Assignment& a)
{
  const auto& vars = a.getVariableNames();
  for (int i = 0; i < static_cast<int>(vars.size()); ++i)
  {
    os << vars[i] << ": " << a.get(vars[i]);
    if (i != static_cast<int>(vars.size()) - 1)
      os << ", ";
  }
  return os;
}

FactorTable::FactorTable(
    const std::map<Assignment, double, AssignmentCompare>& table)
{
  checkInput(table);
  table_ = table;
}

FactorTable::FactorTable(
    const std::unordered_map<Assignment, double, Assignment::Hash>& table)
{
  checkInput(table);
  table_ = std::map<Assignment, double, AssignmentCompare>(table.cbegin(),
                                                           table.cend());
}

FactorTable::FactorTable(const std::vector<Assignment>& assignments,
                         const std::vector<double>& probabilities)
{
  if (assignments.size() != probabilities.size())
    throw std::invalid_argument(
        "Each assignment in a FactorTable must have exactly one probability");

  if (!assignments.empty())
    for (auto i = 0; i < static_cast<int>(assignments.size()); ++i)
      set(assignments[i], probabilities[i]);
}

const FactorTable::Value& FactorTable::get(const Key& key) const
{
  return table_.at(key);
}

bool FactorTable::contains(const Key& key) const
{
  return table_.contains(key);
}

void FactorTable::set(const Key& key, Value value)
{
  if (!table_.empty() && variable_names_ != key.getVariableNames())
    throw std::invalid_argument(
        "All Assignments in a FactorTable must have the same keys");

  if (table_.empty())
    variable_names_ = key.getVariableNames();

  table_[key] = value;
}

std::vector<FactorTable::Key> FactorTable::getAssignments() const
{
  return algo_dm::getKeys(table_);
}

std::vector<Assignment::Key> FactorTable::getVariableNames() const
{
  return variable_names_;
}

bool operator==(const FactorTable& lhs, const FactorTable& rhs)
{
  const auto& assignments = lhs.getAssignments();
  if (assignments != rhs.getAssignments())
    return false;

  for (const auto& a : assignments)
    if (std::abs(lhs.get(a) - rhs.get(a)) > 1e-12)
      return false;

  return true;
}

void FactorTable::normalize()
{
  auto sum_value = [](auto sum, const auto& p) { return sum + p.second; };
  auto z = std::accumulate(table_.cbegin(), table_.cend(), 0.0, sum_value);
  std::ranges::for_each(table_, [&z](auto& p) { p.second /= z; });
}

bool FactorTable::AssignmentCompare::operator()(const Assignment& lhs,
                                                const Assignment& rhs) const
{
  const auto& vars = lhs.getVariableNames();

  for (const auto& var : vars)
  {
    const auto lhs_val = lhs.get(var);
    const auto rhs_val = rhs.get(var);

    if (lhs_val < rhs_val)
      return true;
    else if (lhs_val > rhs_val)
      return false;
  }

  return false;
}

std::ostream& operator<<(std::ostream& os, const FactorTable& table)
{
  const auto& vars = table.getVariableNames();
  const auto& assignments = table.getAssignments();

  for (int row = -1; row < static_cast<int>(assignments.size()); ++row)
  {
    for (int col = 0; col < static_cast<int>(vars.size()); ++col)
    {
      if (row < 0)
        os << vars[col];
      else
        os << assignments[row].get(vars[col]);
      os << "\t";

      if (col == static_cast<int>(vars.size()) - 1)
      {
        if (row < 0)
          os << "P";
        else
          os << table.get(assignments[row]);
      }
    }
    if (row != static_cast<int>(assignments.size()) - 1)
      os << "\n";
  }

  return os;
}

Factor::Factor(const std::vector<Variable>& variables, const FactorTable& table)
    : variables_(variables), table_(table)
{
}

Factor::Factor(std::vector<Variable>&& variables, FactorTable&& table)
    : variables_(std::move(variables)), table_(std::move(table))
{
}

const std::vector<Variable>& Factor::getVariables() const { return variables_; }

const FactorTable& Factor::getFactorTable() const { return table_; }

std::ostream& operator<<(std::ostream& os, const Factor& f)
{
  const auto& vars = f.getVariables();
  const auto& table = f.getFactorTable();

  os << "Variables: ";
  for (int i = 0; i < static_cast<int>(vars.size()); ++i)
  {
    os << vars[i];
    if (i != static_cast<int>(vars.size()) - 1)
      os << ", ";
  }

  os << "\n\n" << table;

  return os;
}

void Factor::normalize() { table_.normalize(); }

bool operator==(const Factor& lhs, const Factor& rhs)
{
  auto equal = lhs.getFactorTable() == rhs.getFactorTable();
  equal = equal && isPermutation(lhs.getVariables(), rhs.getVariables());
  return equal;
}

bool isInScope(const Factor& factor, const Variable::Name name)
{
  auto is_named = [](const auto& var_name) {
    return [&var_name](const auto& var) { return var.getName() == var_name; };
  };

  return std::ranges::any_of(factor.getVariables(), is_named(name));
}

Factor operator*(const Factor& lhs, const Factor& rhs)
{
  const auto& lhs_vars = lhs.getVariables();
  const auto& rhs_vars = rhs.getVariables();

  const auto& lhs_table = lhs.getFactorTable();
  const auto& rhs_table = rhs.getFactorTable();

  const auto rhs_names = getNames(rhs_vars);

  const auto rhs_only_vars = setDiff(rhs_vars, lhs_vars);
  const auto rhs_only_assignments = rhs_only_vars.empty()
                                        ? std::vector({Assignment()})
                                        : assign(rhs_only_vars);

  auto merged_table = FactorTable();
  for (const auto& lhs_assignment : lhs_table.getAssignments())
  {
    for (const auto& rhs_only_assignment : rhs_only_assignments)
    {
      auto merged_assignment = rhs_only_assignment;
      for (const auto& lhs_var_name : lhs_assignment.getVariableNames())
        merged_assignment.set(lhs_var_name, lhs_assignment.get(lhs_var_name));

      const auto rhs_assignment = select(merged_assignment, rhs_names);
      const auto merged_probability =
          lhs_table.get(lhs_assignment) * rhs_table.get(rhs_assignment);

      merged_table.set(merged_assignment, merged_probability);
    }
  }

  auto merged_vars = lhs_vars;
  std::ranges::copy(rhs_only_vars, std::back_inserter(merged_vars));

  return Factor(merged_vars, merged_table);
}

std::vector<Variable::Name> getNames(const std::vector<Factor>& factors)
{
  auto get_names = [](const auto& names, const auto& factor) {
    return setUnion(names, getNames(factor.getVariables()));
  };
  return std::accumulate(factors.cbegin(),
                         factors.cend(),
                         std::vector<Variable::Name>(),
                         get_names);
}

BayesianNetwork::BayesianNetwork(
    const std::unordered_map<Variable::Name, Factor>& nodes,
    const AdjacencyList<Variable::Name>& graph)
    : nodes_(nodes), graph_(graph)
{
}

std::vector<Variable> BayesianNetwork::getVariables() const
{
  auto get_variables = [](const auto& f) { return f.getVariables(); };
  auto variables_view = std::ranges::views::common(
      nodes_ | std::views::values | std::views::transform(get_variables));
  auto join = [](const auto& lhs, const auto& rhs) {
    return setUnion(lhs, rhs);
  };
  return std::accumulate(variables_view.begin(),
                         variables_view.end(),
                         std::vector<Variable>(),
                         join);
}

std::vector<Factor> BayesianNetwork::getFactors() const
{
  auto values_view = std::ranges::views::common(nodes_ | std::views::values);
  return std::vector(values_view.begin(), values_view.end());
}

const std::unordered_map<Variable::Name, Factor>&
BayesianNetwork::getNodes() const
{
  return nodes_;
}

const AdjacencyList<Variable::Name>& BayesianNetwork::getGraph() const
{
  return graph_;
}

const Factor& BayesianNetwork::getFactor(const Variable::Name& node) const
{
  return nodes_.at(node);
}

Assignment select(const Assignment& assignment,
                  const std::vector<Variable::Name>& variable_names)
{
  auto sub_assignment = Assignment();
  auto sub_assign = [&assignment, &sub_assignment](const Variable::Name& name) {
    if (!assignment.contains(name))
      throw std::invalid_argument(
          "Cannot select variable names not contained by Assignment");
    sub_assignment.set(name, assignment.get(name));
  };
  std::ranges::for_each(variable_names, sub_assign);

  return sub_assignment;
}

std::vector<Variable> select(const std::vector<Variable>& variables,
                             const std::vector<Variable::Name>& variable_names)
{
  auto is_selected = [& n = variable_names](const Variable& v) {
    return contains(n, v.getName());
  };
  auto vars_view = variables | std::views::filter(is_selected);
  return std::vector(vars_view.begin(), vars_view.end());
}

std::vector<Variable> erase(const std::vector<Variable>& variables,
                            const Variable::Name& name)
{
  auto is_not_named = [](const auto& var_name) {
    return [&var_name](const auto& var) { return var.getName() != var_name; };
  };

  auto reduced_view = std::views::filter(variables, is_not_named(name));
  return std::vector<Variable>(reduced_view.begin(), reduced_view.end());
}

// TODO: handle 0 num values
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
      assignment.set(names[i], value[i]);
    assignments.push_back(assignment);
  }

  return assignments;
}

double computeProbability(const std::vector<FactorTable>& tables,
                          const Assignment& assignment)
{
  auto get_probability = [&assignment](const FactorTable& table) {
    auto vars = table.getVariableNames();
    auto sub_assignment = select(assignment, vars);
    return table.contains(sub_assignment) ? table.get(sub_assignment) : 0.0;
  };

  auto accumulate_probability = [&get_probability](double p,
                                                   const FactorTable& t) {
    return p * get_probability(t);
  };

  return std::accumulate(
      tables.cbegin(), tables.cend(), 1.0, accumulate_probability);
}
} // namespace algo_dm