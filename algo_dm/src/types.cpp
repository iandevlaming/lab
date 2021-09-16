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

std::optional<Assignment::Value> Assignment::get(const Key& key) const
{
  if (assignment_.contains(key))
    return assignment_.at(key);
  return {};
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
    os << vars[i] << ": " << a.get(vars[i]).value();
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

std::optional<FactorTable::Value> FactorTable::get(const Key& key) const
{
  if (table_.contains(key))
    return table_.at(key);
  return {};
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
    if (std::abs(lhs.get(a).value() - rhs.get(a).value()) > 1e-12)
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
    const auto lhs_val = lhs.get(var).value();
    const auto rhs_val = rhs.get(var).value();

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
        os << assignments[row].get(vars[col]).value();
      os << "\t";

      if (col == static_cast<int>(vars.size()) - 1)
      {
        if (row < 0)
          os << "P";
        else
          os << table.get(assignments[row]).value();
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

AdjacencyList::AdjacencyList(int num_nodes)
    : num_nodes_(num_nodes),
      graph_(std::vector<std::unordered_set<int>>(num_nodes,
                                                  std::unordered_set<int>()))
{
}

int AdjacencyList::getNumNodes() const { return num_nodes_; }

void AdjacencyList::addEdge(int from_node, int to_node)
{
  if (from_node >= num_nodes_ || to_node >= num_nodes_)
    throw std::invalid_argument("Edges must be between nodes on the graph");

  graph_[from_node].insert(to_node);
}
const std::unordered_set<int>& AdjacencyList::getEdges(int node) const
{
  if (node >= num_nodes_)
    throw std::invalid_argument("Cannot get edges for nodes off of the graph");

  return graph_[node];
}

BayesianNetwork::BayesianNetwork(const std::vector<Variable>& variables,
                                 const std::vector<Factor>& factors,
                                 const AdjacencyList& graph)
    : variables_(variables), factors_(factors), graph_(graph)
{
  if (variables_.size() != factors.size())
    throw std::invalid_argument("Every node must have a factor");
}

const std::vector<Variable>& BayesianNetwork::getVariables() const
{
  return variables_;
}

const std::vector<Factor>& BayesianNetwork::getFactors() const
{
  return factors_;
}

const AdjacencyList& BayesianNetwork::getGraph() const { return graph_; }

Assignment select(const Assignment& assignment,
                  const std::vector<Variable::Name>& variable_names)
{
  auto sub_assignment = Assignment();
  auto sub_assign = [&assignment, &sub_assignment](const Variable::Name& name) {
    auto value = assignment.get(name);
    if (!value.has_value())
      throw std::invalid_argument(
          "Cannot select variable names not contained by Assignment");
    sub_assignment.set(name, assignment.get(name).value());
  };
  std::ranges::for_each(variable_names, sub_assign);

  return sub_assignment;
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
    auto p = table.get(sub_assignment);
    return p.has_value() ? p.value() : 0.0;
  };

  auto accumulate_probability = [&get_probability](double p,
                                                   const FactorTable& t) {
    return p * get_probability(t);
  };

  return std::accumulate(
      tables.cbegin(), tables.cend(), 1.0, accumulate_probability);
}
} // namespace algo_dm