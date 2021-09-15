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
  auto variable_names = std::vector<Variable::Name>();
  variable_names.reserve(variables.size());

  auto select_names = [](const Variable& v) { return v.getName(); };
  std::ranges::transform(
      variables, std::back_inserter(variable_names), select_names);

  return variable_names;
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

std::vector<Assignment::Key> Assignment::getKeys() const
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

FactorTable::FactorTable(
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
  if (!table_.empty() && variable_names_ != key.getKeys())
    throw std::invalid_argument(
        "All Assignments in a FactorTable must have the same keys");

  if (table_.empty())
    variable_names_ = key.getKeys();

  table_[key] = value;
}

std::vector<FactorTable::Key> FactorTable::getKeys() const
{
  return algo_dm::getKeys(table_);
}

std::vector<Assignment::Key> FactorTable::getAssignmentKeys() const
{
  return variable_names_;
}

bool FactorTable::operator==(const FactorTable& rhs) const
{
  return table_ == rhs.table_;
}

void FactorTable::normalize()
{
  auto sum_value = [](auto sum, const auto& p) { return sum + p.second; };
  auto z = std::accumulate(table_.cbegin(), table_.cend(), 0.0, sum_value);
  std::ranges::for_each(table_, [&z](auto& p) { p.second /= z; });
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

void Factor::normalize() { table_.normalize(); }

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
    auto vars = table.getAssignmentKeys();
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