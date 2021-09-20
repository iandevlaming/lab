#include <algo_dm/alias.hpp>
#include <algo_dm/inference.hpp>

#include <gtest/gtest.h>

#include <vector>

TEST(ProductTest, Example3P1)
{
  const auto x = ad::Variable("x", 2);
  const auto y = ad::Variable("y", 2);
  const auto z = ad::Variable("z", 2);

  const auto f1_table =
      ad::FactorTable(ad::assign({x, y}), std::vector({0.3, 0.4, 0.2, 0.1}));
  const auto f1 = ad::Factor({x, y}, f1_table);

  const auto f2_table =
      ad::FactorTable(ad::assign({y, z}), std::vector({0.2, 0.0, 0.3, 0.5}));
  const auto f2 = ad::Factor({y, z}, f2_table);

  auto expected_table = ad::FactorTable(
      ad::assign({x, y, z}),
      std::vector({0.06, 0.00, 0.12, 0.20, 0.04, 0.00, 0.03, 0.05}));
  auto expected_product = ad::Factor({x, y, z}, expected_table);

  auto product = f1 * f2;

  EXPECT_EQ(product, expected_product);
  EXPECT_EQ(product, f2 * f1);
}

TEST(MarginalizeTest, Example3P2)
{
  const auto x = ad::Variable("x", 2);
  const auto y = ad::Variable("y", 2);
  const auto z = ad::Variable("z", 2);

  const auto table = ad::FactorTable(
      ad::assign({x, y, z}),
      std::vector({0.08, 0.31, 0.09, 0.37, 0.01, 0.05, 0.02, 0.07}));
  const auto factor = ad::Factor({x, y, z}, table);

  const auto marginalized_factor = ad::marginalize(factor, "y");

  const auto expected_table = ad::FactorTable(
      ad::assign({x, z}), std::vector({0.17, 0.68, 0.03, 0.12}));
  const auto expected_factor = ad::Factor({x, z}, expected_table);

  EXPECT_EQ(marginalized_factor, expected_factor);
}

TEST(ConditionTest, Example3P3)
{
  const auto x = ad::Variable("x", 2);
  const auto y = ad::Variable("y", 2);
  const auto z = ad::Variable("z", 2);

  const auto table = ad::FactorTable(
      ad::assign({x, y, z}),
      std::vector({0.08, 0.31, 0.09, 0.37, 0.01, 0.05, 0.02, 0.07}));
  const auto factor = ad::Factor({x, y, z}, table);

  const auto conditioned_factor = ad::condition(factor, "y", 1);

  const auto expected_table = ad::FactorTable(
      ad::assign({x, z}), std::vector({0.09, 0.37, 0.02, 0.07}));
  const auto expected_factor = ad::Factor({x, z}, expected_table);

  EXPECT_EQ(conditioned_factor, expected_factor);
}

class Exercise3P4 : public ::testing::Test
{
protected:
  Exercise3P4()
  {
    const auto c = ad::Variable(
        "c", 3); // c (class) in {vehicle (0), pedestrian (1), ball (2)}
    const auto s =
        ad::Variable("s", 3); // s (size) in {small (0), medium (1), large (2)}
    const auto v = ad::Variable(
        "v", 3); // v (velocity) in {slow (0), moderate (1), fast(2)}
    variables_ = {c, s, v};

    const auto p_c =
        ad::FactorTable(ad::assign({c}), std::vector{0.80, 0.19, 0.01});
    const auto node_c = ad::Factor({c}, p_c);

    const auto p_s_given_c = ad::FactorTable(
        ad::assign({c, s}),
        std::vector{
            0.001, 0.009, 0.990, 0.200, 0.75, 0.050, 0.800, 0.199, 0.001});
    const auto node_s = ad::Factor({c, s}, p_s_given_c);

    const auto p_v_given_c = ad::FactorTable(
        ad::assign({c, v}),
        std::vector({0.2, 0.2, 0.6, 0.5, 0.4, 0.1, 0.4, 0.4, 0.2}));
    const auto node_v = ad::Factor({c, v}, p_v_given_c);

    auto nodes = std::unordered_map<ad::Variable::Name, ad::Factor>();
    nodes[c.getName()] = node_c;
    nodes[s.getName()] = node_s;
    nodes[v.getName()] = node_v;

    auto graph = ad::AdjacencyList<ad::Variable::Name>();
    graph.addEdge(c.getName(), s.getName());
    graph.addEdge(c.getName(), v.getName());

    bayesian_network_ = ad::BayesianNetwork(nodes, graph);

    query_.push_back(c.getName());

    evidence_[s.getName()] = 1;
    evidence_[v.getName()] = 0;

    const auto p_0 = 0.80 * 0.009 * 0.2;
    const auto p_1 = 0.19 * 0.750 * 0.5;
    const auto p_2 = 0.01 * 0.199 * 0.4;

    expected_table_ = ad::FactorTable(ad::assign({c}), {p_0, p_1, p_2});
    expected_table_.normalize();
  }

  std::vector<ad::Variable> variables_;
  ad::BayesianNetwork bayesian_network_;
  std::vector<ad::Variable::Name> query_;
  std::unordered_map<ad::Variable::Name, ad::Assignment::Value> evidence_;
  ad::FactorTable expected_table_;
};

TEST_F(Exercise3P4, ExactInferenceTest)
{
  const auto method = ad::ExactInference();
  auto factor_c = ad::infer(method, bayesian_network_, query_, evidence_);

  EXPECT_EQ(factor_c.getFactorTable(), expected_table_);
}

TEST_F(Exercise3P4, SumProductVariableEliminatioInferenceTest)
{
  auto get_name = [](const auto& v) { return v.getName(); };
  auto ordering_view =
      variables_ | std::views::transform(get_name) | std::views::reverse;
  const auto ordering = std::vector(ordering_view.begin(), ordering_view.end());

  const auto method = ad::VariableElimination({ordering});
  const auto factor_c = ad::infer(method, bayesian_network_, query_, evidence_);

  EXPECT_EQ(factor_c.getFactorTable(), expected_table_);
}

TEST_F(Exercise3P4, DirectSampling)
{
  const auto method = ad::DirectSampling({10000});
  const auto factor_c = ad::infer(method, bayesian_network_, query_, evidence_);

  const auto& table = factor_c.getFactorTable();
  const auto& assignments = table.getAssignments();

  auto get_prob_of = [](const auto& t) {
    return [&t](const auto& a) { return t.get(a).value(); };
  };
  auto expected_prob_view =
      assignments | std::views::transform(get_prob_of(expected_table_));
  auto expected_prob =
      std::vector(expected_prob_view.begin(), expected_prob_view.end());

  auto prob_view = assignments | std::views::transform(get_prob_of(table));
  auto prob = std::vector(prob_view.begin(), prob_view.end());

  auto sorted_idx_view =
      std::views::iota(0, static_cast<int>(assignments.size()));
  auto sorted_idx = std::vector(sorted_idx_view.begin(), sorted_idx_view.end());
  auto sort_indices_using = [](const auto& c) {
    return [&c](auto i, auto j) { return c[i] < c[j]; };
  };
  std::ranges::sort(sorted_idx, sort_indices_using(expected_prob));

  auto access = [](const auto& c) { return [&c](auto i) { return c[i]; }; };
  auto sorted_prob = sorted_idx | std::views::transform(access(prob));

  EXPECT_TRUE(std::ranges::is_sorted(sorted_prob));
}

int main(int argc, char** argv)
{
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}