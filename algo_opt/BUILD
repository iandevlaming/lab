load("@rules_cc//cc:defs.bzl", "cc_binary", "cc_library", "cc_test")

package(default_visibility = ["//visibility:public"])

copts = [
    "--std=c++20",
    "-Wunused",
    "-Wextra",
]

# libraries
cc_library(
    name = "alias_cc",
    hdrs = ["include/algo_opt/alias.hpp"],
    copts = copts,
    strip_include_prefix = "include",
    deps = [
        "@eigen//:eigen-lib",
    ],
)

cc_library(
    name = "bracketing_cc",
    hdrs = ["include/algo_opt/bracketing.hpp"],
    copts = copts,
    strip_include_prefix = "include",
    deps = [
        ":alias_cc",
    ],
)

cc_library(
    name = "local_descent_cc",
    hdrs = ["include/algo_opt/local_descent.hpp"],
    copts = copts,
    strip_include_prefix = "include",
    deps = [
        ":alias_cc",
        ":bracketing_cc",
        "@eigen//:eigen-lib",
    ],
)

cc_library(
    name = "first_order_cc",
    hdrs = ["include/algo_opt/first_order.hpp"],
    copts = copts,
    strip_include_prefix = "include",
    deps = [
        ":alias_cc",
        ":local_descent_cc",
    ],
)

cc_library(
    name = "second_order_cc",
    hdrs = ["include/algo_opt/second_order.hpp"],
    copts = copts,
    strip_include_prefix = "include",
    deps = [
        ":alias_cc",
        "@eigen//:eigen-lib",
    ],
)

cc_library(
    name = "direct_methods_cc",
    hdrs = ["include/algo_opt/direct_methods.hpp"],
    copts = copts,
    strip_include_prefix = "include",
    deps = [
        ":alias_cc",
        ":local_descent_cc",
        "@eigen//:eigen-lib",
    ],
)

cc_library(
    name = "stochastic_methods_cc",
    hdrs = ["include/algo_opt/stochastic_methods.hpp"],
    copts = copts,
    strip_include_prefix = "include",
    deps = [
        ":alias_cc",
        "@eigen//:eigen-lib",
    ],
)

cc_library(
    name = "lco_cc",
    hdrs = ["include/algo_opt/lco.hpp"],
    copts = copts,
    strip_include_prefix = "include",
    deps = [
        ":alias_cc",
        "@eigen//:eigen-lib",
    ],
)

# tests
cc_test(
    name = "bracketing_test",
    srcs = ["test/bracketing_test.cpp"],
    copts = copts,
    deps = [
        ":bracketing_cc",
        "@com_google_googletest//:gtest",
    ],
)

cc_test(
    name = "lco_test",
    srcs = ["test/lco_test.cpp"],
    copts = copts,
    deps = [
        ":lco_cc",
        "@com_google_googletest//:gtest",
    ],
)

# viz tests
cc_binary(
    name = "bracketing_viz",
    srcs = ["test/bracketing_viz.cpp"],
    copts = copts,
    deps = [
        ":bracketing_cc",
        "@boost//:tuple",
        "@gnuplot-iostream//:gnuplot-iostream-lib",
    ],
)

cc_binary(
    name = "local_descent_viz",
    srcs = ["test/local_descent_viz.cpp"],
    copts = copts,
    deps = [
        ":local_descent_cc",
        "@eigen//:eigen-lib",
        "@gnuplot-iostream//:gnuplot-iostream-lib",
    ],
)

cc_binary(
    name = "first_order_viz",
    srcs = ["test/first_order_viz.cpp"],
    copts = copts,
    deps = [
        ":alias_cc",
        ":first_order_cc",
        "@gnuplot-iostream//:gnuplot-iostream-lib",
    ],
)

cc_binary(
    name = "second_order_viz",
    srcs = ["test/second_order_viz.cpp"],
    copts = copts,
    deps = [
        ":alias_cc",
        ":local_descent_cc",
        ":second_order_cc",
        "@gnuplot-iostream//:gnuplot-iostream-lib",
    ],
)

cc_binary(
    name = "direct_methods_viz",
    srcs = ["test/direct_methods_viz.cpp"],
    copts = copts,
    deps = [
        ":alias_cc",
        ":direct_methods_cc",
        ":local_descent_cc",
        "@gnuplot-iostream//:gnuplot-iostream-lib",
    ],
)

cc_binary(
    name = "stochastic_methods_viz",
    srcs = ["test/stochastic_methods_viz.cpp"],
    copts = copts,
    deps = [
        ":alias_cc",
        ":local_descent_cc",
        ":stochastic_methods_cc",
        "@gnuplot-iostream//:gnuplot-iostream-lib",
    ],
)
