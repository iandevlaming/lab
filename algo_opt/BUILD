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
