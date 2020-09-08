load("@rules_cc//cc:defs.bzl", "cc_binary", "cc_library", "cc_test")

package(default_visibility = ["//visibility:public"])

copts = [
    "--std=c++20",
    "-Wunused",
    "-Wextra",
]

# libraries
cc_library(
    name = "bracketing_cc",
    hdrs = ["include/algo_opt/bracketing.hpp"],
    copts = copts,
    strip_include_prefix = "include",
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
