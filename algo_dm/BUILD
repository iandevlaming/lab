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
    hdrs = ["include/algo_dm/alias.hpp"],
    strip_include_prefix = "include",
    copts = copts,
)

cc_library(
    name = "types_cc",
    hdrs = ["include/algo_dm/types.hpp"],
    strip_include_prefix = "include",
    copts = copts,
    deps = [
        "@boost//:functional",
    ],
)

# binaries

# tests
cc_test(
    name = "types_test",
    srcs = ["test/types_test.cpp"],
    copts = copts,
    deps = [
        ":alias_cc",
        ":types_cc",
        "@com_google_googletest//:gtest",
    ],
)