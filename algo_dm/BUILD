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
    copts = copts,
    strip_include_prefix = "include",
)

cc_library(
    name = "types_cc",
    srcs = ["src/types.cpp"],
    hdrs = ["include/algo_dm/types.hpp"],
    copts = copts,
    strip_include_prefix = "include",
    deps = [
        "@boost//:functional",
    ],
)

cc_library(
    name = "inference_cc",
    srcs = ["src/inference.cpp"],
    hdrs = ["include/algo_dm/inference.hpp"],
    copts = copts,
    strip_include_prefix = "include",
    deps = [
        ":types_cc",
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

cc_test(
    name = "inference_test",
    srcs = ["test/inference_test.cpp"],
    copts = copts,
    deps = [
        ":alias_cc",
        ":inference_cc",
        "@com_google_googletest//:gtest",
    ],
)
