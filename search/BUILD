load("@rules_cc//cc:defs.bzl", "cc_binary", "cc_library", "cc_test")

package(default_visibility = ["//visibility:public"])

copts = [
    "--std=c++20",
    "-Wunused",
    "-Wextra",
]

# libraries
cc_library(
    name = "dijkstra_cc",
    hdrs = ["include/search/dijkstra.hpp"],
    copts = copts,
    strip_include_prefix = "include",
)

# binaries
cc_binary(
    name = "search_viz",
    srcs = ["test/search_viz.cpp"],
    copts = copts,
    deps = [
        ":dijkstra_cc",
    ],
)
