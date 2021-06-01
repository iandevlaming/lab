load("@rules_cc//cc:defs.bzl", "cc_binary", "cc_library", "cc_test")

package(default_visibility = ["//visibility:public"])

copts = [
    "--std=c++20",
    "-Wunused",
    "-Wextra",
]

# libraries
cc_library(
    name = "breadth_first_search_cc",
    hdrs = ["include/search/breadth_first_search.hpp"],
    copts = copts,
    strip_include_prefix = "include",
)

cc_library(
    name = "depth_first_search_cc",
    hdrs = ["include/search/depth_first_search.hpp"],
    copts = copts,
    strip_include_prefix = "include",
)

cc_library(
    name = "dijkstra_cc",
    hdrs = ["include/search/dijkstra.hpp"],
    copts = copts,
    strip_include_prefix = "include",
)

cc_library(
    name = "search_cc",
    hdrs = ["include/search/search.hpp"],
    copts = copts,
    strip_include_prefix = "include",
)

# binaries
cc_binary(
    name = "search_viz",
    srcs = ["test/search_viz.cpp"],
    copts = copts,
    deps = [
        ":breadth_first_search_cc",
        ":depth_first_search_cc",
        ":dijkstra_cc",
        ":search_cc",
    ],
)
