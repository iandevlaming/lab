load("@rules_cc//cc:defs.bzl", "cc_library")

package(default_visibility = ["//visibility:public"])

cc_library(
    name = "gnuplot-iostream-lib",
    hdrs = ["gnuplot-iostream.h"],
    include_prefix = "gnuplot-iostream",
    deps = [
        "@boost//:filesystem",
        "@boost//:iostreams",
        "@boost//:tuple",
        "@boost//:utility",
        "@boost//:version",
    ],
)
