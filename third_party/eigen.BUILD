load("@rules_cc//cc:defs.bzl", "cc_library")

package(default_visibility = ["//visibility:public"])

cc_library(
    name = "eigen-lib",
    hdrs = glob(
        [
            "Eigen/**",
        ],
        exclude = ["Eigen/**/CMakeLists.txt"],
    ),
    defines = [
        "EIGEN_MAX_CPP_VER=20",
    ],
    includes = ["."],
)
