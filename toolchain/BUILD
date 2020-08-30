package(default_visibility = ["//visibility:public"])

load(":cc_toolchain_config.bzl", "cc_toolchain_config")

filegroup(
    name = "empty",
)

cc_toolchain_config(
    name = "gcc_k8_toolchain_config",
)

cc_toolchain(
    name = "gcc_k8_toolchain",
    all_files = ":empty",
    compiler_files = ":empty",
    dwp_files = ":empty",
    linker_files = "empty",
    objcopy_files = ":empty",
    strip_files = ":empty",
    supports_param_files = 0,
    toolchain_config = ":gcc_k8_toolchain_config",
    toolchain_identifier = "gcc-k8-toolchain",
)

cc_toolchain_suite(
    name = "gcc_suite",
    toolchains = {
        "k8|gcc": ":gcc_k8_toolchain",
    },
)
