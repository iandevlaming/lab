workspace(name = "lab")

load("@bazel_tools//tools/build_defs/repo:http.bzl", "http_archive")
load("@bazel_tools//tools/build_defs/repo:git.bzl", "git_repository", "new_git_repository")

# boost
git_repository(
    name = "com_github_nelhage_rules_boost",
    commit = "1e3a69bf2d5cd10c34b74f066054cd335d033d71",
    remote = "https://github.com/nelhage/rules_boost",
    shallow_since = "1591047380 -0700",
)

load("@com_github_nelhage_rules_boost//:boost/boost.bzl", "boost_deps")

boost_deps()

# Google test repository
http_archive(
    name = "com_google_googletest",
    sha256 = "94c634d499558a76fa649edb13721dce6e98fb1e7018dfaeba3cd7a083945e91",
    strip_prefix = "googletest-release-1.10.0",
    url = "https://github.com/google/googletest/archive/release-1.10.0.zip",
)

# gnuplot-iostream repo
http_archive(
    name = "gnuplot-iostream",
    build_file = "@//third_party:gnuplot_iostream.BUILD",
    sha256 = "671df031cf24cb8623a8dd52e4c22ddedaf495633168debce86b39d1a3a99505",
    strip_prefix = "gnuplot-iostream-master",
    urls = ["https://github.com/dstahlke/gnuplot-iostream/archive/master.zip"],
)

# eigen
new_git_repository(
    name = "eigen",
    build_file = "@//third_party:eigen.BUILD",
    commit = "21ae2afd4edaa1b69782c67a54182d34efe43f9c",  # tag 3.3.7
    remote = "https://gitlab.com/libeigen/eigen",
    shallow_since = "1544551075 +0100",
)
