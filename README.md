# Lab
This is a testbed for algorithms, libraries, etc. It's not pretty, and despite the pretty visualizations, it's probably not correct
* There are no (few) tests
  * This may change if I get sufficiently annoyed at a bug, can't visualize code, or if the wind blows in a particular direction
* There are no (few) comments
  * This will almost certainly change as soon as I forget what my old code did, which is today and every day
  * Some of my comments look like doxy if you keep your computer screen several feet away
* Commits into master are made as pull requests even though I'm the only one using this repo
  * I've decided to die on this hill
* There is no CI
  * No tests to run... _yet_
* It's mostly C++
  * C++ is what I use for development, and doing things the hard way is a great learning experience. Someday I'll give up and start using python and thoroughly regret all the time I spent writing C++.

I've tried and mostly failed to make code generic. Some of the generic code is _too_ generic. Look, I know the code isn't all that pretty, If you want to use it (or you'd just like to know what's going on) send me an angry email.

## Algo Opt
This library is _mostly_ algorithms for solving optimization problems from the book _Algorithms For Optimization_. You'll also find supporting code. The original code was written in Julia (which implemented all the previously mentioned supporting code for free), and it is painstakingly translated into C++ here because I enjoy pain. Visualizations are produced using a gnuplot wrapper. See previous sentence for motivation behind that horrible decision.

# Running Code
This code is built using Bazel and compiled with gcc 10 for C++20 language support. There are other system dependecies (e.g. gnuplot), but instead of installing all those yourself, why not head on over to other repo, `rocket`, for a docker file with everything set up for free?
