
This repo contains our branch-and-bound solver to compute the twinwidth of graphs.
It was created for the [PACE challenge 2023](https://pacechallenge.org/2023/).

# Building

We use CMake and have no dependencies apart from C++20. 
To build the project, run the following commands:

```bash
mkdir build
cmake ..
make
```

# Usage

The executable will be `./build/main`. 
It takes a graph in the [PACE format](https://pacechallenge.org/2023/io/) and outputs a contraction sequence of minimal width to stdout.
Logging info is written to stderr.
