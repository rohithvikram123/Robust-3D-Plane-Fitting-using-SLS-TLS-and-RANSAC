#!/usr/bin/env bash
set -e
cd "$(dirname "$0")"
g++ -std=c++17 -I /usr/include/eigen3 main.cpp surfaceFitting/surfaceFitting.cpp -I surfaceFitting -o main `pkg-config --cflags --libs opencv4` $(python3-config --includes) $(python3-config --embed --ldflags)
env -u LD_LIBRARY_PATH ./main
