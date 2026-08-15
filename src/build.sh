#!/usr/bin/env bash
set -e
echo "[1/3] deps"; sudo apt-get install -y -q libprimesieve-dev g++ cmake git
echo "[2/3] sieve128"
g++ -O3 -march=native -fopenmp -std=c++17 sieve128.cpp -o sieve128 -lprimesieve
echo "[3/3] primecount (128-bit)"
[ -d primecount ] || git clone -q --depth 1 https://github.com/kimwalisch/primecount
cmake -S primecount -B primecount/build -DCMAKE_BUILD_TYPE=Release > /dev/null
cmake --build primecount/build -j"$(nproc)" > /dev/null
echo "done. quick check:"
./sieve128 count 100000000000000000000 100000000000000010000
./primecount/build/primecount 1e13
echo "now run:  python3 ramanujan128.py validate --upto 15 --D 1e8 --primecount ./primecount/build/primecount"
