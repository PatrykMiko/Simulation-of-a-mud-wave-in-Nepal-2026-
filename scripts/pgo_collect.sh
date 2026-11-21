#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "$0")/.."
cmake -S . -B build-pgo -DSCIDDICA_PGO_GENERATE=ON -DSCIDDICA_PGO_USE=OFF
cmake --build build-pgo -j4
cp build-pgo/sciddica_ss3_cpp.cpython-315-x86_64-linux-gnu.so .
.venv/bin/python -m sciddica_ss3.cli data/nepal_ss3/alos --headless --steps 200 --use-cpp --use-float32 || true
cmake -S . -B build-pgo -DSCIDDICA_PGO_GENERATE=OFF -DSCIDDICA_PGO_USE=ON
cmake --build build-pgo -j4
cp build-pgo/sciddica_ss3_cpp.cpython-315-x86_64-linux-gnu.so .
