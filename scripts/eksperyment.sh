#!/usr/bin/env bash
set -euo pipefail
cd -- "$(dirname -- "${BASH_SOURCE[0]}")/.."
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1
export MPLCONFIGDIR="$PWD/.local-tools/matplotlib" PYTHONPATH="$PWD"
exec "$PWD/.venv/bin/python" -Bu -m eksperyment.cli "$@"
