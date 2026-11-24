#!/usr/bin/env bash
set -euo pipefail
cd -- "$(dirname -- "${BASH_SOURCE[0]}")/.."
case "$PWD" in *[!a-zA-Z0-9/_.-]*) echo 'Sklonuj projekt w katalogu bez spacji i znaków specjalnych.' >&2; exit 1;; esac
if [[ -f wyniki/checkpoint/state.json ]]; then
    echo 'W katalogu są obliczenia. Wznowienie: ./scripts/eksperyment.sh start' >&2
    exit 1
fi
sudo apt-get update
sudo apt-get install -y build-essential curl unzip grass-core grass-dev gdal-bin libgdal-dev \
    python3-venv python3-pil ffmpeg zlib1g-dev util-linux
/usr/bin/python3 -m venv .venv
.venv/bin/python -m pip install -r requirements.txt
.venv/bin/python -B -m eksperyment.build
./scripts/install-dmtcp.sh
./scripts/eksperyment.sh verify
./scripts/eksperyment.sh zasoby
