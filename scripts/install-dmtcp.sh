#!/usr/bin/env bash
set -euo pipefail
cd -- "$(dirname -- "${BASH_SOURCE[0]}")/.."
prefix="$PWD/.local-tools/dmtcp-4.2.0"
project_root="$PWD"
build_dir="$(mktemp -d /tmp/nepal-dmtcp-install.XXXXXXXX)"
trap 'rm -rf -- "$build_dir"' EXIT
curl -fL --retry 2 https://github.com/dmtcp/dmtcp/archive/refs/tags/v4.2.0.tar.gz -o "$build_dir/source.tar.gz"
echo '043410566fd7c09f21e0ec485cf72e9481722629ea877650142b7505e0d7f2d2  '"$build_dir/source.tar.gz" | sha256sum -c -
tar -xzf "$build_dir/source.tar.gz" -C "$build_dir"
cd "$build_dir/dmtcp-4.2.0"
./configure --prefix="$prefix"
make -j2
make install
"$project_root/scripts/build-checkpoint-plugin.sh"
echo "Zainstalowano lokalnie: $prefix. Nie uruchomiono modelu."
