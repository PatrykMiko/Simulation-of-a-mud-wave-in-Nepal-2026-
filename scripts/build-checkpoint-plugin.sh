#!/usr/bin/env bash
set -euo pipefail
cd -- "$(dirname -- "${BASH_SOURCE[0]}")/.."
test -f .local-tools/dmtcp-4.2.0/include/dmtcp.h || {
    echo 'Najpierw wykonaj ./scripts/install-dmtcp.sh' >&2
    exit 1
}
cc -shared -fPIC -Wall -Wextra -Werror -O2 \
    '-DDMTCP_PACKAGE_VERSION="4.2.0"' \
    -I .local-tools/dmtcp-4.2.0/include scripts/dmtcp-file-policy.c \
    -o .local-tools/libnepal_checkpoint_files.so.tmp
mv -- .local-tools/libnepal_checkpoint_files.so.tmp .local-tools/libnepal_checkpoint_files.so
echo 'Zbudowano ochronę plików środowiska dla checkpointów.'
