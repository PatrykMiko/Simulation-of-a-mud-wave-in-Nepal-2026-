#!/usr/bin/env bash
# Run from Bash on Linux/macOS or Git Bash on Windows.
set -euo pipefail
cd -- "$(dirname -- "${BASH_SOURCE[0]}")/.."
# All absolute paths passed to Multipass are guest paths. Host files use
# relative paths so that Git Bash and native Windows agree.
export MSYS2_ARG_CONV_EXCL='*'
machine=nepal30
guest=/home/ubuntu/nepal30
command=${1:-help}
if [[ "$command" == help ]]; then
    echo 'Użycie: eksperyment-vm.sh install CPU RAM DYSK | zasoby | start [1-5] | status | save-stop | recover | retry | export'
    exit 0
fi
command -v multipass >/dev/null || { echo 'Najpierw zainstaluj Canonical Multipass. Patrz SYMULACJA.md.' >&2; exit 1; }
vm_exec() { multipass exec "$machine" -- "$@"; }
case "$command" in
install)
    [[ $# == 4 ]] || { echo 'Przykład dla jednej próby: install 2 16G 80G' >&2; exit 1; }
    if ! multipass info "$machine" >/dev/null 2>&1; then
        multipass launch 24.04 --name "$machine" --cpus "$2" --memory "$3" --disk "$4"
    else
        multipass start "$machine"
    fi
    # Use only an empty project directory, or retry the exact same source
    # version after an interrupted transfer/install. Never overlay a study.
    revision=$(git rev-parse HEAD)
    vm_exec bash -eu -c '
        d=$1; rev=$2
        if test -e "$d/wyniki"; then echo "Istnieją wyniki. Użyj status lub start." >&2; exit 1; fi
        if test -e "$d" && ! test "$(cat "$d/.vm-revision" 2>/dev/null)" = "$rev"; then
            echo "Istnieje inna kopia projektu w maszynie nepal30." >&2; exit 1
        fi
        mkdir -p "$d"; echo "$rev" > "$d/.vm-revision"
    ' sh "$guest" "$revision"
    # Git preserves Linux line endings and executable bits on every host.
    archive=$(mktemp .nepal-wejscie.XXXXXXXX)
    trap 'rm -f -- "$archive"' EXIT
    git archive --format=tar HEAD > "$archive"
    multipass transfer "$archive" "$machine:nepal30-input.tar"
    vm_exec tar -xf /home/ubuntu/nepal30-input.tar -C "$guest"
    vm_exec rm -f /home/ubuntu/nepal30-input.tar
    vm_exec "$guest/scripts/eksperyment-install.sh"
    ;;
start|recover|retry)
    multipass start "$machine"
    if [[ "$command" == start && $# == 2 ]]; then
        [[ $2 =~ ^[1-5]$ ]] || { echo 'Wybierz od 1 do 5 prób.' >&2; exit 1; }
        vm_exec env "NEPAL_WORKERS=$2" "$guest/scripts/eksperyment.sh" start
    else
        vm_exec "$guest/scripts/eksperyment.sh" "$command"
    fi
    ;;
save-stop)
    vm_exec "$guest/scripts/eksperyment.sh" save-stop
    multipass stop "$machine"
    echo 'Zapis ukończony, maszyna wirtualna zatrzymana. Możesz wyłączyć komputer.'
    ;;
status|zasoby)
    vm_exec "$guest/scripts/eksperyment.sh" "$command"
    ;;
export)
    mkdir -p wyniki
    part=$(mktemp wyniki/paczka.XXXXXXXX)
    trap 'rm -f -- "$part"' EXIT
    vm_exec test -f "$guest/wyniki/finished.json"
    multipass transfer "$machine:$guest/wyniki/paczka.zip" "$part"
    if [[ -e wyniki/paczka.zip ]] && ! cmp -s "$part" wyniki/paczka.zip; then
        echo 'W katalogu gospodarza jest inne wyniki/paczka.zip. Przenieś je i ponów export.' >&2
        exit 1
    fi
    mv -- "$part" wyniki/paczka.zip
    echo 'Gotowe: wyniki/paczka.zip'
    ;;
*) echo "Nieznane polecenie: $command" >&2; exit 1;;
esac
