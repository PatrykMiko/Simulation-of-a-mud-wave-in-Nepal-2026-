"""Uruchamianie r.avaflow 4.0G z pliku konfiguracji YAML.

Z konfiguracji generuje skrypt startowy (taki jak ca.avaflow.start.sh z paczki
Cascadia), uruchamia go przez `grass <mapset> --exec sh <skrypt>`, mierzy czas
i wczytuje plik *_evaluation.txt.

Użycie (z katalogu głównego projektu, poza sesją GRASS i najlepiej bez condy):
    python -m ravaflow_tools.runner configs/ca_exp2.yaml
"""
from __future__ import annotations

import subprocess
import os
import re
import shlex
import signal
import sys
import time
from pathlib import Path

import yaml


def load_config(path: str | Path) -> dict:
    with open(path, encoding="utf-8") as f:
        cfg = yaml.safe_load(f)
    for key in ("workdir", "grass_mapset", "inputs", "params"):
        if key not in cfg:
            raise ValueError(f"Brak klucza '{key}' w konfiguracji {path}")
    if "prefix" not in cfg["params"]:
        raise ValueError("params musi zawierać 'prefix'")
    return cfg


def _fmt(value) -> str:
    """Listy -> 'a,b,c'; reszta -> str. Bez spacji (wymóg wywołania r.avaflow)."""
    if isinstance(value, (list, tuple)):
        return ",".join(str(v) for v in value)
    return str(value)


def build_start_script(cfg: dict) -> str:
    """Zwraca tekst skryptu startowego (import rastrów, region, r.avaflow, region domyślny)."""
    if 'hentrmax' in cfg['params']:
        raise ValueError('r.avaflow 4.0G Rev.7 ignoruje hentrmax; użyj hentrmax1 i rhentrmax1')
    lines = ["#!/bin/sh", "set -eu", "# Wygenerowano przez ravaflow_tools.runner - nie edytuj ręcznie", ""]
    for grass_name, file_path in cfg["inputs"].items():
        override = " -o" if cfg.get("allow_projection_override", False) else ""
        lines.append(f"r.in.gdal{override} --overwrite input={shlex.quote(str(file_path))} output={shlex.quote(grass_name)}")
    lines += ["", f"g.region -s rast={shlex.quote(cfg.get('region_from', 'elevation'))}", ""]

    # opcjonalne polecenia pomocnicze (np. r.lakefill, r.mapcalc) przed modelem
    for cmd in cfg.get("pre_commands", []):
        lines.append(cmd)
    if cfg.get("pre_commands"):
        lines.append("")

    flags = cfg.get("flags", "").strip()
    args = [shlex.quote(f"{k}={_fmt(v)}") for k, v in cfg["params"].items()]
    call = " ".join(["r.avaflow.40G"] + ([flags] if flags else []) + args)
    lines += [call, "", "g.region -d", ""]
    return "\n".join(lines)


def parse_evaluation(path: str | Path) -> dict:
    """Parsuje [prefix]_evaluation.txt: linie 'opis<TAB>wartość'.

    Zwraca słownik {pełny opis: wartość(float)}. Procenty zamieniane na ułamki
    (16.7% -> 0.167). -9999 oznacza brak danych (zostawiamy bez zmian).
    """
    out: dict[str, float] = {}
    for raw in Path(path).read_text(encoding="utf-8").splitlines():
        if "\t" not in raw:
            continue
        label, value = raw.rsplit("\t", 1)
        value = value.strip()
        try:
            if value.endswith("%"):
                out[label.strip()] = float(value[:-1]) / 100.0
            else:
                out[label.strip()] = float(value)
        except ValueError:
            continue
    return out


def run(cfg: dict, quiet: bool = False) -> dict:
    """Uruchamia jeden przebieg. Zwraca słownik z czasami, ścieżkami i oceną."""
    workdir = Path(cfg["workdir"]).expanduser().resolve()
    prefix = cfg["params"]["prefix"]
    if not re.fullmatch(r"[A-Za-z][A-Za-z0-9_]*", prefix):
        raise ValueError("prefix: wymagane litery/cyfry/podkreślenia, pierwsza litera")
    mapset = str(Path(cfg["grass_mapset"]).expanduser())
    workdir.mkdir(parents=True, exist_ok=True)
    if (workdir / f"{prefix}_results").exists():
        raise FileExistsError(f"Wyniki {prefix} już istnieją; użyj nowego katalogu lub prefiksu")

    script_path = workdir / f"{prefix}.run.sh"
    script_path.write_text(build_start_script(cfg), encoding="utf-8")
    log_path = workdir / f"{prefix}.run.log"

    lock = Path(mapset) / ".gislock"
    if lock.exists():
        raise RuntimeError(
            f"Mapset {mapset} jest zablokowany ({lock}). Sprawdź `pgrep -a grass`: jeśli "
            f"coś działa, poczekaj lub zamknij to; jeśli nic nie działa, usuń plik blokady."
        )

    if not quiet:
        print(f"[{prefix}] start (może potrwać kilka minut), log: {log_path}", flush=True)
    t0 = time.perf_counter()
    tempdir = workdir / '.grass_tmp'
    tempdir.mkdir(exist_ok=True)
    with open(log_path, "w", encoding="utf-8") as log:
        proc = subprocess.Popen(
            [cfg.get("grass_cmd", "grass"), mapset, "--exec", "sh", script_path.name],
            cwd=workdir, stdout=log, stderr=subprocess.STDOUT, text=True,
            stdin=subprocess.DEVNULL,  # GRASS nie może czekać na terminal (stan 'T')
            start_new_session=True,
            env={**grass_environment(), 'TMPDIR': str(tempdir)},
        )
        try:
            proc.wait(timeout=cfg.get("timeout_seconds"))
        except BaseException:
            terminate_group(proc)
            raise
    wall = time.perf_counter() - t0

    results_dir = workdir / f"{prefix}_results"
    files_dir = results_dir / f"{prefix}_files"
    result = {
        "prefix": prefix,
        "returncode": proc.returncode,
        "wall_seconds": wall,
        "log": str(log_path),
        "results_dir": str(results_dir),
        "ascii_dir": str(results_dir / f"{prefix}_ascii"),
    }
    time_file = files_dir / f"{prefix}_time.txt"
    if time_file.exists():
        result["model_seconds"] = float(time_file.read_text().split()[0])
    eval_file = files_dir / f"{prefix}_evaluation.txt"
    if eval_file.exists():
        result["evaluation"] = parse_evaluation(eval_file)

    if proc.returncode != 0 or not results_dir.exists():
        raise RuntimeError(
            f"Przebieg {prefix} nie powiódł się (kod {proc.returncode}). Zobacz log: {log_path}"
        )
    validate_success(files_dir, prefix)
    normalize_output_headers(results_dir / f'{prefix}_ascii', prefix)
    if not quiet:
        print(f"[{prefix}] OK, czas całkowity {wall:.1f} s, czas modelu {result.get('model_seconds', float('nan')):.1f} s")
        csi = result.get("evaluation", {}).get("Critical success index CSI")
        if csi is not None:
            print(f"[{prefix}] CSI (Jaccard) wg r.avaflow: {csi:.3f}")
    return result


def grass_environment() -> dict:
    """System GRASS must not load Python/Qt/Conda libraries from the caller."""
    env = dict(os.environ)
    for name in ("PYTHONHOME", "PYTHONPATH", "LD_LIBRARY_PATH", "QT_PLUGIN_PATH"):
        env.pop(name, None)
    env["PATH"] = "/usr/bin:/bin:" + env.get("PATH", "")
    env["GRASS_PYTHON"] = "/usr/bin/python3"
    env["GRASS_GUI"] = "text"
    env["OMP_NUM_THREADS"] = "1"
    return env


def terminate_group(proc) -> None:
    """Stop GRASS and its solver, including on timeout or Ctrl+C."""
    for sig in (signal.SIGTERM, signal.SIGKILL):
        try:
            os.killpg(proc.pid, sig)
        except ProcessLookupError:
            break
        try:
            proc.wait(timeout=5)
        except subprocess.TimeoutExpired:
            pass
    proc.wait()


def validate_success(files_dir: Path, prefix: str) -> None:
    # r.avaflow can exit with code 0 even after numerical failure.
    marker = files_dir / f"{prefix}_nout1.txt"
    lines = marker.read_text().splitlines() if marker.exists() else []
    if len(lines) < 2 or lines[1].strip() != "1":
        raise RuntimeError(f"Brak potwierdzenia poprawnego zakończenia solvera: {marker}")


def read_output_ascii(path):
    """Read a GRASS r.avaflow 4.0G Rev.7 raster with its actual corner coordinates.

    WITHGRASS foutasc writes cellhd.west/south under xllcenter/yllcenter.
    This workaround is specific to this solver, never generic ESRI ASCII input.
    """
    from . import metrics as M
    data, meta = M.read_ascii(path)
    with open(path) as f:
        keys = {f.readline().split()[0].lower() for _ in range(6)}
    if 'xllcenter' in keys:
        meta['xll'] += meta['cellsize'] / 2
    if 'yllcenter' in keys:
        meta['yll'] += meta['cellsize'] / 2
    return data, meta


def normalize_output_headers(directory: Path, prefix: str) -> None:
    """Fix the Rev.7 GRASS center/corner label in newly generated solver rasters.

    Equal-length substitutions preserve the entire numeric raster without copying it.
    Only this run's output rasters are modified, never source data.
    """
    for path in directory.glob(f'{prefix}_*.asc'):
        with path.open('r+b') as f:
            header = b''.join(f.readline() for _ in range(6))
            corrected = header.replace(b'xllcenter', b'xllcorner').replace(b'yllcenter', b'yllcorner')
            if corrected != header:
                f.seek(0)
                f.write(corrected)


if __name__ == "__main__":
    if len(sys.argv) != 2:
        sys.exit("Użycie: python -m ravaflow_tools.runner configs/plik.yaml")
    run(load_config(sys.argv[1]))
