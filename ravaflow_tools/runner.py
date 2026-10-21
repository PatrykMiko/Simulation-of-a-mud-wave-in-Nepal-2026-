"""Uruchamianie r.avaflow 4.0G z pliku konfiguracji YAML.

Z konfiguracji generuje skrypt startowy (taki jak ca.avaflow.start.sh z paczki
Cascadia), uruchamia go przez `grass <mapset> --exec sh <skrypt>`, mierzy czas
i wczytuje plik *_evaluation.txt.

Użycie (z katalogu głównego projektu, poza sesją GRASS i najlepiej bez condy):
    python -m ravaflow_tools.runner configs/ca_exp2.yaml
"""
from __future__ import annotations

import subprocess
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
    lines = ["# Wygenerowano przez ravaflow_tools.runner - nie edytuj ręcznie", ""]
    for grass_name, file_path in cfg["inputs"].items():
        lines.append(f"r.in.gdal -o --overwrite input={file_path} output={grass_name}")
    lines += ["", f"g.region -s rast={cfg.get('region_from', 'elevation')}", ""]

    # opcjonalne polecenia pomocnicze (np. r.lakefill, r.mapcalc) przed modelem
    for cmd in cfg.get("pre_commands", []):
        lines.append(cmd)
    if cfg.get("pre_commands"):
        lines.append("")

    flags = cfg.get("flags", "").strip()
    args = [f"{k}={_fmt(v)}" for k, v in cfg["params"].items()]
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
    mapset = str(Path(cfg["grass_mapset"]).expanduser())

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
    with open(log_path, "w", encoding="utf-8") as log:
        proc = subprocess.run(
            ["grass", mapset, "--exec", "sh", script_path.name],
            cwd=workdir, stdout=log, stderr=subprocess.STDOUT, text=True,
            stdin=subprocess.DEVNULL,  # GRASS nie może czekać na terminal (stan 'T')
            timeout=cfg.get("timeout_seconds"),
        )
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
    if not quiet:
        print(f"[{prefix}] OK, czas całkowity {wall:.1f} s, czas modelu {result.get('model_seconds', float('nan')):.1f} s")
        csi = result.get("evaluation", {}).get("Critical success index CSI")
        if csi is not None:
            print(f"[{prefix}] CSI (Jaccard) wg r.avaflow: {csi:.3f}")
    return result


if __name__ == "__main__":
    if len(sys.argv) != 2:
        sys.exit("Użycie: python -m ravaflow_tools.runner configs/plik.yaml")
    run(load_config(sys.argv[1]))
