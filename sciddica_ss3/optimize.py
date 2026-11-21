"""Optuna-based hyperparameter optimization for SCIDDICA-SS3 on Nepal data."""

import json
import logging
from dataclasses import dataclass
from pathlib import Path
from typing import Any

import numpy as np
import optuna

from .automaton import SciddicaSS3
from .constants import MIN_THICKNESS
from .scenario import jaccard_fitness
from .state import MediumParameters, ModelParameters

logger = logging.getLogger(__name__)


def _trial_callback(study: optuna.Study, trial: optuna.Trial) -> None:
    value = trial.value if trial.value is not None else float("nan")
    logger.info(
        "trial=%d value=%.6f best=%.6f params=%s",
        trial.number,
        value,
        study.best_value if study.best_trial is not None else float("nan"),
        trial.params,
    )


@dataclass(frozen=True)
class NepalScene:
    altitude: np.ndarray
    thickness: np.ndarray
    soil_depth: np.ndarray
    domain: np.ndarray
    reference: np.ndarray
    pc: float
    water_level: float | None
    source_cellsize: float
    source_dir: str

    def __post_init__(self) -> None:
        shapes = [
            np.asarray(self.altitude).shape,
            np.asarray(self.thickness).shape,
            np.asarray(self.soil_depth).shape,
            np.asarray(self.domain).shape,
            np.asarray(self.reference).shape,
        ]
        if len({s for s in shapes}) != 1:
            raise ValueError("altitude, thickness, soil_depth, domain, reference must share shape")
        if not np.isfinite(self.pc) or self.pc <= 0:
            raise ValueError("pc must be finite and > 0")
        ref_arr = np.asarray(self.reference)
        dom_arr = np.asarray(self.domain)
        if np.any(ref_arr & ~dom_arr):
            raise ValueError("reference must be a subset of domain")


def validate_ss3_conversion(out_dir: Path) -> dict:
    out_dir = Path(out_dir)
    if not out_dir.exists() or not out_dir.is_dir():
        raise FileNotFoundError(f"{out_dir} does not exist or is not a directory; run ravaflow_tools.to_ss3 first")

    npz_path = out_dir / "ss3_arrays.npz"
    meta_path = out_dir / "ss3_meta.json"
    if not npz_path.exists():
        raise FileNotFoundError(f"missing {npz_path}; run ravaflow_tools.to_ss3 first")
    if not meta_path.exists():
        raise FileNotFoundError(f"missing {meta_path}; run ravaflow_tools.to_ss3 first")

    try:
        data = np.load(npz_path)
    except Exception as exc:
        raise ValueError(f"failed to load {npz_path}: {exc}") from exc

    required = ["altitude", "thickness", "soil_depth", "domain", "reference"]
    missing = [k for k in required if k not in data]
    if missing:
        raise ValueError(f"missing keys in {npz_path}: {missing}")

    altitude = np.asarray(data["altitude"])
    thickness = np.asarray(data["thickness"])
    soil_depth_arr = np.asarray(data["soil_depth"])
    domain = np.asarray(data["domain"])
    reference = np.asarray(data["reference"])

    shapes = [altitude.shape, thickness.shape, soil_depth_arr.shape, domain.shape, reference.shape]
    if len({s for s in shapes}) != 1 or len(shapes[0]) != 2:
        raise ValueError(f"arrays must share the same 2-D shape; got {[s for s in shapes]}")

    if not np.all(np.isfinite(altitude)):
        raise ValueError("altitude must be finite")
    if altitude.size == 0:
        raise ValueError("altitude must be non-empty")
    if not np.all(np.isfinite(thickness)) or np.any(thickness < 0):
        raise ValueError("thickness must be finite and >= 0")
    if not np.all(np.isfinite(soil_depth_arr)) or np.any(soil_depth_arr < 0):
        raise ValueError("soil_depth must be finite and >= 0")

    domain_bool = domain.astype(bool)
    reference_bool = reference.astype(bool)
    if np.any(reference_bool & ~domain_bool):
        raise ValueError("reference must be a subset of domain")
    if not np.any(domain_bool):
        raise ValueError("domain contains no cells")

    try:
        with open(meta_path) as f:
            meta = json.load(f)
    except Exception as exc:
        raise ValueError(f"failed to load {meta_path}: {exc}") from exc

    if "cellsize" not in meta:
        raise ValueError(f"missing 'cellsize' in {meta_path}")
    pc = float(meta["cellsize"])
    if not np.isfinite(pc) or pc <= 0:
        raise ValueError("cellsize must be finite and > 0")
    water_level = meta.get("water_level")
    if water_level is not None:
        water_level = float(water_level)
        if not np.isfinite(water_level):
            raise ValueError("water_level must be finite or None")
    source_cellsize = float(meta.get("source_cellsize", np.nan))
    source_dir = str(meta.get("source_dir", ""))

    return {
        "altitude": altitude.astype(np.float64),
        "thickness": thickness.astype(np.float64),
        "soil_depth": soil_depth_arr.astype(np.float64),
        "domain": domain_bool,
        "reference": reference_bool,
        "pc": pc,
        "water_level": water_level,
        "source_cellsize": source_cellsize,
        "source_dir": source_dir,
    }


def load_nepal_scene(root: Path | str = "data/nepal", source: str = "alos", data_dir: Path | str | None = None) -> NepalScene:
    if data_dir is not None:
        out_dir = Path(data_dir)
    else:
        root_p = Path(root)
        out_dir = root_p / f"{source}_ss3"
    data = validate_ss3_conversion(out_dir)
    logger.info("loaded SS3 scene from %s (shape=%s, pc=%.3f)", out_dir, data["altitude"].shape, data["pc"])
    return NepalScene(
        altitude=data["altitude"],
        thickness=data["thickness"],
        soil_depth=data["soil_depth"],
        domain=data["domain"],
        reference=data["reference"],
        pc=data["pc"],
        water_level=data["water_level"],
        source_cellsize=data["source_cellsize"],
        source_dir=data["source_dir"],
    )


def params_from_trial(trial: optuna.Trial, pc: float) -> ModelParameters:
    pt = 2.0
    rho = 2000.0

    air = MediumParameters(
        adh=trial.suggest_float("adh", 0.015, 0.025),
        pf=trial.suggest_float("pf", 0.0, 0.5),
        pmt=trial.suggest_float("pmt", 0.0, 5.0),
        per=trial.suggest_float("per", 0.0, 0.2),
        ptd=trial.suggest_float("ptd", 0.1, 0.7),
        ped=trial.suggest_float("ped", 0.0, 0.5),
        shift_friction=trial.suggest_float("shift_friction", 0.02, 0.3),
    )

    return ModelParameters(
        pc=pc,
        pt=pt,
        rho=rho,
        air=air,
        water=MediumParameters(),
        pml=trial.suggest_float("pml", 0.0, 0.5),
        pwr=trial.suggest_float("pwr", 0.0, 1.0),
        alpha_scale=trial.suggest_float("alpha_scale", 0.5, 2.0),
    )


def run_model(
    scene: NepalScene,
    params: ModelParameters,
    steps: int | None = None,
    soil_depth: float = 1.0,
    stable_tolerance: float | None = 1e-4,
    stable_window: int = 5,
    max_steps: int = 10000,
    log_every: int = 100,
    check_every: int = 100,
    use_float32: bool = False,
    use_numba: bool | None = None,
    use_cpp: bool | None = None,
) -> SciddicaSS3:
    from .dtypes import compute_dtype
    dtype = compute_dtype(use_float32)
    soil = np.asarray(scene.soil_depth, dtype=dtype).copy()
    model = SciddicaSS3(
        altitude=np.asarray(scene.altitude, dtype=dtype),
        thickness=np.asarray(scene.thickness, dtype=dtype),
        soil_depth=soil,
        params=params,
        water_level=scene.water_level,
        use_float32=use_float32,
        use_numba=use_numba,
        use_cpp=use_cpp,
    )

    previous_mass = float(model.total_volume())
    previous_mask = model.thickness > MIN_THICKNESS
    stable_count = 0
    step = 0

    while steps is None or step < steps:
        if step >= max_steps:
            break
        model.step()
        step += 1

        if check_every and step % check_every == 0:
            current_mass = float(model.total_volume())
            current_mask = model.thickness > MIN_THICKNESS

            mass_change = abs(current_mass - previous_mass) / max(previous_mass, 1e-12)
            extent_change = 1.0 - jaccard_fitness(previous_mask, current_mask)
            current_iou = jaccard_fitness(scene.reference, current_mask)

            if log_every and step % log_every == 0:
                logger.info(
                    "step=%d mass=%.6f mass_change=%.3e extent_change=%.3e iou=%.6f active_cells=%d",
                    step,
                    current_mass,
                    mass_change,
                    extent_change,
                    current_iou,
                    int(current_mask.sum()),
                )

            if stable_tolerance is not None and mass_change < stable_tolerance and extent_change < stable_tolerance:
                stable_count += 1
                if stable_count >= stable_window:
                    break
            else:
                stable_count = 0

            previous_mass = current_mass
            previous_mask = current_mask

    return model


def run_simulation(*args, **kwargs) -> np.ndarray:
    model = run_model(*args, **kwargs)
    return model.thickness > MIN_THICKNESS


def iou_for_params(
    scene: NepalScene,
    params: ModelParameters,
    steps: int | None = None,
    stable_tolerance: float | None = 1e-4,
    use_float32: bool = True,
    use_numba: bool | None = True,
    use_cpp: bool | None = None,
) -> float:
    simulated = run_simulation(scene, params, steps, stable_tolerance=stable_tolerance, use_float32=use_float32, use_numba=use_numba, use_cpp=use_cpp)
    return jaccard_fitness(scene.reference, simulated)


def save_best_visualization(scene: NepalScene, model: SciddicaSS3, trial: optuna.Trial, output_dir: Path | str = "best_results") -> Path:
    import matplotlib.pyplot as plt
    from matplotlib.patches import Patch

    out = Path(output_dir)
    out.mkdir(parents=True, exist_ok=True)
    path = out / f"best_trial_{trial.number:04d}_iou.png"

    simulated = model.thickness > MIN_THICKNESS
    ref = np.asarray(scene.reference, dtype=bool)
    sim = np.asarray(simulated, dtype=bool)

    tp = ref & sim  # true positive
    fp = ~ref & sim  # false positive
    fn = ref & ~sim  # false negative
    tn = ~ref & ~sim  # true negative

    comp = np.ones(ref.shape + (3,), dtype=float)  # white background
    comp[tp] = [0.0, 0.8, 0.0]  # green
    comp[fp] = [1.0, 0.0, 0.0]  # red
    comp[fn] = [0.0, 0.4, 1.0]  # blue

    iou = jaccard_fitness(ref, sim)
    fig, ax = plt.subplots(1, 1, figsize=(8, 6))
    ax.imshow(comp, interpolation="nearest")
    ax.set_title(f"Simulated vs Real - Trial {trial.number} (IoU={iou:.3f})")
    ax.axis("off")

    legend_elements = [
        Patch(facecolor="green", label="True Positive (overlap)"),
        Patch(facecolor="red", label="False Positive (sim only)"),
        Patch(facecolor="blue", label="False Negative (real only)"),
        Patch(facecolor="white", edgecolor="black", label="True Negative"),
    ]
    ax.legend(handles=legend_elements, loc="lower right", fontsize=8, framealpha=0.9)
    fig.tight_layout()
    fig.savefig(path, dpi=150)
    plt.close(fig)
    logger.info("saved best comparison map: %s", path)
    return path


def make_objective(scene: NepalScene, steps: int | None = None, stable_tolerance: float | None = 1e-4, output_dir: Path | str = "best_results", use_float32: bool = True, use_numba: bool | None = True, use_cpp: bool | None = None, check_every: int = 100):
    best = {"value": -1.0}

    def objective(trial: optuna.Trial) -> float:
        params = params_from_trial(trial, pc=scene.pc)
        try:
            model = run_model(scene, params, steps, stable_tolerance=stable_tolerance, use_float32=use_float32, use_numba=use_numba, use_cpp=use_cpp, check_every=check_every)
            simulated = model.thickness > MIN_THICKNESS
            iou = jaccard_fitness(scene.reference, simulated)
            trial.set_user_attr("iou", iou)
            if iou > best["value"]:
                best["value"] = iou
                save_best_visualization(scene, model, trial, output_dir)
            return iou
        except Exception as exc:  # pragma: no cover - defensive
            trial.set_user_attr("error", str(exc))
            return 0.0

    return objective


def optimize_nepal(
    data_root: Path | str = "data/nepal",
    source: str = "alos",
    data_dir: Path | str | None = None,
    steps: int | None = None,
    n_trials: int = 50,
    n_jobs: int = 1,
    study_name: str = "sciddica_ss3_nepal",
    stable_tolerance: float | None = 1e-4,
    storage: str | None = "sqlite:///sciddica_ss3_nepal.db",
    callbacks: bool = True,
    output_dir: Path | str = "best_results",
    use_float32: bool = True,
    use_numba: bool | None = None,
    use_cpp: bool | None = None,
    check_every: int = 300,
) -> optuna.Study:
    if not logging.getLogger().handlers:
        logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s %(message)s")
    scene = load_nepal_scene(data_root, source, data_dir=data_dir)
    study = optuna.create_study(
        study_name=study_name,
        direction="maximize",
        sampler=optuna.samplers.TPESampler(seed=42),
        storage=storage,
        load_if_exists=True,
    )
    callback_list = [_trial_callback] if callbacks else None
    study.optimize(
        make_objective(scene, steps, stable_tolerance, output_dir, use_float32=use_float32, use_numba=use_numba, use_cpp=use_cpp, check_every=check_every),
        n_trials=n_trials,
        n_jobs=n_jobs,
        callbacks=callback_list,
    )
    return study
if __name__ == "__main__":
    import argparse
    import logging

    logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s %(message)s")
    p = argparse.ArgumentParser(description="Optimize SCIDDICA-SS3 parameters on Nepal data")
    p.add_argument("root", nargs="?", default="data/nepal", help="Root directory containing <source>_ss3/")
    p.add_argument("--data-root", default=None, help="Root directory containing <source>_ss3/ (alias)")
    p.add_argument("--data-dir", default=None, help="Direct path to SS3 directory")
    p.add_argument("--source", default="alos", help="Source name")
    p.add_argument("--trials", "--n-trials", dest="trials", type=int, default=20, help="Number of optuna trials")
    p.add_argument("--n-jobs", type=int, default=1, help="Number of parallel Optuna workers")
    p.add_argument("--steps", type=int, default=None, help="Steps per trial (default: run to stability)")
    p.add_argument(
        "--use-float32",
        action=argparse.BooleanOptionalAction,
        default=True,
        help="Use float32 compute dtype (default: enabled)",
    )

    p.add_argument(
        "--use-numba",
        action=argparse.BooleanOptionalAction,
        default=True,
        help="Use numba kernels when available (default: enabled)",
    )
    p.add_argument(
        "--use-cpp",
        action=argparse.BooleanOptionalAction,
        default=True,
        help="Use C++ kernels when available (default: enabled)",
    )
    p.add_argument("--check-every", type=int, default=300, help="Interval for early-stopping check/logging")
    p.add_argument("--study-name", default="sciddica_ss3_nepal", help="Optuna study name")
    p.add_argument("--storage", default="sqlite:///sciddica_ss3_nepal.db", help="Optuna storage URL")
    p.add_argument("--output-dir", default="best_results", help="Directory for best-trial visualizations")
    a = p.parse_args()

    data_root = a.data_root if a.data_root is not None else a.root
    study = optimize_nepal(
        data_root=data_root,
        source=a.source,
        data_dir=a.data_dir,
        steps=a.steps,
        n_trials=a.trials,
        n_jobs=a.n_jobs,
        study_name=a.study_name,
        storage=a.storage,
        output_dir=a.output_dir,
        use_float32=a.use_float32,
        use_numba=a.use_numba,
        use_cpp=a.use_cpp,
        check_every=a.check_every,
    )
    print("best:", study.best_params, study.best_value)
