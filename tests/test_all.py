"""Testy automatyczne (pytest): metryki, walidator, generator, scenariusze, parser wyników.

Uruchomienie z katalogu głównego projektu:  python -m pytest -q
Testy nie wymagają GRASS-a ani r.avaflow (używają danych syntetycznych w formacie .asc).
"""
import sys
import tempfile
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from ravaflow_tools import metrics as M  # noqa: E402
from ravaflow_tools import runner, scenarios, synthetic  # noqa: E402
from ravaflow_tools.check_inputs import check_inputs  # noqa: E402


def _synthetic(tmp: str):
    out = Path(tmp) / "syn"
    synthetic.generate(out, fmt="asc", volume=2e6)
    return out


# ------------------------------------------------------------------ metryki
def test_jaccard_identical_and_disjoint():
    a = np.zeros((10, 10), bool)
    a[2:5, 2:5] = True
    ones = np.ones_like(a)
    assert M.score(a, a, ones, 30.0)["jaccard"] == 1.0
    assert M.score(a, ~a, ones, 30.0)["jaccard"] == 0.0


def test_known_counts():
    obs = np.zeros((4, 4), bool)
    obs[0, :] = True                      # 4 komórki
    mod = np.zeros((4, 4), bool)
    mod[0, :2] = True                     # 2 TP
    mod[1, :2] = True                     # 2 FP, więc 2 FN
    s = M.score(mod, obs, np.ones_like(obs), 10.0)
    assert (s["tp"], s["fp"], s["fn"]) == (2, 2, 2)
    assert abs(s["jaccard"] - 2 / 6) < 1e-12
    assert abs(s["area_ratio"] - 1.0) < 1e-12


def test_invalid_cells_are_ignored():
    obs = np.zeros((3, 3), bool)
    mod = np.ones((3, 3), bool)
    valid = np.zeros((3, 3), bool)
    valid[0, 0] = True
    s = M.score(mod, obs, valid, 1.0)
    assert s["fp"] == 1 and s["tn"] == 0


def test_ascii_roundtrip_and_resample():
    with tempfile.TemporaryDirectory() as tmp:
        arr = np.arange(12, dtype=float).reshape(3, 4)
        meta = {"ncols": 4, "nrows": 3, "xll": 100.0, "yll": 200.0, "cellsize": 10.0, "nodata": -9999.0}
        p = Path(tmp) / "a.asc"
        M.write_ascii(p, arr, meta)
        back, m2 = M.read_ascii(p)
        assert np.allclose(back, arr) and M.same_grid(meta, m2)
        coarse = {"ncols": 2, "nrows": 1, "xll": 100.0, "yll": 200.0, "cellsize": 20.0, "nodata": -9999.0}
        out = M.sample_to_grid(arr, meta, coarse)
        assert out.shape == (1, 2)


# ------------------------------------------------------------ walidator danych
def test_synthetic_passes_validator():
    with tempfile.TemporaryDirectory() as tmp:
        rep = check_inputs(_synthetic(tmp), volume=2e6, cellsize=30)
        assert rep.ok, rep.errors


def test_validator_catches_problems():
    with tempfile.TemporaryDirectory() as tmp:
        d = _synthetic(tmp)
        h, m = M.read_ascii(d / "release_height.asc")
        h[0, 0] = -9999.0                                   # nodata zamiast 0
        M.write_ascii(d / "release_height.asc", h, m)
        r, mr = M.read_ascii(d / "reference_extent.asc")
        r[50, 50] = 7.0                                     # niedozwolona wartość maski
        M.write_ascii(d / "reference_extent.asc", r, mr)
        rep = check_inputs(d, volume=5e6, cellsize=10)
        text = " ".join(rep.errors)
        assert not rep.ok
        for needle in ("rozmiar komórki", "bez wartości", "objętość", "reference_extent"):
            assert needle in text, needle


def test_validator_missing_file():
    with tempfile.TemporaryDirectory() as tmp:
        d = _synthetic(tmp)
        (d / "reference_extent.asc").unlink()
        assert any("reference_extent" in e for e in check_inputs(d).errors)


# ------------------------------------------------------------------ scenariusze
def _base_cfg(workdir):
    return {"workdir": str(workdir), "grass_mapset": "x",
            "inputs": {"syn_elev": "elevation.asc", "syn_release": "release_height.asc",
                       "syn_ref": "reference_extent.asc"},
            "params": {"prefix": "syn_base", "elevation": "syn_elev", "hrelease": "syn_release"}}


def test_scale_release_scales_volume():
    with tempfile.TemporaryDirectory() as tmp:
        d = _synthetic(tmp)
        cfg = scenarios.make_variant(_base_cfg(d), {"name": "half", "scale_release": 0.5})
        h0, _ = M.read_raster(d / "release_height.asc")
        h1, _ = M.read_raster(d / cfg["inputs"]["syn_release"])
        assert abs(h1.sum() / h0.sum() - 0.5) < 1e-6
        assert cfg["params"]["prefix"] == "syn_base_half"


def test_barrier_raises_terrain():
    with tempfile.TemporaryDirectory() as tmp:
        d = _synthetic(tmp)
        spec = {"name": "dam", "barrier": {"x1": 342700, "y1": 3115985, "x2": 344000, "y2": 3115985,
                                           "height": 40, "width": 60}}
        cfg = scenarios.make_variant(_base_cfg(d), spec)
        z0, _ = M.read_raster(d / "elevation.asc")
        z1, _ = M.read_raster(d / cfg["inputs"]["syn_elev"])
        diff = z1 - z0
        assert diff.max() > 39.9 and (diff > 0).sum() > 20


# ------------------------------------------------------- runner (bez GRASS-a)
def test_parse_evaluation():
    with tempfile.TemporaryDirectory() as tmp:
        p = Path(tmp) / "ev.txt"
        p.write_text("Criterion\tValue\n\nCritical success index CSI\t0.549\n"
                     "Rate of positive observations rOP\t16.7%\n"
                     "Ratio of deposit rLd\t-9999.00\n", encoding="utf-8")
        d = runner.parse_evaluation(p)
        assert d["Critical success index CSI"] == 0.549
        assert abs(d["Rate of positive observations rOP"] - 0.167) < 1e-9
        assert d["Ratio of deposit rLd"] == -9999.0


def test_start_script_matches_cascadia_example():
    cfg = runner.load_config(Path(__file__).resolve().parents[1] / "configs" / "ca_exp2.yaml")
    script = runner.build_start_script(cfg)
    call = [ln for ln in script.splitlines() if ln.startswith("r.avaflow.40G")][0]
    assert "prefix=ca_exp2" in call and "friction=35,0,0,35,0,0,0,0,0.05" in call
    assert "r.in.gdal -o --overwrite input=DATA/ca_elev.tif output=ca_elev" in script

def test_grid_differences_reports_what_differs():
    """Przypadek z syn_base: r.avaflow zwraca siatkę 201x301 przesuniętą względem referencji 200x300."""
    model = {"ncols": 201, "nrows": 301, "xll": 339975.0, "yll": 3109965.0, "cellsize": 30.0, "nodata": -9999.0}
    ref = {"ncols": 200, "nrows": 300, "xll": 340000.0, "yll": 3110000.0, "cellsize": 30.0, "nodata": -9999.0}
    assert not M.same_grid(model, ref)
    text = " ".join(M.grid_differences(model, ref))
    assert "201x301" in text and "200x300" in text
    assert "dx=+5.0" in text and "dy=-5.0" in text
    assert "rozmiar komórki" not in text          # komórki są takie same, więc nie wolno o nich pisać
    assert M.grid_differences(ref, dict(ref)) == []
    coarse = {**ref, "cellsize": 60.0}
    assert "rozmiar komórki: 30.0 m vs 60.0 m" in " ".join(M.grid_differences(ref, coarse))


def test_shifted_grid_reference_is_resampled_and_edges_ignored():
    ref_meta = {"ncols": 200, "nrows": 300, "xll": 340000.0, "yll": 3110000.0, "cellsize": 30.0, "nodata": -9999.0}
    model_meta = {"ncols": 201, "nrows": 301, "xll": 339975.0, "yll": 3109965.0, "cellsize": 30.0, "nodata": -9999.0}
    ref = np.zeros((300, 200))
    ref[100:150, 80:120] = 1.0                    # 50 x 40 = 2000 komórek zasięgu
    out = M.sample_to_grid(ref, ref_meta, model_meta)
    assert out.shape == (301, 201)
    obs, valid = M.reference_masks(out, {**model_meta, "nodata": ref_meta["nodata"]})
    assert obs.sum() == 50 * 40                   # zasięg nie ginie przy przesunięciu o 5 m
    assert valid.sum() == 200 * 300               # skrajna kolumna i wiersz modelu są poza referencją
