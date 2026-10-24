"""Integration tests for complex sequences."""

import numpy as np

from sciddica_ss3 import SciddicaSS3, ModelParameters


def test_single_step_conserves_solid_within_tol():
    alt = np.array([[10.0, 10.0], [10.0, 10.0]])
    th = np.array([[1.0, 0.0], [0.0, 0.0]])
    soil = np.array([[0.5, 0.5], [0.5, 0.5]])
    p = ModelParameters()
    model = SciddicaSS3(alt, th, soil, p)
    area = p.cell_area
    initial_solid = (model.thickness.sum() + model.soil_depth.sum()) * area
    model.step()
    final_solid = (model.thickness.sum() + model.soil_depth.sum()) * area
    assert abs((final_solid + model.mass_lost) - initial_solid) < 1e-9


def test_multiple_steps_run():
    alt = np.ones((10, 10)) * 100.0
    th = np.zeros((10, 10))
    th[4:6, 4:6] = 2.0
    soil = np.ones((10, 10)) * 0.1
    model = SciddicaSS3(alt, th, soil)
    model.run(5)
    assert model.steps_completed == 5
    assert model.time == 5 * model.params.pt


def test_volume_and_energy_finite():
    alt = np.ones((8, 8)) * 50.0
    th = np.zeros((8, 8))
    th[3, 3] = 1.5
    soil = np.ones((8, 8)) * 0.2
    model = SciddicaSS3(alt, th, soil)
    model.run(10)
    v = model.total_volume()
    ke = model.kinetic_energy()
    te = model.total_energy()
    assert np.isfinite(v)
    assert np.isfinite(ke)
    assert np.isfinite(te)
    assert v >= 0
    assert ke >= 0
    assert te >= 0
