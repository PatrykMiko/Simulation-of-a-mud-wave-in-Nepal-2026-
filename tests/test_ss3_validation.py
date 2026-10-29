"""Output validation against theoretical/empirical expectations."""

import numpy as np

from sciddica_ss3 import SciddicaSS3, ModelParameters, make_albano_like, jaccard_fitness


def test_cell_area_and_spacing():
    p = ModelParameters(pc=12.0)
    assert abs(p.cell_area - 1.5 * np.sqrt(3.0) * 144.0) < 1e-12
    assert abs(p.centre_spacing - np.sqrt(3.0) * 12.0) < 1e-12


def test_submerged_gravity():
    p = ModelParameters(rho=2000.0, rho_water=1000.0)
    assert abs(p.submerged_gravity - 9.81 * 1000.0 / 2000.0) < 1e-12


def test_initial_scenario_mass():
    model, x, y = make_albano_like(n=40)
    # Initial debris volume should be positive
    assert model.total_volume() > 0
    # Solid balance should hold
    area = model.params.cell_area
    solid = (model.thickness.sum() + model.soil_depth.sum()) * area
    assert solid > 0


def test_jaccard_fitness():
    a = np.array([True, False, True])
    b = np.array([True, True, False])
    f = jaccard_fitness(a, b)
    assert abs(f - 1.0 / 3.0) < 1e-12
    assert jaccard_fitness(a, a) == 1.0
    assert jaccard_fitness(np.array([False]), np.array([False])) == 1.0


def test_no_negative_state_after_run():
    alt = np.ones((12, 12)) * 80.0
    th = np.zeros((12, 12))
    th[5, 5] = 3.0
    soil = np.ones((12, 12)) * 1.0
    model = SciddicaSS3(alt, th, soil)
    model.run(20)
    assert np.all(model.thickness >= 0)
    assert np.all(model.soil_depth >= 0)
    assert np.all(model.kinetic_head >= 0)
    assert np.all(np.isfinite(model.thickness))
    assert np.all(np.isfinite(model.kinetic_head))
