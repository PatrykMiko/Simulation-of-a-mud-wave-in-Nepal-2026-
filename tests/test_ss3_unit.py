"""Unit tests for core components."""

import numpy as np

from sciddica_ss3 import (
    MediumParameters,
    ModelParameters,
    SciddicaSS3,
    minimization_algorithm,
    momentum_elevation_correction,
    neighbour_stack,
    shifted,
)


def test_shifted_no_wrap():
    arr = np.array([[1.0, 2.0], [3.0, 4.0]])
    out = shifted(arr, 1, 0, fill_value=0.0)
    assert np.allclose(out, np.array([[2.0, 0.0], [4.0, 0.0]]))


def test_neighbour_stack_shape():
    arr = np.array([[1.0, 2.0], [3.0, 4.0]])
    ns = neighbour_stack(arr)
    assert ns.shape == (6, 2, 2)


def test_medium_params_validation():
    try:
        MediumParameters(adh=-1.0)
        assert False
    except ValueError:
        pass
    try:
        MediumParameters(ptd=1.5)
        assert False
    except ValueError:
        pass


def test_model_params_validation():
    try:
        ModelParameters(pc=0.0)
        assert False
    except ValueError:
        pass
    try:
        ModelParameters(pml=1.5)
        assert False
    except ValueError:
        pass


def test_ma_basic():
    movable = np.array([[1.0]])
    centre = np.array([[0.0]])
    neighbours = np.array([[[0.0]], [[0.0]], [[0.0]], [[0.0]], [[0.0]], [[0.0]]])
    allowed = np.ones((6, 1, 1), dtype=bool)
    out = minimization_algorithm(movable, centre, neighbours, allowed)
    assert out.shape == (6, 1, 1)
    assert np.all(out >= 0)


def test_momentum_correction_shape():
    th = np.ones((3, 3))
    px = np.ones((3, 3))
    py = np.zeros((3, 3))
    p = ModelParameters()
    corr = momentum_elevation_correction(th, px, py, p)
    assert corr.shape == (6, 3, 3)


def test_sciddica_init_and_state():
    alt = np.zeros((5, 5))
    th = np.zeros((5, 5))
    soil = np.ones((5, 5))
    model = SciddicaSS3(alt, th, soil)
    assert model.height == 5
    assert model.width == 5
    assert np.allclose(model.total_volume(), 0.0)
    assert model.steps_completed == 0
