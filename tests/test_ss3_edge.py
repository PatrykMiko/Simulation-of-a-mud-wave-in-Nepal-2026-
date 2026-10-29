"""Edge-case and error handling tests."""

import numpy as np
import pytest

from sciddica_ss3 import SciddicaSS3, ModelParameters, shifted, minimization_algorithm


def test_shifted_invalid_dim():
    try:
        shifted(np.array([1, 2, 3]), 0, 0)
        assert False
    except ValueError:
        pass


def test_ma_invalid_shapes():
    try:
        minimization_algorithm(np.array([1.0]), np.array([0.0]), np.zeros((6, 1, 1)), np.ones((6, 1, 1), dtype=bool))
        assert False
    except ValueError:
        pass


def test_ma_negative_mass_raises():
    try:
        minimization_algorithm(np.array([[-1.0]]), np.array([[0.0]]), np.zeros((6, 1, 1)), np.ones((6, 1, 1), dtype=bool))
        assert False
    except ValueError:
        pass


def test_sciddica_init_errors():
    try:
        SciddicaSS3(np.zeros((3,)), np.zeros((3,)), np.zeros((3,)))
        assert False
    except ValueError:
        pass
    try:
        SciddicaSS3(np.zeros((2, 2)), np.array([[1.0, -1.0], [0.0, 0.0]]), np.zeros((2, 2)))
        assert False
    except ValueError:
        pass
    try:
        SciddicaSS3(np.array([np.nan]), np.array([0.0]), np.array([0.0]))
        assert False
    except ValueError:
        pass


def test_run_negative_steps():
    alt = np.zeros((3, 3))
    th = np.zeros((3, 3))
    soil = np.zeros((3, 3))
    model = SciddicaSS3(alt, th, soil)
    try:
        model.run(-1)
        assert False
    except ValueError:
        pass


def test_water_level_invalid():
    try:
        SciddicaSS3(np.zeros((3, 3)), np.zeros((3, 3)), np.zeros((3, 3)), water_level=np.nan)
        assert False
    except ValueError:
        pass


def test_empty_like_behavior():
    alt = np.array([[5.0]])
    th = np.array([[0.0]])
    soil = np.array([[0.0]])
    model = SciddicaSS3(alt, th, soil)
    model.step()
    assert model.total_volume() >= 0
