#!/usr/bin/env python3
"""Unit tests for geometry and neighbourhood operations."""

import numpy as np

from sciddica_ss3 import shifted, neighbour_stack, cell_centres, hex_vertices, NEIGHBOUR_OFFSETS


class TestShifted:
    def test_shift_no_wrap(self):
        arr = np.array([[1, 2, 3], [4, 5, 6], [7, 8, 9]])
        result = shifted(arr, 1, 0, fill_value=0.0)
        expected = np.array([[2, 3, 0], [5, 6, 0], [8, 9, 0]])
        assert np.array_equal(result, expected)

    def test_shift_up(self):
        arr = np.array([[1, 2], [3, 4]])
        result = shifted(arr, 0, -1, fill_value=-1)
        expected = np.array([[-1, -1], [1, 2]])
        assert np.array_equal(result, expected)

    def test_shift_all_zero_outside(self):
        arr = np.ones((2, 2))
        result = shifted(arr, 5, 5, fill_value=0.0)
        assert np.allclose(result, 0.0)

    def test_shift_invalid_ndim(self):
        try:
            shifted(np.array([1, 2, 3]), 0, 0)
            assert False, "Should raise ValueError"
        except ValueError:
            pass


class TestNeighbourStack:
    def test_shape(self):
        arr = np.ones((5, 6))
        stack = neighbour_stack(arr)
        assert stack.shape == (6, 5, 6)

    def test_covers_all_offsets(self):
        arr = np.arange(9).reshape(3, 3)
        stack = neighbour_stack(arr, fill_value=0.0)
        assert len(NEIGHBOUR_OFFSETS) == 6
        for i in range(6):
            assert stack[i].shape == (3, 3)


class TestCellCentres:
    def test_positive_dimensions(self):
        x, y = cell_centres(3, 4, 1.0)
        assert x.shape == (3, 4)
        assert y.shape == (3, 4)

    def test_invalid_dimensions(self):
        try:
            cell_centres(0, 5, 1.0)
            assert False
        except ValueError:
            pass

    def test_spacing_increases_with_side(self):
        x1, y1 = cell_centres(2, 2, 1.0)
        x2, y2 = cell_centres(2, 2, 2.0)
        # Distance between cells scales with side length
        d1 = np.hypot(x1[0, 1] - x1[0, 0], y1[0, 1] - y1[0, 0])
        d2 = np.hypot(x2[0, 1] - x2[0, 0], y2[0, 1] - y2[0, 0])
        assert abs(d2 / d1 - 2.0) < 1e-10


class TestHexVertices:
    def test_shape(self):
        x, y = cell_centres(2, 3, 1.0)
        v = hex_vertices(x, y, 1.0)
        assert v.shape == (6, 6, 2)

    def test_six_vertices_per_hex(self):
        x = np.array([[0.0]])
        y = np.array([[0.0]])
        v = hex_vertices(x, y, 1.0)
        assert v.shape[1] == 6
        assert v.shape[2] == 2
