from BasicTools.Numerical.Pod import *

import numpy as np
import numpy.testing as npt


def test_pod():
    snapshots = np.zeros((3, 3), dtype=np.double)
    snapshots[0,0] = 2.0
    snapshots[1,1] = 1.0
    snapshots[2,2] = 3.0

    full_basis = np.zeros((3, 3), dtype=np.double)
    full_basis[0,2] = 1.0
    full_basis[1,0] = 1.0
    full_basis[2,1] = 1.0

    actual = pod_basis(snapshots.tolist())
    expected = full_basis
    npt.assert_allclose(actual, expected)

    actual = pod_basis(snapshots)
    expected = full_basis
    npt.assert_allclose(actual, expected)

    actual = pod_basis(snapshots, rank_max=2)
    expected = full_basis[0:2, :]
    npt.assert_allclose(expected, actual)

    actual = pod_basis(snapshots, truncation_tol=1.0)
    expected = full_basis[0:0, :]
    npt.assert_allclose(actual, expected)

    actual = pod_basis(snapshots, truncation_tol=np.sqrt(1.0/3.0))
    expected = full_basis[0:2, :]
    npt.assert_allclose(actual, expected)

    actual = pod_basis(snapshots, alternate_truncation_tol=0.3334)
    expected = full_basis[0:2, :]
    npt.assert_allclose(actual, expected)

    actual = pod_basis(snapshots, alternate_truncation_tol=0.6667)
    expected = full_basis[0:1, :]
    npt.assert_allclose(actual, expected)

    actual = pod_basis(snapshots, rank_max=2, truncation_tol=np.sqrt(0.5))
    expected = full_basis[0:1, :]
    npt.assert_allclose(actual, expected)

    actual = pod_basis(snapshots, rank_max=1, truncation_tol=np.sqrt(0.5))
    expected = full_basis[0:1, :]
    npt.assert_allclose(actual, expected)

    actual = pod_basis(snapshots, rank_max=1, truncation_tol=np.sqrt(1.0/3.0))
    expected = full_basis[0:1, :]
    npt.assert_allclose(actual, expected)

    actual = pod_basis(snapshots, rank_max=2, alternate_truncation_tol=0.6667)
    expected = full_basis[0:1, :]
    npt.assert_allclose(actual, expected)

    actual = pod_basis(snapshots, rank_max=1, alternate_truncation_tol=0.3334)
    expected = full_basis[0:1, :]
    npt.assert_allclose(actual, expected)


def test_discarded_energy_fractions():
    actual = discarded_energy_fractions([3.0, 2.0, 1.0])
    expected = np.array([1.0, np.sqrt(5.0/14.0), np.sqrt(1.0/14.0)])
    npt.assert_allclose(actual, expected)

    actual = discarded_energy_fractions(np.array([3.0, 2.0, 1.0]))
    expected = np.array([1.0, np.sqrt(5.0/14.0), np.sqrt(1.0/14.0)])
    npt.assert_allclose(actual, expected)

    actual = discarded_energy_fractions(np.array([1.0, 0.1, 0.01]))
    expected = np.array([1.0, 0.09999504987253119, 0.009949879346007117])
    npt.assert_allclose(actual, expected)
