import numpy as np
import numpy.linalg as la


def pod_basis(snapshots, rank_max=None, truncation_tol=None):
    """
    Generate a POD basis from a set of snapshots

    Parameters
    ----------
    snapshots : (n, N) array-like
        A multivector of n snapshots of size N.
    rank_max : int, optional
        An upper bound on the size of the POD basis. Truncation will be
        performed if necessary.  By default, do not perform truncation.
    truncation_tol : float, optional
        An energy-based criterion controlling truncation. Keep as many vectors
        in the POD basis such that the energy fraction of the discarded
        components is lower than truncation_tol. By default, do not perform
        truncation.

    Returns
    -------
    (m, N) ndarray
        An orthonormal multivector representing the POD basis.

    Notes
    -----
        A multivector is an (n, N)-array where:
          n denotes the number of vectors
          N denotes the number of entries in each vector
    """
    result, _ = pod_basis_and_singular_values(snapshots, rank_max, \
            truncation_tol)
    return result


def pod_basis_and_singular_values(snapshots, rank_max=None, truncation_tol=None):
    """
    Generate a POD basis and the associated singular values from a set of
    snapshots

    Parameters
    ----------
    snapshots : (n, N) array-like
        A multivector of n snapshots of size N.
    rank_max : int, optional
        An upper bound on the size of the POD basis. Truncation will be
        performed if necessary.  By default, do not perform truncation.
    truncation_tol : float, optional
        An energy-based criterion controlling truncation. Keep as many vectors
        in the POD basis such that the energy fraction of the discarded
        components is lower than truncation_tol. By default, do not perform
        truncation.

    Returns
    -------
    (m, N) ndarray
        An orthonormal multivector representing the POD basis.
    (m,) ndarray
        A 1-d vector of singular values

    Notes
    -----
        A multivector is an (n, N)-array where:
          n denotes the number of vectors
          N denotes the number of entries in each vector
    """
    snapshot_count = snapshots.shape[0]

    if rank_max is None:
        rank_max = snapshot_count

    v, s, u = la.svd(snapshots, full_matrices=False, compute_uv=True)

    if truncation_tol is not None:
        truncation_thresholds = discarded_energy_fractions(s)[::-1]
        drop_count = \
                min(np.searchsorted(truncation_thresholds, truncation_tol), \
                snapshot_count - 1)
        rank_max = min(snapshot_count - drop_count, rank_max)

    return (u, s) if rank_max >= snapshot_count else (u[:rank_max, :], s[:rank_max])


def discarded_energy_fractions(singular_values):
    """
    Compute the fraction of energy that would be discarded if truncation was
    performed right after the corresponding singular value.

    Parameters
    ----------
    singular_values : (n,) array-like
        Singular values, sorted in decreasing order

    Returns
    -------
    (n,) ndarray
        Energy fractions
    """
    modal_energies_squared = np.square(singular_values)

    reversed_modal_energies_squared = modal_energies_squared[::-1]
    reversed_cumulated_energies_squared = \
            np.cumsum(reversed_modal_energies_squared, out=reversed_modal_energies_squared)
    cumulated_energies_squared = reversed_cumulated_energies_squared[::-1]
    cumulated_energies = np.sqrt(cumulated_energies_squared, \
            out=cumulated_energies_squared)

    result = np.empty_like(singular_values)
    result[-1] = 0.0
    total_energy = cumulated_energies[0]
    np.divide(cumulated_energies[1:], total_energy, out=result[:-1])

    return result


def PODSnapshot(a, epsilon):
    eigenValues, eigenVectors = np.linalg.eigh(a,UPLO='L')

    idx = eigenValues.argsort()[::-1]
    eigenValues = eigenValues[idx]
    eigenVectors = eigenVectors[:,idx]

    id_max = 0
    bound = (epsilon**2)*eigenValues[0]
    for e in eigenValues:
      if e > bound:
        id_max += 1
    id_max2 = 0
    bound = (1-epsilon**2)*np.sum(eigenValues)
    temp = 0
    for e in eigenValues:
      temp += e
      if temp < bound:
        id_max2 += 1
    id_max = max(id_max, id_max2)

    return eigenValues[0:id_max], eigenVectors[:,0:id_max]


def CheckIntegrity():
    a = np.random.rand(10,10)
    a = np.dot(a.T,a)
    PODSnapshot(a, 1.e-6)
    pod_basis(a, None, 1.e-6)
    return 'ok'


if __name__ == '__main__':
    print(CheckIntegrity())# pragma: no cover
