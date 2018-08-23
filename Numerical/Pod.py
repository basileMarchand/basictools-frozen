import numpy as np
import numpy.linalg as la


def pod_basis(snapshots, \
        rank_max=None, truncation_tol=None, alternate_truncation_tol=None):
    """
    Generate a POD basis from a set of snapshots

    Parameters
    ----------
    snapshots : (n, N) array-like
        A multivector of n snapshots of size N.
    rank_max : int, optional
        An upper bound on the size of the POD basis. Truncation will be
        performed if necessary. By default, do not perform truncation.
    truncation_tol : float, optional
        An energy-based criterion controlling truncation. Keep as few vectors
        in the POD basis as necessary such that the cumulated energy fraction
        of the discarded components is lower than or equal to truncation_tol.
        By default, do not perform truncation.
    alternate_truncation_tol: float, optional
        An energy-based criterion controlling truncation. Discard from the POD
        basis any vector such that the ratio of the associated energy to the
        energy of the leading vector is lower than or equal to
        alternate_truncation_tol. By default, do not perform truncation.

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
    result, __ = pod_basis_and_singular_values(snapshots, \
            rank_max, truncation_tol, alternate_truncation_tol)
    return result


def pod_basis_and_singular_values(snapshots, \
        rank_max=None, truncation_tol=None, alternate_truncation_tol=None):
    """
    Generate a POD basis and the associated singular values from a set of
    snapshots

    Parameters
    ----------
    snapshots : (n, N) array-like
        A multivector of n snapshots of size N.
    rank_max : int, optional
        An upper bound on the size of the POD basis. Truncation will be
        performed if necessary. By default, do not perform truncation.
    truncation_tol : float, optional
        An energy-based criterion controlling truncation. Keep as few vectors
        in the POD basis as necessary such that the cumulated energy fraction
        of the discarded components is lower than or equal to truncation_tol.
        By default, do not perform truncation.
    alternate_truncation_tol: float, optional
        An energy-based criterion controlling truncation. Discard from the POD
        basis any vector such that the ratio of the associated energy to the
        energy of the leading vector is lower than or equal to
        alternate_truncation_tol. By default, do not perform truncation.

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
    __, s, u = la.svd(snapshots, full_matrices=False, compute_uv=True)
    rank = truncated_pod_basis_rank(s, \
            rank_max, truncation_tol, alternate_truncation_tol)
    return (u[:rank, :], s[:rank])


def truncated_pod_basis_rank(singular_values, \
        rank_max=None, truncation_tol=None, alternate_truncation_tol=None):
    snapshot_count = singular_values.size
    result = min(snapshot_count, max(rank_max, 0)) \
            if rank_max is not None else snapshot_count

    """
    if truncation_tol is not None:
        truncation_thresholds = discarded_energy_fractions(s)[::-1]
        drop_count = np.searchsorted( \
                truncation_thresholds, truncation_tol, side='right')
        rank_max = min(min(snapshot_count,s.shape[0]) - drop_count, rank_max)

    if alternate_truncation_tol is not None and snapshot_count > 0:
        leading_vector_energy = s[0]
        cutoff = leading_vector_energy * alternate_truncation_tol
        drop_count = np.searchsorted(s[::-1], cutoff, side='right')
        rank_max = min(min(snapshot_count,s.shape[0]) - drop_count, rank_max)

    #return (u, s) if rank_max >= snapshot_count else (u[:rank_max, :], s[:rank_max])
    return (u, s) if rank_max >= snapshot_count else (u[:rank_max, :], s)"""

    if truncation_tol is not None and result > 0:
        fractions = discarded_energy_fractions(singular_values)
        result -= truncation_level(fractions[:result], truncation_tol)

    if alternate_truncation_tol is not None and result > 0:
        leading_vector_energy = singular_values[0]
        threshold = leading_vector_energy * alternate_truncation_tol
        result -= truncation_level(singular_values[:result], threshold)

    return result



def discarded_energy_fractions(singular_values):
    """
    Compute the fraction of energy that would be discarded if truncation was
    performed just before the corresponding singular value.

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
    result[0] = 1.0
    total_energy = cumulated_energies[0]
    np.divide(cumulated_energies[1:], total_energy, out=result[1:])

    return result


def truncation_level(values, threshold):
    return np.searchsorted(values[::-1], threshold, side='right')


def PODSnapshot(a, epsilon, nLinearSize = None):
    
    if nLinearSize is not None:#convert a vector into a lower triangular matrix
      ac = np.copy(a)
      a = np.zeros((nLinearSize,nLinearSize))
      a[np.tril_indices(nLinearSize)] = ac

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

    return eigenValues, eigenValues[0:id_max], eigenVectors[:,0:id_max]


def CheckIntegrity():
    import BasicTools.Numerical.test.test_pod as tests
    test_functions = \
            (t for n, t in tests.__dict__.items() if n.startswith("test"))
    for f in test_functions:
        f()

    a = np.random.rand(10,10)
    a = np.dot(a.T,a)
    PODSnapshot(a, 1.e-6)
    pod_basis(a, None, 1.e-6)

    return 'ok'


if __name__ == '__main__':
    print(CheckIntegrity()) # pragma: no cover
