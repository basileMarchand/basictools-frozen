import numpy as np
import numpy.linalg as la
import scipy.linalg as sla


def deim(basis):
    """
    Generate empirical sample entries using the standard DEIM algorithm [1].

    Parameters
    ----------
    basis : array-like
        An orthonormal multivector.

    Returns
    -------
    ndarray
        An 1-d array containing the sample entries.

    Notes
    -----
        A multivector is (n, N)-array where n is the number of vectors and N
        the numbers of dofs.

    References
    ----------
        [1] S. Chaturantabut & D. Sorensen, "Nonlinear Model Reduction via
        Discrete Empirical Interpolation", SIAM Journal on Scientific
        Computing, vol. 32, num. 5, pp. 2737-2764, 2010.
    """
    basis_rank = basis.shape[0]
    sample_entries = np.empty(basis_rank, dtype=np.int_)

    if basis_rank > 0:
        sample_entries[0] = np.argmax(np.fabs(basis[0, :]))

        sampled_basis = np.empty((basis_rank, basis_rank - 1), dtype=np.double)

    for k in xrange(1, basis_rank):
        added_basis_slice = sampled_basis[:, k - 1]
        np.copyto(added_basis_slice, basis[:, sample_entries[k - 1]], casting='no')

        s = slice(0, k)
        interpolation_matrix = np.transpose(sampled_basis[s, s])
        rhs = sampled_basis[k, s]
        components = la.solve(interpolation_matrix, rhs)
        np.negative(components, out=components)
        residual = np.dot(components, basis[s, :])
        residual += basis[k, :]

        sample_entries[k] = np.argmax(np.fabs(residual))

    return sample_entries


def q_deim(basis):
    """
    Generate empirical sample entries using the Q-DEIM algorithm [1].

    Parameters
    ----------
    basis : array-like
        An orthonormal multivector.

    Returns
    -------
    ndarray
        An 1-d array containing the sample entries.

    Notes
    -----
        A multivector is (n, N)-array where n is the number of vectors and N
        the numbers of dofs.

    References
    ----------
        [1] Z. Drmac and S. Gugercin. A New Selection Operator for the Discrete
        Empirical Interpolation Method : Improved a priori error bound and
        extensions, 2015. URL: http://arxiv.org/abs/1505.00370.
    """
    basis_rank = basis.shape[0]

    if basis_rank == 0:
        return np.empty(basis_rank, dtype=np.int_)

    Q, R, P = sla.qr(basis, pivoting=True)
    return P[:basis_rank]
