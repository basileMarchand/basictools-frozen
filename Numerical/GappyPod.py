import numpy as np
import numpy.linalg as la
import scipy.linalg as sla
import math


# Generation of gappy masks via greedy algorithms
def deim(basis):
    """
    Generate empirical sample entries using the standard DEIM algorithm [1].

    Parameters
    ----------
    basis : array-like
        An orthonormal multivector.

    Returns
    -------
    ndarray of ints
        An 1-d array containing the sample entries.

    Notes
    -----
        A multivector is an (n, N)-array where n is the number of vectors and N
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
    ndarray of ints
        An 1-d array containing the sample entries.

    Notes
    -----
        A multivector is an (n, N)-array where n is the number of vectors and N
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


def greedy_gappy_mask(basis, \
        mask_size_max=None, \
        error_bound_max=None, \
        initial_mask=None, \
        acceleration=False):
    """
    Generate sample entries defining a gappy POD mask according to a greedy
    algorithm aiming at limiting both mask size and error due to the gappy
    approximation.

    Parameters
    ----------
    basis : (n, N) ndarray
        An orthonormal multivector.
    mask_size_max : int or None, optional
        An upper bound on the size of the mask. None implies that the option is
        inactive. Default is None.
    error_bound_max : scalar or None, optional
        Specifies the threshold below which the gappy aliasing error bound
        associated with the mask should be. None implies that the option is
        inactive. Default is None.
    initial_mask : (m,) ndarray of ints or None, optional
        Makes the first m mask entries match those in initial_mask. None lets
        the algorithm determine all mask entries. Default is None.
    acceleration: bool, optional
        Replaces the exact computation of the spectrum by a faster
        approximation based on the secular equation adapted from [1]. Default
        is False.

    Returns
    -------
    (k,) ndarray of ints
        The sample entries defining the mask, ranked according to their
        selection order.

    Notes
    -----
        A multivector is an (n, N)-array where n is the number of vectors and N
        the numbers of dofs.

    References
    ----------
        [1] R. Zimmermann & K. Willcox, "An Accelerated Greedy Missing Point
        Estimation Procedure", SIAM Journal on Scientific Computing, vol.38,
        num. 5, pp 2827--2850, 2016.
    """
    sg = ExactSpectrumGenerator(basis) if acceleration is False \
            else AcceleratedSpectrumGenerator(basis)
    ec = eigenvalue_harmonic_mean
    sc = ErrorBoundCriterion(error_bound_max)
    return greedy_mpe_mask(sg, ec, sc, \
        sample_size_max=mask_size_max, \
        initial_sample=initial_mask)


# Error bounds and estimators
def aliasing_error_bound(basis, sampling):
    """
    Compute a tight bound for the aliasing error (in the Euclidian 2-norm) for
    the gappy projection compared to the complete orthognal projection.

    Parameters
    ----------
    basis : (n, N) ndarray
        An orthonormal multivector.
    sampling: 1-d ndarray of ints
        The indices of the sample entries defining the mask.

    Returns
    -------
    float
        The aliasing error bound.

    Notes
    -----
        A multivector is an (n, N)-array where n is the number of vectors and N
        the numbers of dofs.
    """
    sampled_basis = basis[:, sampling]
    sigma_min = la.norm(sampled_basis, -2)
    result_squared = (1.0 / sigma_min)**2 - 1.0
    # Numerical errors could lead to taking the sqrt of a negative value
    result_squared_fixed = max(result_squared, 0.0)
    return math.sqrt(result_squared_fixed)


def expected_aliasing_error_bound(basis, sampling):
    bound = aliasing_error_bound(basis, sampling)
    ratio = expected_value_norm_bound_ratio(basis.shape[-1], basis.shape[0])
    return bound * ratio


def root_mean_square_aliasing_error(basis, sampling):
    sampled_basis = basis[:, sampling]
    singular_vals = sla.svdvals(sampled_basis)
    spectrum = np.square(singular_vals, out=singular_vals)
    spectrum_reciprocal = np.reciprocal(spectrum, out=spectrum)
    frobenius_norm_squared = np.sum(spectrum_reciprocal) - spectrum.shape[0]
    # Numerical errors could lead to taking the sqrt of a negative value
    frobenius_norm_squared_fixed = max(frobenius_norm_squared, 0.0)
    return math.sqrt(frobenius_norm_squared_fixed / basis.shape[-1])


def expected_value_norm_bound_ratio(dim, rank):
    import scipy.special
    return scipy.special.beta(0.5 * dim, 0.5) / scipy.special.beta(0.5 * rank, 0.5)


# Empirical/gappy operator assembly
def empirical_operator(basis, sample_entries):
    sampled_basis = basis[:, sample_entries]
    nonsingular = sampled_basis.shape[0] == sampled_basis.shape[1]
    matrixT = la.inv(sampled_basis) if nonsingular else la.pinv(sampled_basis)
    return lambda x: np.dot(x, matrixT)


# Greedy gappy mask: Implementation details
## Common driver implementation
def greedy_mpe_mask( \
        spectrum_generator, enrichment_criterion, stopping_criterion, \
        sample_size_max=None, initial_sample=None):
    greedy_generator = GreedyGappyMaskGenerator( \
            spectrum_generator, enrichment_criterion, initial_sample)

    if sample_size_max is None:
        sample_size_max = greedy_generator.mask_size_max()
    else:
        sample_size_max = min(sample_size_max, \
                greedy_generator.mask_size_max())

    while greedy_generator.mask_size() < sample_size_max \
            and not stopping_criterion(greedy_generator.spectrum()):
        greedy_generator.mask_size_inc()

    if sample_size_max < greedy_generator.mask_size():
        return greedy_generator.mask()[:sample_size_max]
    else:
        return greedy_generator.mask()


## Greedy enrichment criteria
def lowest_eigenvalue(spectra, current_spectrum=None, iteration_index=None):
    return spectra[:, 0]


def lowest_eigenvalue_with_periodic_switch(period=None):
    if period is None:
        period = 3
    def criterion(spectra, current_spectrum=None, iteration_index=None):
        switching = np.remainder(iteration_index + 1, period) == 0
        return spectra[:, int(switching)]
    return criterion


def eigenvalue_harmonic_mean(spectra, current_spectrum=None, iteration_index=None):
    return harmonic_mean(spectra)


def lowest_eigenvalue_with_potential(threshold=None):
    if threshold is None:
        threshold = 0.05
    index = lowest_eigenvalue_index_with_potential(threshold)
    def criterion(spectra, current_spectrum, iteration_index=None):
        return spectra[:, index(current_spectrum)]
    return criterion


def lowest_eigenvalue_with_potential_original(threshold=None):
    if threshold is None:
        threshold = 0.05
    index = lowest_eigenvalue_index_with_potential_original(threshold)
    def criterion(spectra, current_spectrum, iteration_index=None):
        return spectra[:, index(current_spectrum)]
    return criterion


### Frobenius criterion
def harmonic_mean(a, axis=-1):
    masked_a = np.ma.masked_less_equal(a, 0.0)
    np.ma.set_fill_value(masked_a, 1.0)
    filled_a = np.ma.filled(masked_a)
    result = np.mean(np.reciprocal(filled_a, out=filled_a), axis=-1)
    np.reciprocal(result, out=result)
    if np.ma.is_masked(masked_a):
        result_mask = np.any(np.ma.getmask(masked_a), axis=-1)
        masked_result = np.ma.array(result, mask=result_mask)
        np.ma.set_fill_value(masked_result, 0.0)
        result = np.ma.filled(masked_result)
    return result


### Criteria based on current spectrum
def lowest_eigenvalue_index_with_potential(threshold):
    def rank(spectrum):
        potential = np.empty_like(spectrum)
        # General case: All but last entries
        potential_view = potential[..., :-1]
        # Goal: np.diff(spectrum, axis=-1, out=potential_view)
        # But np.diff does not accept an out parameter
        np.subtract(spectrum[..., 1:], spectrum[..., :-1], out=potential_view)
        potential_view /= spectrum[..., 1:]
        # Special case: last entry
        np.subtract(1.0, spectrum[..., -1], out=potential[..., -1])
        invalid = potential < threshold
        idx = first_false_index(invalid, axis=-1, tmp=invalid)
        return idx % potential.shape[-1]
    return rank


def lowest_eigenvalue_index_with_potential_original(threshold):
    def rank(spectrum):
        if spectrum.shape[-1] <= 2:
            return 0
        potential = np.diff(spectrum, axis=-1)
        potential /= spectrum[..., 1:]
        invalid = potential <= threshold
        idx = first_false_index(invalid, axis=-1, tmp=invalid)
        return idx % potential.shape[-1]
    return rank


def first_false_index(a, axis=-1, tmp=None):
    return np.sum(np.logical_and.accumulate(a, axis=axis, out=tmp), axis=axis)


## Stopping criterion
class ErrorBoundCriterion(object):
    def __init__(self, target=None):
        if target is None:
            target = 0.0
        self.lambda_min_target = 1.0 / (target**2 + 1.0)

    def __call__(self, spectrum):
        if spectrum.size == 0:
            return True
        lambda_min = spectrum[0]
        return lambda_min >= self.lambda_min_target


## Spectrum generators
class AbstractSpectrumGenerator(object):
    def basis_rank(self):
        pass

    def atom_count(self):
        pass

    def current_spectrum(self):
        pass

    def selection_spectra(self):
        pass

    def atom_add(self, index):
        pass

    def mask_add(self, indices):
        pass


class ExactSpectrumGenerator(object):
    def __init__(self, basis):
        assert basis.ndim == 2
        # TODO: Handle case of higher dimensional atoms
        self.contributions = \
            np.expand_dims(basis.T, axis=-2) * \
            np.expand_dims(basis.T, axis=-1)
        # Note: The above statement is slightly faster than the alternative:
        # entry_contributions = np.einsum('ji,ki->ijk', basis, basis)
        atom_count = self.contributions.shape[0]
        basis_rank = self.contributions.shape[1]
        self.correlation = np.zeros((basis_rank, basis_rank), dtype=basis.dtype)
        self.spectrum = np.zeros((basis_rank,), dtype=basis.dtype)
        self.candidate_correlations = np.empty_like(self.contributions)
        self.candidate_spectra = np.empty((atom_count, basis_rank), dtype=basis.dtype)
        self.is_candidate_data_uptodate = False

    def basis_rank(self):
        return self.contributions.shape[1]

    def atom_count(self):
        return self.contributions.shape[0]

    def current_spectrum(self):
        return self.spectrum

    def selection_spectra(self):
        if not self.is_candidate_data_uptodate:
            self._update_candidate_data()
        return self.candidate_spectra

    def atom_add(self, index):
        if self.is_candidate_data_uptodate:
            np.copyto(self.correlation, self.candidate_correlations[index])
            np.copyto(self.spectrum, self.candidate_spectra[index])
        else:
            self.correlation += self.contributions[index]
            self.spectrum = la.eigvalsh(self.correlation)
            self.spectrum.sort()
        self.contributions[index].fill(0.0)
        self.is_candidate_data_uptodate = False

    def mask_add(self, indices):
        if indices.size == 0:
            return
        self.correlation += np.sum(self.contributions[indices], axis=0)
        self.spectrum = la.eigvalsh(self.correlation)
        self.spectrum.sort()
        self.contributions[indices] = 0.0
        self.is_candidate_data_uptodate = False

    def _update_candidate_data(self):
        np.add(self.contributions, self.correlation, \
                out=self.candidate_correlations)
        for a, c in enumerate(self.candidate_correlations):
            self.candidate_spectra[a, :] = la.eigvalsh(c)
        self.candidate_spectra.sort()
        self.is_candidate_data_uptodate = True


class AcceleratedSpectrumGenerator(object):
    def __init__(self, basis):
        assert basis.ndim == 2
        atom_count = basis.shape[1]
        basis_rank = basis.shape[0]
        self.basis = np.array(basis)
        self.spectrum = np.zeros((basis_rank,), dtype=basis.dtype)
        self.mask_data = np.empty((atom_count,), dtype=np.int_)
        #TODO mask_data is a duplicated variable (cf algorithm driver function)
        self.mask = self.mask_data[:0]
        self.candidate_spectra = np.empty((atom_count, basis_rank), dtype=basis.dtype)
        self.is_candidate_data_uptodate = False

    def basis_rank(self):
        return self.basis.shape[0]

    def atom_count(self):
        return self.basis.shape[1]

    def current_spectrum(self):
        return self.spectrum

    def selection_spectra(self):
        if not self.is_candidate_data_uptodate:
            self._update_candidate_data()
        return self.candidate_spectra

    def atom_add(self, index):
        s = self.mask.shape[0]
        self.mask_data[s] = index
        self.mask = self.mask_data[:s+1]
        self._update_spectrum()
        self.is_candidate_data_uptodate = False

    def mask_add(self, indices):
        if indices.size == 0:
            return
        s = self.mask.shape[0]
        c = indices.shape[0]
        np.copyto(self.mask_data[s:s+c], indices)
        self.mask = self.mask_data[:s+c]
        self._update_spectrum()
        self.is_candidate_data_uptodate = False

    def _update_spectrum(self):
        sampled_basis = self.basis[:, self.mask]
        singular_vals = sla.svdvals(sampled_basis)
        np.square(singular_vals, out=singular_vals)
        np.copyto( \
                self.spectrum[-singular_vals.shape[0]:], \
                singular_vals[::-1])

    def _update_candidate_data(self):
        updated_spectrum_lower_bound(self.basis, self.mask, \
                out=self.candidate_spectra[:, ::-1])


## Underlying stateful object used by drivers
class GreedyGappyMaskGenerator(object):
    def __init__(self, spectrum_generator, enrichment_criterion, initial_mask=None):
        self.spectrum_generator = spectrum_generator
        self.enrichment_criterion = enrichment_criterion
        mask_size_max = self.spectrum_generator.atom_count()
        basis_rank = self.spectrum_generator.basis_rank()
        self.mask_data = np.empty((mask_size_max,), dtype=np.int_)
        if initial_mask is None:
            self.current_mask = self.mask_data[:0]
            self.current_spectrum = np.zeros((basis_rank,), dtype=np.double)
        else:
            assert initial_mask.shape[0] <= mask_size_max
            initial_mask_size = initial_mask.shape[0]
            self.current_mask = self.mask_data[:initial_mask_size]
            np.copyto(self.current_mask, initial_mask)
            self.spectrum_generator.mask_add(self.current_mask)
            self.current_spectrum = np.empty((basis_rank,), dtype=np.double)
            np.copyto(self.current_spectrum, \
                    self.spectrum_generator.current_spectrum())
        self.enrichment_iteration_index = 0

    def mask(self):
        return self.current_mask

    def mask_size(self):
        return self.current_mask.shape[0]

    def mask_size_max(self):
        return self.mask_data.shape[0]

    def spectrum(self):
        return self.current_spectrum

    def mask_size_inc(self):
        if self.mask_size() == self.mask_size_max():
            raise Exception('Reached maximum mask size')
        k = self.mask_size()
        basis_rank = self.spectrum_generator.basis_rank()
        rank_deficiency_step = 1 # TODO
        current_rank_deficiency = max(basis_rank - (rank_deficiency_step * k), 0)
        eigenvalue_index_min = max(current_rank_deficiency - rank_deficiency_step, 0)

        spectra = self.spectrum_generator.selection_spectra()
        fitnesses = self.enrichment_criterion( \
                spectra[:, eigenvalue_index_min:], \
                current_spectrum=self.current_spectrum[eigenvalue_index_min:], \
                iteration_index=self.enrichment_iteration_index)
        selected = np.argmax(fitnesses)

        self.spectrum_generator.atom_add(selected)
        self.mask_data[k] = selected
        self.current_mask = self.mask_data[:k+1]
        np.copyto(self.current_spectrum, \
                self.spectrum_generator.current_spectrum())

        has_enrichment_started = current_rank_deficiency == 0
        if has_enrichment_started:
            self.enrichment_iteration_index += 1


## Approximate spectrum computed via the secular equation
def updated_spectrum_lower_bound(basis, sampling, out=None):
    # basis -> (n, N), sampling -> (s,), with s >= n
    if out is None:
        result = np.empty_like(basis.T)
    else:
        assert out.shape == basis.T.shape and out.dtype == basis.dtype
        result = out
    # result -> (N, n)

    if sampling.shape[0] == 0:
        result.fill(0.0)
        vsq = np.square(basis.T)
        result[:, 0] = np.sum(vsq, axis=1)
        return result

    sampled_basis = basis[:, sampling]
    # sampled_basis -> (n, s)
    deficient = sampled_basis.shape[1] < sampled_basis.shape[0]
    phi, sigma, __ = \
            la.svd(sampled_basis, full_matrices=deficient, compute_uv=True)
    # phi -> (n, n), sigma -> (k,), __ -> (k, s) with k = min(n, s)
    assert phi.shape[0] == phi.shape[1] == basis.shape[0]
    d = np.zeros((basis.shape[0],), dtype=basis.dtype)
    d[:sigma.shape[0]] = np.square(sigma)
    # d -> (n,)
    v = np.dot(basis.T, phi)
    # v -> (N, n)
    vsq = np.square(v, out=v)
    # vsq -> (N, n)

    super_deficient = sampled_basis.shape[1] + 1 < sampled_basis.shape[0]
    if super_deficient:
        # Updated rank will be at most r = (s+1) < n
        s = sampled_basis.shape[1]
        vsq[:, s] = np.sum(vsq[:, s:], axis=1)
        d = d[:s+1]
        vsq = vsq[:, :s+1]
        result.fill(0.0)
        res = result[:, :s+1]
    else:
        res = result
    # d -> (r,), vsq -> (N, r), res -> (N, r) with r = min(n, s+1)

    dd = np.expand_dims(d, axis=0)
    # dd -> (1, r)
    d_diffs = np.ma.array( \
            dd.T - dd, \
            mask=np.eye(d.shape[0]) + np.eye(d.shape[0], k=-1), \
            fill_value=1.0)
    # d_diffs -> (r[poles], r[bounds])
    cup_vec = np.ma.divide(1.0, d_diffs)
    cup_vec = np.ma.filled(cup_vec, fill_value=0.0)
    # cup_vec -> (r[poles], r[bounds])
    cup_array = np.dot(vsq, cup_vec[:, :-1])
    cup_array += 1.0
    # cup_array -> (N, r-1)
    cup_sign = np.sign(cup_array)
    # cup_sign -> (N, r-1)
    cup_array_inv = np.reciprocal(cup_array, out=cup_array)
    # cup_array_inv -> (N, r-1)
    g1 = dd[:, :-1] + dd[:, 1:]
    g2 = dd[:, :-1] * dd[:, 1:]
    # g1, g2 -> (1, r-1)
    alpha1 = (vsq[:, 1:] + vsq[:, :-1]) * cup_array_inv + g1
    alpha2 = (dd[:, :-1] * vsq[:, 1:] + dd[:, 1:] * vsq[:, :-1]) * cup_array_inv + g2
    # alpha1, alpha2 -> (N, r-1)
    sqrt_arg = 0.25 * np.square(alpha1) - alpha2
    # sqrt_arg -> (N, r-1)
    np.maximum(sqrt_arg, 0.0, out=sqrt_arg)
    res[:, 1:] = 0.5 * alpha1 - cup_sign * np.sqrt(sqrt_arg)

    v_norm_sq = np.sum(vsq, axis=1)
    # v_norm_sq -> (N, )
    bound0 = np.add(d[0], v_norm_sq, out=v_norm_sq)
    # bound0 -> (N, )
    d_diffs0 = dd[:, 1:] - np.expand_dims(bound0, axis=1)
    # d_diffs0 -> (N, r-1)
    cup_vec0 = np.reciprocal(d_diffs0, out=d_diffs0)
    # cup_vec0 -> (N, r-1)
    cup_vec_fracs0 = np.multiply(cup_vec0, vsq[:, 1:], out=cup_vec0)
    # cup_vec_fracs0 -> (N, r-1)
    cup_array0 = np.sum(cup_vec_fracs0, axis=1, out=cup_array[:, 0])
    # cup_array0 -> (N,)
    cup_array0 += 1.0
    np.divide(vsq[:, 0], cup_array0, out=result[:, 0])
    res[:, 0] += d[0]

    res[sampling, :] = dd
    return result


# Efficient greedy step (Zimmermann/Willcox)
def efficient_greedy_step(basis, sampling):
    lambda_cup = updated_lowest_eigenvalue_lower_bound(basis, sampling)
    return np.argmax(lambda_cup)

def updated_lowest_eigenvalue_lower_bound(basis, sampling):
    # basis -> (n, N), sampling -> (s,), with s >= n
    # basis[:, sampling] -> (n, s), with s >= n
    phi, sigma, __ = la.svd(basis[:, sampling] , full_matrices=False, compute_uv=True)
    # phi -> (n, r), sigma -> (r,), __ -> (r, s) where r = min(n, s) = n
    d = np.square(sigma, out=sigma)
    # d -> (r,) where r = min(n, s) = n
    v = np.dot(phi.T, basis)
    # v -> (n, N)
    vsq = np.square(v, out=v)
    # vsq -> (n, N)
    cup_vec = 1.0 / (d[:-2] - d[-2])
    # cup_vec -> (r-2,) where r = min(n, s) = n
    cup_array = np.dot(cup_vec, vsq[:-2]) + 1.0
    # cup_array -> (N,)
    cup_array_inv = np.reciprocal(cup_array, out=cup_array)
    # cup_array_inv -> (N,)
    g1 = d[-2] + d[-1]
    g2 = d[-2] * d[-1]
    # g1, g2 -> (1,)
    alpha1 = (vsq[-1] + vsq[-2]) * cup_array_inv + g1
    alpha2 = (d[-2] * vsq[-1] + d[-1] * vsq[-2]) * cup_array_inv + g2
    # alpha1, alpha2 -> (N,)
    sqrt_arg = 0.25 * np.square(alpha1) - alpha2
    np.maximum(sqrt_arg, 0.0, out=sqrt_arg)
    result = 0.5 * alpha1 - np.sqrt(sqrt_arg)
    result[sampling] = d[-1]
    # result -> (N,)
    return result


# Testing support
def updated_lowest_eigenvalue(basis, sampling):
    contributions = \
            np.expand_dims(basis.T, axis=-2) * np.expand_dims(basis.T, axis=-1)
    base_correlation = np.sum(contributions[sampling], axis=0)
    contributions[sampling] = 0.0
    correlations  = np.add(contributions, base_correlation, out=contributions)
    result = np.empty(correlations.shape[0])
    for i, c in enumerate(correlations):
        result[i] = la.norm(c, ord=-2)
    return result

def updated_spectrum(basis, sampling):
    contributions = \
            np.expand_dims(basis.T, axis=-2) * np.expand_dims(basis.T, axis=-1)
    base_correlation = np.sum(contributions[sampling], axis=0)
    contributions[sampling] = 0.0
    correlations = np.add(contributions, base_correlation, out=contributions)
    result = np.empty_like(basis.T)
    for i, c in enumerate(correlations):
        result[i, :] = la.eigvalsh(c)
    result.sort()
    return result[:, ::-1]


# Didactic algorithm equivalent to QDEIM
def fake_q_deim(basis):
    basis_rank = basis.shape[0]
    sample_entries = np.empty(basis_rank, dtype=np.int_)
    sampled_basis = np.empty((basis_rank, basis_rank), dtype=basis.dtype)

    for k in xrange(0, basis_rank):
        if k > 0:
            phi, __1, __2 = \
                la.svd(sampled_basis[:, :k], full_matrices=True, compute_uv=True)
            v = np.dot(basis.T, phi[:, k:])
        else:
            v = np.array(basis.T)
        vsq = np.square(v, out=v)
        fitnesses = np.sum(vsq, axis=1)
        sample_entries[k] = np.argmax(fitnesses)
        np.copyto(sampled_basis[:, k], basis[:, sample_entries[k]], casting='no')

    return sample_entries
