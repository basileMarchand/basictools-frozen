from BasicTools.Numerical.GappyPod import *
import BasicTools.Numerical.Pod as Pod

import numpy as np
import numpy.testing as npt


def assert_sample(algorithm, basis, expected_entries):
    actual = algorithm(basis); actual.sort()
    expected = np.array(expected_entries); expected.sort()
    npt.assert_array_equal(actual, expected)


def empty():
    return np.empty((0, 10), dtype=np.double)


def negative_identity():
    return -np.eye(3, 10, dtype=np.double)


def shifted_identity():
    return np.eye(3, 10, k=1, dtype=np.double)


def staggered():
    return (negative_identity() + (np.sqrt(2.0) * shifted_identity())) \
            / np.sqrt(3.0)


def test_deim():
    algorithm = deim

    assert_sample(algorithm, empty(), [])
    assert_sample(algorithm, negative_identity(), [0, 1, 2])
    assert_sample(algorithm, shifted_identity(), [1, 2, 3])
    assert_sample(algorithm, staggered(), [1, 2, 3])

def test_q_deim():
    algorithm = q_deim

    assert_sample(algorithm, empty(), [])
    assert_sample(algorithm, negative_identity(), [0, 1, 2])
    assert_sample(algorithm, shifted_identity(), [1, 2, 3])
    assert_sample(algorithm, staggered(), [1, 2, 3])

def test_fake_q_deim():
    algorithm = fake_q_deim

    assert_sample(algorithm, empty(), [])
    assert_sample(algorithm, negative_identity(), [0, 1, 2])
    assert_sample(algorithm, shifted_identity(), [1, 2, 3])
    assert_sample(algorithm, staggered(), [1, 2, 3])


greedy_mpe_mask_driver_enrichments = {
        'frobenius': eigenvalue_harmonic_mean,
        'greediest': lowest_eigenvalue,
        'potential': lowest_eigenvalue_with_potential(),
        'potential_original': lowest_eigenvalue_with_potential_original(),
        'switching' : lowest_eigenvalue_with_periodic_switch()
        }


greedy_mpe_mask_driver_evaluators = { \
        'exact': ExactSpectrumGenerator, \
        'accelerated': AcceleratedSpectrumGenerator}

def make_greedy_mpe_mask_drivers():
    def define(enrichment, enrichment_name, evaluator, evaluator_name):
        return lambda \
                basis, \
                sample_size_max=None, \
                error_bound_max=None, \
                initial_sample=None: \
            greedy_mpe_mask( \
                    evaluator(basis), \
                    enrichment, \
                    ErrorBoundCriterion(error_bound_max), \
                    sample_size_max=sample_size_max, \
                    initial_sample=initial_sample)
    es = greedy_mpe_mask_driver_enrichments.items()
    ss = greedy_mpe_mask_driver_evaluators.items()
    return {en: {sn: define(e, en, s, sn) for sn, s in ss} for en, e in es}

greedy_mpe_mask_drivers = make_greedy_mpe_mask_drivers()

frobenius_exact_greedy_mpe_mask = \
        greedy_mpe_mask_drivers['frobenius']['exact']

greediest_exact_greedy_mpe_mask = \
        greedy_mpe_mask_drivers['greediest']['exact']

potential_exact_greedy_mpe_mask = \
        greedy_mpe_mask_drivers['potential']['exact']

potential_original_exact_greedy_mpe_mask = \
        greedy_mpe_mask_drivers['potential_original']['exact']

switching_exact_greedy_mpe_mask = \
        greedy_mpe_mask_drivers['switching']['exact']

frobenius_accelerated_greedy_mpe_mask = \
        greedy_mpe_mask_drivers['frobenius']['accelerated']

greediest_accelerated_greedy_mpe_mask = \
        greedy_mpe_mask_drivers['greediest']['accelerated']

potential_accelerated_greedy_mpe_mask = \
        greedy_mpe_mask_drivers['potential']['accelerated']

potential_original_accelerated_greedy_mpe_mask = \
        greedy_mpe_mask_drivers['potential_original']['accelerated']

switching_accelerated_greedy_mpe_mask = \
        greedy_mpe_mask_drivers['switching']['accelerated']

def test_frobenius_exact_greedy_mpe_mask():
    def algorithm(sample_size_max=None):
        def func(b):
            return frobenius_exact_greedy_mpe_mask(
                    b, 
                    sample_size_max=sample_size_max
                    )
        return func

    assert_sample(algorithm(), empty(), [])
    assert_sample(algorithm(), negative_identity(), [0, 1, 2])
    assert_sample(algorithm(3), shifted_identity(), [1, 2, 3])
    assert_sample(algorithm(2), staggered(), [1, 3])
    assert_sample(algorithm(3), staggered(), [1, 2, 3])
    assert_sample(algorithm(4), staggered(), [0, 1, 2, 3])

def test_greediest_exact_greedy_mpe_mask():
    def algorithm(sample_size_max=None):
        def func(b):
            return greediest_exact_greedy_mpe_mask(
                    b,
                    sample_size_max=sample_size_max
                    )
        return func

    assert_sample(algorithm(), empty(), [])
    assert_sample(algorithm(), negative_identity(), [0, 1, 2])
    assert_sample(algorithm(3), shifted_identity(), [1, 2, 3])
    assert_sample(algorithm(2), staggered(), [1, 3])
    assert_sample(algorithm(3), staggered(), [1, 2, 3])
    assert_sample(algorithm(4), staggered(), [0, 1, 2, 3])

def test_frobenius_accelerated_greedy_mpe_mask():
    def algorithm(sample_size_max=None):
        def func(b):
            return frobenius_accelerated_greedy_mpe_mask(
                    b,
                    sample_size_max=sample_size_max
                    )
        return func

    assert_sample(algorithm(), empty(), [])
    #Following tests deactivated for lack of support for repeated eigenvalues
    #assert_sample(algorithm(), negative_identity(), [0, 1, 2])
    #assert_sample(algorithm(3), shifted_identity(), [1, 2, 3])
    assert_sample(algorithm(2), staggered(), [1, 3])
    assert_sample(algorithm(3), staggered(), [1, 2, 3])
    assert_sample(algorithm(4), staggered(), [0, 1, 2, 3])

def test_frobenius_accelerated_greedy_mpe_mask_with_qdeim():
    def algorithm(sample_size_max=None):
        def func(b):
            return frobenius_accelerated_greedy_mpe_mask(b, \
                    sample_size_max=sample_size_max, \
                    initial_sample=q_deim(b))
        return func

    assert_sample(algorithm(), empty(), [])
    assert_sample(algorithm(), negative_identity(), [0, 1, 2])
    assert_sample(algorithm(3), shifted_identity(), [1, 2, 3])
    assert_sample(algorithm(2), staggered(), [1, 2]) # truncated q_deim
    assert_sample(algorithm(3), staggered(), [1, 2, 3])
    assert_sample(algorithm(4), staggered(), [0, 1, 2, 3])

def test_greediest_accelerated_greedy_mpe_mask():
    def algorithm(sample_size_max=None):
        def func(b):
            return greediest_accelerated_greedy_mpe_mask(b, \
                    sample_size_max=sample_size_max)
        return func

    assert_sample(algorithm(), empty(), [])
    #Following tests deactivated for lack of support for repeated eigenvalues
    #assert_sample(algorithm(), negative_identity(), [0, 1, 2])
    #assert_sample(algorithm(3), shifted_identity(), [1, 2, 3])
    assert_sample(algorithm(2), staggered(), [1, 3])
    assert_sample(algorithm(3), staggered(), [1, 2, 3])
    assert_sample(algorithm(4), staggered(), [0, 1, 2, 3])

def test_greediest_accelerated_greedy_mpe_mask_with_qdeim():
    def algorithm(sample_size_max=None):
        def func(b):
            return greediest_accelerated_greedy_mpe_mask(b, \
                    sample_size_max=sample_size_max, \
                    initial_sample=q_deim(b))
        return func

    assert_sample(algorithm(), empty(), [])
    assert_sample(algorithm(), negative_identity(), [0, 1, 2])
    assert_sample(algorithm(3), shifted_identity(), [1, 2, 3])
    assert_sample(algorithm(2), staggered(), [1, 2]) # truncated q_deim
    assert_sample(algorithm(3), staggered(), [1, 2, 3])
    assert_sample(algorithm(4), staggered(), [0, 1, 2, 3])

def test_aliasing_error_bound():
    def naive_aliasing_error_bound(basis, sampling):
        sampled_basis = basis[:, sampling]
        correlation_matrix = np.dot(sampled_basis, sampled_basis.T)
        error_matrix = la.inv(correlation_matrix) - np.eye(correlation_matrix.shape[0])
        return np.sqrt(la.norm(error_matrix, ord=2))

    snapshots = np.random.rand(50, 100)
    basis = Pod.pod_basis(snapshots, 20, 1.0e-6)
    mask = np.random.permutation(basis.shape[-1])[:50]

    actual = aliasing_error_bound(basis, mask)
    expected = naive_aliasing_error_bound(basis, mask)
    npt.assert_allclose(expected, actual)


def check_empirical_operator(sampling_algorithm):
    snapshots = np.random.rand(50, 100)
    basis = Pod.pod_basis(snapshots, 20, 1.0e-6)
    sample_points = sampling_algorithm(basis)
    operator = empirical_operator(basis, sample_points)
    npt.assert_allclose(
            operator(basis[:, sample_points]),
            np.eye(basis.shape[0]), atol=5.0e-15
            )
    npt.assert_allclose(
            operator(basis[0, sample_points]),
            np.eye(basis.shape[0])[0], atol=5.0e-15
            )

def test_interpolating_empirical_operator():
    check_empirical_operator(deim)

def test_leastsquares_empirical_operator():
    sampling_algorithm = lambda b: greedy_gappy_mask(b, 30)
    check_empirical_operator(sampling_algorithm)


def test_lowest_eigenvalue_index_with_potential():
    cases = np.array([ 
            [0.1, 0.3, 0.4, 0.5],
            [0.1, 0.1, 0.3, 1.0],
            [0.2, 0.3, 0.4, 1.0],
            [0.2, 0.3, 0.4, 0.5],
            [0.2, 0.3, 0.4, 0.6]
            ])

    expected = np.array([0, 1, 2, 3, 0])

    criterion = lowest_eigenvalue_index_with_potential(0.5)
    actual = criterion(cases)
    npt.assert_array_equal(actual, expected)

def test_lowest_eigenvalue_index_with_potential_original():
    cases = np.array([
            [0.1, 0.3, 0.4, 0.5],
            [0.1, 0.1, 0.3, 1.0],
            [0.1, 0.2, 0.4, 1.0],
            [0.1, 0.2, 0.3, 0.4]
            ])

    expected = np.array([0, 1, 2, 0])

    criterion = lowest_eigenvalue_index_with_potential_original(0.5)
    actual = criterion(cases)
    npt.assert_array_equal(actual, expected)

def test_first_false_index():
    cases = np.array([
            [False, False, False],
            [False, False, True ],
            [False, True,  False],
            [False, True,  True ],
            [True, False,  False],
            [True, False,  True ],
            [True, True,   False],
            [True, True,   True ]
            ])
    expected = np.array([0, 0, 0, 0, 1, 1, 2, 3])

    actual = first_false_index(cases)
    npt.assert_array_equal(actual, expected)

def test_updated_lowest_eigenvalue():
    basis = np.eye(3, 4, dtype=np.double)
    sampling = np.arange(2)
    expected = np.array([0., 0., 1., 0.])

    actual = updated_lowest_eigenvalue(basis, sampling)
    npt.assert_array_equal(actual, expected)

def test_updated_spectrum():
    basis = np.eye(3, 4, dtype=np.double)
    sampling = np.arange(2)
    expected = np.array([[1., 1., 0.], [1., 1., 0.], [1., 1., 1.], [1., 1., 0.]])

    actual = updated_spectrum(basis, sampling)
    npt.assert_array_equal(actual, expected)

def test_updated_lowest_eigenvalue_trivial():
    basis = np.array([[0.1, 0.0, 0.3, 0.0, 0.0], [0.0, 0.2, 0.0, 0.1, 0.0]])
    sampling = np.arange(2)
    expected = np.square(np.array([0.1, 0.1, 0.2, 0.1, 0.1]))

    actual = updated_lowest_eigenvalue(basis, sampling)
    npt.assert_allclose(actual, expected)

def test_updated_lowest_eigenvalue_lower_bound_trivial():
    basis = np.array([[0.1, 0.0, 0.3, 0.0, 0.0], [0.0, 0.2, 0.0, 0.1, 0.0]])
    sampling = np.arange(2)
    expected = np.square(np.array([0.1, 0.1, 0.2, 0.1, 0.1]))

    actual = updated_lowest_eigenvalue_lower_bound(basis, sampling)
    npt.assert_allclose(actual, expected)

def test_updated_lowest_eigenvalue_lower_bound_compare():
    snapshots = np.random.rand(20, 100)
    basis = Pod.pod_basis(snapshots, 10, 1.0e-6)
    mask = q_deim(basis)

    actual = updated_lowest_eigenvalue(basis, mask)
    bound = updated_lowest_eigenvalue_lower_bound(basis, mask)
    npt.assert_allclose(bound[mask], actual[mask])
    npt.assert_allclose(bound, np.minimum(actual, bound))
    npt.assert_allclose(bound, np.maximum(bound, bound[mask[0:1]]))

def test_updated_spectrum_lower_bound_trivial():
    basis = np.array([[0.1, 0.0, 0.3, 0.0, 0.0], [0.0, 0.2, 0.0, 0.1, 0.0]])
    sampling = np.arange(2)
    expected = np.array([[0.04, 0.01], [0.04, 0.01], [0.04, 0.04], \
            [0.05, 0.01], [0.04, 0.01]])
    actual = updated_spectrum_lower_bound(basis, sampling)
    npt.assert_allclose(actual, expected)

def test_updated_spectrum_lower_bound_matches_lowest():
    snapshots = np.random.rand(20, 100)
    basis = Pod.pod_basis(snapshots, 10, 1.0e-6)
    mask = q_deim(basis)

    spectrum = updated_spectrum_lower_bound(basis, mask)
    lambda_min = updated_lowest_eigenvalue_lower_bound(basis, mask)
    npt.assert_allclose(spectrum[mask, -1], lambda_min[mask])
    npt.assert_allclose(spectrum[:, -1], lambda_min)

def test_updated_spectrum_lower_bound_compare():
    snapshots = np.random.rand(20, 100)
    basis = Pod.pod_basis(snapshots, 10, 1.0e-6)
    mask = q_deim(basis)

    actual = updated_spectrum(basis, mask)
    bound = updated_spectrum_lower_bound(basis, mask)
    npt.assert_allclose(bound[mask], actual[mask])
    npt.assert_allclose(bound[:, -1], np.minimum(actual[:, -1], bound[:, -1]))
    npt.assert_allclose(bound, np.minimum(actual, bound))
    npt.assert_allclose(bound, np.maximum(bound, bound[mask[0:1]]))
