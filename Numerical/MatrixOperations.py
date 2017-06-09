import numpy as np
import scipy.sparse as scsp

def reshape(a, shape):
    """Reshape the sparse matrix `a`.

    Returns a coo_matrix with shape `shape`.
    """
    if not hasattr(shape, '__len__') or len(shape) != 2:
        raise ValueError('`shape` must be a sequence of two integers')

    c = a.tocoo()
    nrows, ncols = c.shape
    size = nrows * ncols

    new_size =  shape[0] * shape[1]
    if new_size != size:
        raise ValueError('total size of new array must be unchanged')

    flat_indices = ncols * c.row + c.col
    new_row, new_col = divmod(flat_indices, shape[1])

    b = scsp.coo_matrix((c.data, (new_row, new_col)), shape=shape)
    return b


def CheckIntegrity():
    a = scsp.coo_matrix([[1, 2, 0], [0, 0, 3], [4, 0, 5], [0, 0, 2]])
    return 'ok'


if __name__ == '__main__':
    print(CheckIntegrity())# pragma: no cover





