
from scipy.sparse import coo_matrix
def convert_dict_to_sparse_matrix(data_dict,num_rows,num_cols,dtype=None):
    # Get the row, column, and data arrays
    rows, cols, data = zip(*[(k[0], k[1], v) for k, v in data_dict.items()])

    # Create the sparse matrix using the coo_matrix function
    if dtype is None:
        sparse_matrix = coo_matrix((data, (rows, cols)), shape=(num_rows, num_cols))
    else:
        sparse_matrix = coo_matrix((data, (rows, cols)), shape=(num_rows, num_cols),dtype=dtype)
    return sparse_matrix
def safe_divide(numerator,denominator):
    denominator[denominator==0] = 1
    return numerator/denominator
def safe_divide_sparse(numerator,denominator):
    denominator[denominator==0] = 1
    return numerator.multiply(1/denominator)
import numpy as np
import scipy.sparse as sp

def sp_unique(sp_matrix, axis=0):
    ''' Returns a sparse matrix with the unique rows (axis=0)
    or columns (axis=1) of an input sparse matrix sp_matrix'''
    if axis == 1:
        sp_matrix = sp_matrix.T

    old_format = sp_matrix.getformat()
    dt = np.dtype(sp_matrix)
    ncols = sp_matrix.shape[1]

    if old_format != 'lil':
        sp_matrix = sp_matrix.tolil()

    _, ind = np.unique(sp_matrix.data + sp_matrix.rows, return_index=True)
    rows = sp_matrix.rows[ind]
    data = sp_matrix.data[ind]
    nrows_uniq = data.shape[0]

    sp_matrix = sp.lil_matrix((nrows_uniq, ncols), dtype=dt)  #  or sp_matrix.resize(nrows_uniq, ncols)
    sp_matrix.data = data
    sp_matrix.rows = rows

    ret = sp_matrix.asformat(old_format)
    if axis == 1:
        ret = ret.T        
    return ret