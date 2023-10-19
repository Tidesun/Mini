
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