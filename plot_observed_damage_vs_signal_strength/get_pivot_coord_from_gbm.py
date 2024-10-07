import numpy as np
from scipy.stats import gmean
import pandas as pd

#Code written by Daniel Svensson: https://github.com/danisven

def PLR_single(mat, i, D):
    """
    Function to calculate pivot logratios for a single variable.

    Divides the variable at the i:th index with the geometric mean of the rest of the variables in the matrix (mat).
    Then calculates the log of that ratio (logratio). Then multiplies the logratio with a normalizing constant (norm).

    INPUT
    mat: a numpy 2D array with variables as columns and observations as rows
    i: the index of the column (0-based) for which pivot logratios are to be calculated

    RETURNS
    plr: a numpy array containing the pivot logratios of the column at the i:th index in mat.
    """

    # The variable to calculate pivot logratios for (the column "to the left")
    pivot_col = mat[i]

    # All the other columns in the matrix (the columns "to the right")
    r_cols = np.delete(mat, i, axis=0)

    # The geometric mean of the other columns
    r_cols_gm = gmean(r_cols)

    # The normalizing constant
    norm = np.sqrt((D - 1) / (D))

    # The plr
    plr = norm * np.log(pivot_col / r_cols_gm)

    return plr


def PLR(mat):
    """
    Function to create pivot logratios of all variables in a 2D array.

    For each variable in the input matrix, PLR calls PLR_single to calculate the
    pivot logratios.

    INPUT
    mat: a numpy 2D array with variables as columns and observations as rows

    OUTPUT
    mat_plr: a numpy 2D array ordered the same way and with the same dimensions
             as the input matrix (mat).
    """

    # The number of variables
    D = mat.shape[0]

    # The list to hold all the variable's pivot logratios
    plr_list = []

    # Loop over each variable in the input matrix and calculate its pivot logratios, then append to plr_list
    for i in range(0, D):
        plr_i = PLR_single(mat, i, D)
        plr_list.append(plr_i)

    # Return the result as an nparray
    return np.array(plr_list).T

input_file = '../kraken_reports/all_weeks_count_values_removed_zero_inflated_gbm.csv'
output_file = '../kraken_reports/all_weeks_count_values_removed_zero_inflated_gbm_pivot_coord.csv'
# Read the data
df_data = pd.read_csv(input_file, sep=',')
tax_ids = df_data.columns
# Convert the dataframe to a 2D numpy array
df_arr = df_data.T.to_numpy()
# Run the PLR transformation for the whole matrix
plr_arr = PLR(df_arr)
# Convert to pandas dataframe if you want
plr_df = pd.DataFrame(plr_arr)
plr_df.columns = tax_ids
#Save file
plr_df.to_csv(output_file, sep=',')