"""Flow cytometry helper functions

"""
import cytoflow as flow

def subset_exp(exp, subset_list = []):
    """Applies subsetting function to experiment if subset is provided.

    Parameters
    ----------
    exp : Cytoflow experiment
        Experiment to be subsetted.
    subset_list : List, optional
        List of boolean statements to be ANDed together
    """
    # If the list has anything in it... and them all together
    if len(subset_list) != 0:
        subset_expression = ' and '.join(subset_list)
        exp = exp.query(subset_expression)

    return exp
