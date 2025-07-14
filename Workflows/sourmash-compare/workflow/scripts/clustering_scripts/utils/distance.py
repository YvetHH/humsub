import pandas as pd
import numpy as np


def pairewise_id2matrix(pairewise_id, fillna=0):
    """
        This functions turns a pairewise genome distance series [genome1, genome2, id...]
        In to a matrix [genome1 genome2] of the values of column.
        When ANI values are symetrical (with minimal error),
        usually only one halve of NxN possibilities values are calculated.


        Diagonal values are set to 1

    """

    assert pairewise_id.min()>=0
    assert pairewise_id.max()<=1

    ID = pairewise_id.unstack()


    all_indexes = ID.index.union(ID.columns)

    ID = ID.reindex(index=all_indexes, columns=all_indexes)
    ID = ID.fillna(0)
    ID = ID + ID.T
    ID.values[np.eye(ID.shape[0], dtype=bool)] = 1
    return ID.replace(0, fillna)
