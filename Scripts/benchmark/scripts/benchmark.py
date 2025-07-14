import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
from scipy.spatial import distance
from statsplot import transformations
from sklearn import metrics
from typing import Union
from skbio import OrdinationResults
from skbio.stats.distance import DistanceMatrix
from gemelli.matrix_completion import MatrixCompletion
from gemelli import rpca
from scipy.linalg import svd


def standardize_input(
    test: pd.Series, reference: pd.Series, subspecies=True
) -> pd.DataFrame:
    """
    This function standardizes names of of test and reference series for further use.
    Returns:
        df: pd.Dataframe containing with shape (n,2), n-number of taxa.
    
    Arguments:
        test: pd.Series. Contains same taxonomic ids as in reference 
        as index and abundances as values.
        
        reference: pd.Series. Contains same taxonomic ids as in test
        as index and abundances as values.
        
        subspecies: boolean. True if want to calculate on subspecies level.
    """
    test.name = "test"
    reference.name = "reference"
    df = pd.concat([test, reference], axis=1).fillna(0)
    if subspecies == False:
        df.index = [x[:4] for x in df.index]
        df = df.T.groupby(level=0, axis=1).sum().T
    return df


def get_metrics(
    test: pd.Series, reference: pd.Series, subspecies=True, abund_threshold=0
):
    """
    This function takes abundances obtained using the method getting tested
    and the reference abundance.
    Returns:
        precision: float
        recall: float
        F1 score: float
    
    Arguments:
        test: pd.Series. Contains same taxonomic ids as in reference 
        as index and abundances as values.
        
        reference: pd.Series. Contains same taxonomic ids as in test
        as index and abundances as values.
        
        subspecies: boolean. True if want to calculate on subspecies level.
        
        abund_threshold: float. Abundance threshold on which to calculate metrics.
        It will exclude all taxa with abundance < abund_threshold in test data.
        Reference values will stay the same.
    """
    df = standardize_input(test, reference, subspecies=subspecies)
    df.loc[df.test < abund_threshold, "test"] = 0
    TP = ((df.test > 0) & (df.reference > 0)).sum()
    FP = ((df.test > 0) & (df.reference == 0)).sum()
    FN = ((df.test == 0) & (df.reference > 0)).sum()
    prec = TP / (TP + FP)
    recall = TP / (TP + FN)
    F1 = 2 * ((prec * recall) / (prec + recall))
    if np.isnan(prec):
        prec = 1.0
    if np.isnan(recall):
        recall = 1.0
    if np.isnan(F1):
        F1 = 1.0
    return (prec, recall, F1)

def l2_distance(test: pd.Series, reference: pd.Series, subspecies=True):
    
    df = standardize_input(test, reference, subspecies=subspecies)
    # Calculate the Euclidean (L2) distance between the two columns
    df = df.loc[df.sum(1) > 0]
    l2_dist = distance.euclidean(df['reference'], df['test'])
    return l2_dist

def calculate_auc(
    test: pd.Series,
    reference: pd.Series,
    subspecies=True,
    abundance_step=0.001,
    upper_abund_limit=0.1,
    plot=True,
) -> float:
    """
    Calculate area under the curve of precision-recall plot.
    
    Returns:
        auc: float. Area under the curve of of precision-recall plot.
    
    Arguments:
        test: pd.Series. Contains same taxonomic ids as in reference 
        as index and abundances as values.
        
        reference: pd.Series. Contains same taxonomic ids as in test
        as index and abundances as values.
        
        subspecies: boolean. True if want to calculate on subspecies level.
        
        abundance_step: float. How big are steps for changing the abundance threshold between 0 and 1.

        upper_abund_limit: float. Upper limit of abundance tested.
        
        plot: boolean. If True, plot of precision-recall curve.
    """
    df = standardize_input(test, reference, subspecies=subspecies)
    prec_list = []
    recall_list = []
    for abund_threshold in np.arange(0, float(upper_abund_limit), abundance_step):
        prec, recall, _ = get_metrics(
            test=df["test"],
            reference=df["reference"],
            abund_threshold=abund_threshold,
            subspecies=subspecies,
        )
        prec_list.append(prec)
        recall_list.append(recall)
    if plot:
        ax = sns.lineplot(
            x=list(reversed(recall_list)), y=list(reversed(prec_list))
        )
        ax.set_title("Precision-recall curve")
        ax.set_xlabel("Recall")
        ax.set_ylabel("Precision")
    auc = metrics.auc(recall_list, prec_list)
    return auc