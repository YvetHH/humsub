import pandas as pd
import numpy as np
import matplotlib.pylab as plt
import seaborn as sns
import os, sys


TAXONMIC_LEVELS=['kindom','phylum','class','order','family','genus','species']

def tax2table(Taxonomy_Series,split_character=';',remove_prefix=False):
    """
        Transforms (green_genes) taxonomy to a table
        Expect the following input format:
        d__Bacteria;p__Bacteroidota;c__Bacteroidia;f__

        Replaces empty values and can remove prefix 'c__'

    """

    if Taxonomy_Series.isnull().any():
        warnings.warn("Some samples have no taxonomy asigned based on checkm. Samples:\n"+ \
                    ', '.join(Taxonomy_Series.index[Taxonomy_Series.isnull()])
                    )
        Taxonomy_Series= Taxonomy_Series.dropna().astype(str)

    Tax= pd.DataFrame(list(  Taxonomy_Series.apply(lambda s: s.split(split_character))),
                       index=Taxonomy_Series.index)


    Tax.columns= TAXONMIC_LEVELS[:len(Tax.columns)]

    if remove_prefix:
        Tax=Tax.applymap(lambda s: s[3:]).replace('',pd.NA)
    else:
        Tax[Tax.applymap(len)==3]=np.nan

    return Tax




def load_subspecies_attribution(species_rep):
    """
    opens subspecies pandas Series from relative path.
    
    Input: species representative id, e.g. GUT_GENOME0000001
    """

    Subspecies = pd.read_table(f'../data/subspecies_tables/{species_rep}.tsv',index_col=0).squeeze()

    Subspecies=Subspecies.str.rsplit(expand=True).iloc[:,-1]
    Subspecies.name= 'Subspecies'
    return Subspecies









## Statistic scripts


from scipy import stats
def evaluate_count_table(CountMatrix):
    g, p, dof, Expected =stats.chi2_contingency(CountMatrix)
    return p


# TODO: for now this script only displays the tables, display is a function specific to ipython, I think.
def subevaluation(CountMatrix):

    from itertools import combinations
    for headers in combinations(CountMatrix.columns,2):
        for index in combinations(CountMatrix.index,2):

            sub_matrix= CountMatrix.loc[index,headers]

            #if all(sub_matrix.sum(1)>0) and all(sub_matrix.sum(0)>0):
            
            _,p= stats.fisher_exact(sub_matrix) #evaluate_count_table(sub_matrix)

            if p<0.05:

                display(sub_matrix)
                print(f"p={p:.2g}")

                    
                    
## shannon
# this is copied from scikit bio


def _validate(counts, suppress_cast=False):
    """Validate and convert input to an acceptable counts vector type.
    Note: may not always return a copy of `counts`!
    """
    counts = np.asarray(counts)

    if not suppress_cast:
        counts = counts.astype(int, casting='safe', copy=False)

    if counts.ndim != 1:
        raise ValueError("Only 1-D vectors are supported.")
    elif (counts < 0).any():
        raise ValueError("Counts vector cannot contain negative values.")

    return counts



def shannon(counts, base=2):
    """Calculate Shannon entropy of counts (H), default in bits.
    Parameters
    ----------
    counts : 1-D array_like, int
        Vector of counts.
    base : scalar, optional
        Logarithm base to use in the calculations.
    Returns
    -------
    double
        Shannon diversity index H.
    Notes
    -----
    The implementation here is based on the description given in the SDR-IV
    online manual [1]_, except that the default logarithm base used here is 2
    instead of :math:`e`.
    References
    ----------
    .. [1] http://www.pisces-conservation.com/sdrhelp/index.html
    """
    counts = _validate(counts,suppress_cast=True)
    freqs = counts / counts.sum()
    nonzero_freqs = freqs[freqs.nonzero()]
    return -(nonzero_freqs * np.log(nonzero_freqs)).sum() / np.log(base)



                    
                    
                    
## Script run at start

#load taxonomy
#genome_metadata= pd.read_table('../data/genomes-nr_metadata.tsv',index_col=0,low_memory=False)
#genome_metadata.eval('QualityScore= Completeness - 5* Contamination',inplace=True)
# get tax table
#Tax= tax2table(genome_metadata.Lineage,remove_prefix=True)






# A function to calculate some statistics 

from itertools import combinations
def N(x): return np.sum(~np.isnan(x))
def Q1(x): return np.percentile(x,25)
def Q3(x): return np.percentile(x,75)


def Sumarize_Table(Data,grouping_variable, order=None ,
                   aggregation_functions=[N,'mean','std',Q1,'median', Q3 ],
                   calculate_diff_for=['median','mean'] ):
    """ Sumarize a datatable grouped by a grouping variable using the functions in aggregation function

        eg. Row 1 mean std Q1 median Q2

        If a table with Pairwise p-values is given it is added to the table
        and the differences for the variables in the list calculate_diff_for are calculated

    """

    if order is None:
        order=np.unique(grouping_variable)

    Summary=Data.groupby(grouping_variable).aggregate(aggregation_functions).T.unstack()

    Summary= Summary[order]

    for measure in calculate_diff_for:
        assert measure in aggregation_functions, "Need {measure} values for {measure} difference calculations".format(measure=measure)

        for (G1,G2) in combinations(order,2):
            comparison='{} <-> {}'.format(G1,G2)
            Summary[(measure+'_diff',comparison)] = Summary[(G2,measure)]-Summary[(G1,measure)]

    return Summary