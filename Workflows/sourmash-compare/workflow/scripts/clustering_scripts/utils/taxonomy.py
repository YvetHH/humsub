
import pandas as pd
import numpy as np


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



def table2tax(Tax,split_character=';',kindom_name='domain'):
    "Transforms taxonomy to gg taxonomy d__Bacteria;o__ ..."
    save_Tax=Tax.rename(columns={'kindom':kindom_name}).fillna('')

    for level,values in save_Tax.T.iterrows():

        save_Tax [level]= level[0]+'__'+values

    save_Tax = save_Tax.apply(split_character.join,axis=1)

    return save_Tax
