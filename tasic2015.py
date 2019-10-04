from pathlib import Path

#extended libraries
import requests
import numpy as np
import pandas as pd
import scipy as sp
import matplotlib.pyplot as plt

#tree stuff
import anytree
import anytree.util
from anytree.importer import DictImporter
from anytree.exporter import DotExporter

#custom imports
from neurondm import *
from neurondm.models.huang2017 import Genes
from pyontutils.namespaces import ilxtr as pred
from neurondm import phenotype_namespaces as phns
from nifstd.nifstd_tools.utils import ncbigenemapping
from utils import *

###############################################################################

df_clf = pd.read_csv('Data/cell_classification.csv')
df_cls_mtd = pd.read_csv('Data/cluster_metadata.csv')
df_cell_mtd = pd.read_csv('Data/cell_metadata.csv')
df_cre_mtd  = pd.read_excel('Data/tasic_crelines.xlsx')
dendro = pd.read_csv("Data/web/big_tree_final.csv", dtype = {'position':str})

df_clf.rename(index=str, columns={"Unnamed: 0":"cell_index"}, inplace = True)
df_clf = df_clf[['cell_index','coretype', 'primary', 'secondary']]

###############################################################################

#Reorder and drop columns that are not interesting
df_cell_mtd = df_cell_mtd[['short_name','cre','major_class','sub_class',
                           'major_dissection', 'layer_dissectoin']]
df_cls_mtd = df_cls_mtd[['cluster_id','cluster_order','vignette_label',
                         'group','markers_present','markers_sparse',
                         'genes_absent','Tasic_et_al_2016_label']]

df_cls_mtd['cluster_id'] = df_cls_mtd['cluster_id'].astype(str)

#merge
df_types = df_cell_mtd.merge(df_clf, left_on='short_name', right_on='cell_index')
df_types = df_types[['short_name', 'coretype', 'primary', 'secondary',
                     'cre','major_dissection', 'layer_dissectoin']]

###############################################################################

tree_dict = {'label':'root'} #base dictionary

#parses binary node positions into a dictionary with tree structure
for label, bin_str in zip(dendro['cluster'], dendro['position']):
    parse_binary(tree_dict, bin_str, label)

importer = DictImporter()
tree = importer.import_(tree_dict)

for ind, node in enumerate(anytree.LevelOrderIter(tree)):
    node.pos = ind + 1 #node index starting at one
    if node.is_leaf:
        node.name = str(ind + 1) + " " + node.label
    else:
        node.name = str(ind + 1)

for leaf in tree.leaves: #TODO: phenotypes & deal with mismatching
    #edge cases for the last two endothelial cells
    if leaf.label  == 'Endo Tbc1d4':
        leaf.cluster_id = 'f48'
    elif leaf.label == 'Endo Myl9':
        leaf.cluster_id = 'f49'
    else:
        leaf.cluster_id = list(df_cls_mtd[df_cls_mtd['vignette_label'] == leaf.label]['cluster_id'])[0]

###############################################################################

df_types['cluster'] = df_types.apply(cluster_converter, tree = tree, axis = 1)
df_types['markers_present'] = df_types.apply(gene_merge, df = df_cls_mtd, 
                                          index = 'markers_present', axis = 1)
df_types['markers_absent'] = df_types.apply(gene_merge, df = df_cls_mtd,
                                             index = 'genes_absent', axis = 1)

###############################################################################

