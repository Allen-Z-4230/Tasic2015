#!/usr/bin/env python
# coding: utf-8

# In[1]:


#standard libraries
from pathlib import Path
from importlib import reload

#extended libraries
import requests
import rdflib
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
from pyontutils.closed_namespaces import rdf, rdfs, owl
from nifstd.nifstd_tools.utils import ncbigenemapping

#local imports
from utils import *
from tasic2015 import *


# In[2]:


get_ipython().run_line_magic('load_ext', 'autoreload')
get_ipython().run_line_magic('autoreload', '1')
get_ipython().run_line_magic('aimport', 'tasic2015')


# In[4]:


def cellguard(addns=False):
    # unfortunately ipy.hooks['pre_run_code_hook'].add(__cellguard)
    # causes this to be called too frequently :/
    setLocalNames()
    setLocalContext()
    if addns:
        setLocalNames(phns.BBP)


# In[5]:


df_clf = pd.read_csv('Data/cell_classification.csv')
df_cls_mtd = pd.read_csv('Data/cluster_metadata.csv')
df_cell_mtd = pd.read_csv('Data/cell_metadata.csv')

df_cre_mtd  = pd.read_excel('Data/tasic_crelines.xlsx')

dendro = pd.read_csv("Data/web/big_tree_final.csv", dtype = {'position':str})


# In[6]:


"""
Information about the cluster membership of each cell,
including whether the cell is a "core" (unambiguously assigned to a single cluster) 
or "transition" (shares membership between two clusters) cell, 
as well as its membership score (from 0-10) for each cluster (labeled f01 to f49). 
"""
df_clf.head(3)
#TODO: Only including unambiguous cells with f-values = 10 for one cluster only.


# In[7]:


"""
 Information about each data-driven cluster, including its label, 
 the corresponding label in Tasic et al. (Nat. Neuro, 2106), the primary cell class 
 membership, and marker genes (including genes with widespread expression in the cluster, 
 sparse expression in the cluster, and no expression in the cluster). 
"""
df_cls_mtd.head(3)


# In[8]:


"""
Information about each cell profiled, including its nomenclature, 
Cre line of origin, dissection, date of collection and sequencing, 
and read mapping statistics
"""
df_cell_mtd.head(3)


# In[9]:


df_clf.rename(index=str, columns={"Unnamed: 0":"cell_index"}, inplace = True)
df_clf = df_clf[['cell_index','coretype', 'primary', 'secondary']]

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


# In[10]:


#cre line conversion metadata
df_cre_mtd.head(3)


# In[11]:


#final combined table for cluster, layer, & cre
df_types.head(3)


# main dataframe: df_types
# 
# metadata: cre_mtd & cls_mtd

# ## Dendrogram Parsing

# In[12]:


dendro.head()


# In[13]:


tree_dict = {'label':'root'} #base dictionary

#parses binary node positions into a dictionary with tree structure
for label, bin_str in zip(dendro['cluster'], dendro['position']):
    parse_binary(tree_dict, bin_str, label)


# In[14]:


importer = DictImporter()
tree = importer.import_(tree_dict)

for ind, node in enumerate(anytree.LevelOrderIter(tree)):
    node.pos = ind + 1 #node index starting at one
    if node.is_leaf:
        node.name = str(ind + 1) + " " + node.label
    else:
        node.name = str(ind + 1)
        
print(anytree.RenderTree(tree))
DotExporter(tree).to_picture("dendro.png")

#node.name for rendering
#node.pos & node.label used for parsing


# In[15]:


for leaf in tree.leaves: #TODO: phenotypes & deal with mismatching
    #edge cases for the last two endothelial cells
    if leaf.label  == 'Endo Tbc1d4':
        leaf.cluster_id = 'f48'
    elif leaf.label == 'Endo Myl9':
        leaf.cluster_id = 'f49'
    else:
        leaf.cluster_id = list(df_cls_mtd[df_cls_mtd['vignette_label'] == leaf.label]['cluster_id'])[0]
    
    #FIXME: label name matching somehow not working
    """
    genes_present = list(df_cls_mtd[df_cls_mtd['cluster_id'] == leaf.cluster_id]['markers_present'].str.split(","))
    genes_absent = list(df_cls_mtd[df_cls_mtd['cluster_id'] == leaf.cluster_id]['genes_absent'].str.split(","))
    if genes_present == np.nan:
        leaf.genes_present = set()
    if genes_absent == np.nan:
        leaf.genes_absent = set()
    else:
        leaf.genes_present = set(genes_present[0])
        #leaf.genes_absent = set(genes_absent)
    """


# In[16]:


#gene and cluster columns converted from dendrogram structure
df_types['cluster'] = df_types.apply(cluster_converter, tree = tree, axis = 1)
df_types['markers_present'] = df_types.apply(gene_merge, df = df_cls_mtd, index = 'markers_present', axis = 1)
df_types['markers_absent'] = df_types.apply(gene_merge, df = df_cls_mtd, index = 'genes_absent', axis = 1)


# In[17]:


df_types[df_types['coretype'] == 'Transition'].head(3)


# In[18]:


response = requests.get('http://api.brain-map.org/api/v2/data/query.json?criteria='
                        'model::TransgenicLine,rma::options[num_rows$eqall]')
cre_ref = pd.DataFrame(response.json()['msg'])
cre_ref['stock_number'] = pd.to_numeric(cre_ref['stock_number'])

#creates final cre metadata based on a name merge & a stock number merge
cre_df1 = pd.merge(df_cre_mtd, cre_ref,  how='inner',
                   left_on='Driver Line',right_on = 'name')

cre_df2 = pd.merge(df_cre_mtd.dropna(subset = ['Public Repository Stock #']),
                   cre_ref.dropna(subset = ['stock_number']),  how='inner',
                   left_on='Public Repository Stock #',
                   right_on = 'stock_number')

#drop duplicates
cre_df1.drop_duplicates(inplace = True)
cre_df2.drop_duplicates(inplace = True)
cre_df = pd.concat([cre_df1, cre_df2], axis = 0)
cre_df.drop_duplicates(inplace = True)

#reorder
cre_df = cre_df[['Abbreviation','name', 'id', 'stock_number',
                 'transgenic_line_source_name',
                 'transgenic_line_type_name',
                 'url_prefix', 'url_suffix','description']]


# ### Data Visualization

# In[19]:


pd.value_counts(df_types['coretype']).plot.barh();
plt.title("Tasic 2015 Cell Type Clusters");


# In[20]:


pd.value_counts(df_types['layer_dissectoin']).plot.barh(figsize = (10,5));
plt.title("Layer of Dissection distributions");


# In[21]:


pd.value_counts(df_types['cre']).plot.barh(figsize = (10,5));
plt.title("Mouse Cre Line Distributions");


# ### Phenotype Bagging

# In[22]:


#Tables should now be in final version
df_types.head(3)


# In[23]:


df_cre_mtd[df_cre_mtd['Abbreviation'] == 'Ndnf']


# In[24]:


len(df_types['cre'].value_counts())


# In[25]:


len(df_cre_mtd)


# In[26]:


#We're supposed to have 29 cre-line types
len(cre_df)


# In[28]:


pd.diff(cre_df['Abbreviation'], df_types['cre'])


# In[44]:


b = sorted(cre_df['Abbreviation'].items())
b


# In[45]:


a = sorted(df_types['cre'].value_counts().items())
a


# In[47]:


list(zip(a,b))


# ## Writing & Serializing

# In[27]:


ttl_test_path = '/var/host/media/removable/SD Card/Neuron/Tasic/ttl_export'
config = Config("tasic-2015", ttl_export_dir=Path(ttl_test_path))
tb = TasicBagger(data = df_types, **metadata)
ind = 0
for label, bag in tb.bags:
    if ind == 1:
        break
    TasicNeuron(*bag, label = label, override=True)
    ind += 1
config.write()
config.write_python()


# In[ ]:


from neurondm.sheets import Sheet
from neurondm import OntId, OntTerm, Config, NeuronEBM, Neuron
from pyontutils.utils import byCol, relative_path
from pyontutils.namespaces import ilxtr
from pyontutils.closed_namespaces import rdfs

main()


# In[ ]:


metadata = {"cre":cre_df}
tb = TasicBagger(data = df_types, **metadata)
tb.build_transgenic_lines()


# ### Test Code

# In[ ]:


# metada query test
df_cre_mtd[df_cre_mtd['Abbreviation'] == df_types['cre'][12]]


# In[ ]:


def main():
    metadata = {"cre": cre_df}
    ttl_test_path = '/mnt/c/Users/allen/Desktop/Neuron/Tasic2015/ttl_export'
    config = Config("tasic-2015", ttl_export_dir=Path(ttl_test_path))
    tb = TasicBagger(data=df_types, **metadata)
    ind = 0
    for label, bag in tb.bags:
        if ind == 10:
            break
        TasicNeuron(*bag, label=label, override=True)
        ind += 1
    config.write()


# In[ ]:


with Tasic2015:
    print(Mouse)


# In[ ]:


simpleOnt()


# In[ ]:


Phenotype("ilxtr:cluster6", "ilxtr:hasComputedMolecularPhenotype", label = "Tasic2015", check = False)


# In[ ]:


from pyontutils.utils import byCol, relative_path
from pyontutils.namespaces import ilxtr
from pyontutils.closed_namespaces import rdfs


# In[ ]:


with Neuron(phns.Species.Mouse, phns.Regions.V1) as context:
    print(context)
    n11 = Neuron(phns.Layers.L1)
    print(n11)

