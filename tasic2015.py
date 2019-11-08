from pyontutils.config import devconfig
from pyontutils.closed_namespaces import rdf, rdfs, owl
from pyontutils.namespaces import makePrefixes, ilxtr, definition
from pathlib import Path
import requests

# extended libraries
import numpy as np
import pandas as pd
import scipy as sp
import rdflib
import matplotlib.pyplot as plt

# tree stuff
import anytree
import anytree.util
from anytree.importer import DictImporter
from anytree.exporter import DotExporter

# custom imports
from neurondm import *
from pyontutils.core import simpleOnt
from neurondm.models.huang2017 import Genes  # FIXME: This is failing
from pyontutils.namespaces import ilxtr as pred
from neurondm import phenotype_namespaces as phns
from nifstd.nifstd_tools.utils import ncbigenemapping
from utils import *

###############################################################################

#Read in dataframes
df_clf = pd.read_csv('Data/cell_classification.csv')
df_cls_mtd = pd.read_csv('Data/cluster_metadata.csv')
df_cell_mtd = pd.read_csv('Data/cell_metadata.csv')
df_cre_mtd = pd.read_excel('Data/tasic_crelines.xlsx')
dendro = pd.read_csv("Data/web/big_tree_final.csv", dtype={'position': str})

df_clf.rename(index=str, columns={"Unnamed: 0": "cell_index"}, inplace=True)
df_clf = df_clf[['cell_index', 'coretype', 'primary', 'secondary']]

###############################################################################

# Reorder and drop columns that are not interesting
df_cell_mtd = df_cell_mtd[['short_name', 'cre', 'major_class', 'sub_class',
                           'major_dissection', 'layer_dissectoin']]
df_cls_mtd = df_cls_mtd[['cluster_id', 'cluster_order', 'vignette_label',
                         'group', 'markers_present', 'markers_sparse',
                         'genes_absent', 'Tasic_et_al_2016_label']]

df_cls_mtd['cluster_id'] = df_cls_mtd['cluster_id'].astype(str)

# merge
df_types = df_cell_mtd.merge(df_clf, left_on='short_name', right_on='cell_index')
df_types = df_types[['short_name', 'coretype', 'primary', 'secondary',
                     'cre', 'major_dissection', 'layer_dissectoin']]

###############################################################################

# create breath-first binary search tree with tasic clusters
tree_dict = {'label': 'root'}  # base dictionary

# parses binary node positions into a dictionary with tree structure
for label, bin_str in zip(dendro['cluster'], dendro['position']):
    parse_binary(tree_dict, bin_str, label)

importer = DictImporter()
tree = importer.import_(tree_dict)

for ind, node in enumerate(anytree.LevelOrderIter(tree)):
    node.pos = ind + 1  # node index starting at one
    if node.is_leaf:
        node.name = str(ind + 1) + " " + node.label
    else:
        node.name = str(ind + 1)

DotExporter(tree).to_picture("dendro.png")

for leaf in tree.leaves:  # TODO: phenotypes & deal with mismatching
    # edge cases for the last two endothelial cells
    if leaf.label == 'Endo Tbc1d4':
        leaf.cluster_id = 'f48'
    elif leaf.label == 'Endo Myl9':
        leaf.cluster_id = 'f49'
    else:
        leaf.cluster_id = list(
            df_cls_mtd[df_cls_mtd['vignette_label'] == leaf.label]['cluster_id'])[0]

###############################################################################

# Reads in computed molecular phenotypes (cluster & genes)
df_types['cluster'] = df_types.apply(cluster_converter, tree=tree, axis=1)
df_types['markers_present'] = df_types.apply(gene_merge, df=df_cls_mtd,
                                             index='markers_present', axis=1)
df_types['markers_absent'] = df_types.apply(gene_merge, df=df_cls_mtd,
                                            index='genes_absent', axis=1)

###############################################################################

# Extracts cre line information from a variety of tables
response = requests.get('http://api.brain-map.org/api/v2/data/query.json?criteria='
                        'model::TransgenicLine,rma::options[num_rows$eqall]')
cre_ref = pd.DataFrame(response.json()['msg'])
cre_ref['stock_number'] = pd.to_numeric(cre_ref['stock_number'])

# creates final cre metadata based on a name merge & a stock number merge
cre_df1 = pd.merge(df_cre_mtd, cre_ref, how='inner',
                   left_on='Driver Line', right_on='name')

cre_df2 = pd.merge(df_cre_mtd.dropna(subset=['Public Repository Stock #']),
                   cre_ref.dropna(subset=['stock_number']), how='inner',
                   left_on='Public Repository Stock #',
                   right_on='stock_number')

# drop duplicates
cre_df1.drop_duplicates(inplace=True)
cre_df2.drop_duplicates(inplace=True)
cre_df = pd.concat([cre_df1, cre_df2], axis=0)
cre_df.drop_duplicates(inplace=True)

# reorder
cre_df = cre_df[['Abbreviation', 'name', 'id', 'stock_number',
                 'transgenic_line_source_name',
                 'transgenic_line_type_name',
                 'url_prefix', 'url_suffix', 'description']]

###############################################################################


class Tasic2015(Genes, phns.Species,
                phns.Regions, phns.Layers):
    branch = devconfig.neurons_branch

    # TODO: Ambiguous, most likely post-clustering layers
    L6a = Phenotype(ilxtr.TasicL6a, ilxtr.hasSomaLocatedIn,
                    label="Tasic Layer VI - A", override=True)
    L6b = Phenotype(ilxtr.TasicL6b, ilxtr.hasSomaLocatedIn,
                    label="Tasic Layer VI - B", override=True)

    # Aggregate Layer Phenotypes
    with phns.Layers:
        upper = LogicalPhenotype(OR, L1, L2, L3)
        lower = LogicalPhenotype(OR, L4, L5, L6)
        All = LogicalPhenotype(OR, upper, lower)

    # Tasic Hiearchical Cluster Position (1-indexed, breadth first)
    cmp = pred.hasComputedMolecularPhenotype


class TasicBagger:
    # from Allen Cell Types
    branch = devconfig.neurons_branch
    prefixes = {**{'JAX': 'http://jaxmice.jax.org/strain/',
                   'MMRRC': 'http://www.mmrrc.org/catalog/getSDS.jsp?mmrrc_id=',
                   'AllenTL': 'http://api.brain-map.org/api/v2/data/TransgenicLine/'},
                **makePrefixes('definition', 'ilxtr', 'owl')}

    def __init__(self, data=None, **metadata):
        """
        Initializes a Tasic Neuron Bagger Object

        data: Pandas DataFrame
        **metadata: any relevant metadata DataFrames
        """
        self.data = data
        self.ns = {k: rdflib.Namespace(v) for k, v in self.prefixes.items()}
        if metadata:
            for key, dataframe in metadata.items():
                setattr(self, key, dataframe)

    @staticmethod
    def layer_parse(key):
        """method that parses the layer labels for Tasic 2015 cells.
        returns the layer phenotype.

        key: str, layer label
        """
        if key == 'L2/3':
            layer_phn = L23
        else:
            layer_phn = Tasic2015[key]

        return layer_phn

    @staticmethod
    def gene_parse(gene_set, mode='present'):
        """method that parses a set of gene markers and returns
        a list of phenotypes. Genes that are not mapped properly are printed.

        gene_set: set, set of gene markers in dtype str
        mode: str, "present" or "absent" markers.
        """
        gene_phns = []
        undef = set()
        f = Phenotype if mode == 'present' else NegPhenotype
        for gene in gene_set:
            if gene in Genes.__dict__:
                gene_phns.append(f(Genes[gene]))
            else:
                undef.add(gene)
        if len(undef) > 0:  # FIXME: some genes cannot be mapped
            mappings, to_add, errors = ncbigenemapping(undef)
            for gene_name, iri in mappings.items():
                gene_phns.append(f(iri, ilxtr.hasExpressionPhenotype,
                                   label=gene_name, override=True))
        if len(gene_phns) != len(gene_set):
            print(errors)
        return gene_phns

    @staticmethod
    def cluster_parse(pos):
        """Tasic Computed Gene-based hiearchical cluster
        pos: int
        """
        cmp = pred.hasComputedMolecularPhenotypeFromRNA
        cluster_phn = Phenotype("ilxtr:cluster" + str(pos), cmp,
                                label=str(pos))
        return cluster_phn

    # def transgenic_parse(self, row):
    #     phenotypes = []
    #     pred = 'ilxtr:hasDriverExpressionPhenotype'
    #     cre = self.cre_metadata[self.cre_metadata['Abbreviation'] == row['cre']]
    #     cre_phn = Phenotype('PR:000013502', pred.hasExpressionPhenotype)
    #     for tl in cell_line['donor']['transgenic_lines']:
    #     prefix = self.cre_metadata['Driver Line']
    #     suffix = self.cre_metadata['Public Repository Stock #'] if tl['stock_number'] else str(
    #         tl['id'])
    #     line_names = []
    #     if prefix and suffix and prefix in ['AIBS', 'MMRRC', 'JAX']:
    #     if prefix == 'AIBS':
    #     prefix = 'AllenTL'
    #     iri = self.ns[prefix][suffix]
    #     phenotypes.append(Phenotype(iri, pred))
    #     return cre_phn

    def build_transgenic_lines(self):
        triples = []
        for ind, tl in self.cre.iterrows():
            _id = tl['stock_number'] if tl['stock_number'] else tl['id']
            prefix = tl['transgenic_line_source_name']
            line_type = tl['transgenic_line_type_name']
            if prefix not in ['JAX', 'MMRRC', 'AIBS']:
                print('WARNING:', 'unknown prefix')
                continue
            elif prefix == 'AIBS':
                prefix = 'AllenTL'

            _class = self.ns[prefix][str(_id)]
            triples.append((_class, rdf.type, owl.Class))
            triples.append((_class, rdfs.label, rdflib.Literal(tl['name'])))
            triples.append((_class, definition, rdflib.Literal(tl['description'])))
            triples.append((_class, rdfs.subClassOf, ilxtr.transgenicLine))
            triples.append((_class, ilxtr.hasTransgenicType, ilxtr[line_type + 'Line']))

        # TODO aspects.ttl?
        transgenic_lines = simpleOnt(filename='tasic-transgenic-lines',
                                     path='ttl/generated/',
                                     prefixes=self.prefixes,
                                     triples=triples,
                                     comment='Tasic transgenic lines for cell types',
                                     branch=self.branch)

        transgenic_lines._graph.write()

        return transgenic_lines

    @property
    def bags(self):
        with Tasic2015:
            # Every neuron sampled in this paper is from V1
            with Neuron(phns.Species.Mouse, phns.Regions.V1) as context:
                self.build_transgenic_lines()
                for row in self.data.itertuples():
                    label = str(row.short_name)  # change this
                    cluster_phn = TasicBagger.cluster_parse(row.cluster)
                    layer_phn = TasicBagger.layer_parse(row.layer_dissectoin)
                    cre_phn = None  # FIX after building the cre lines

                    present_phns = TasicBagger.gene_parse(row.markers_present, mode='present')
                    absent_phns = TasicBagger.gene_parse(row.markers_absent, mode='absent')
                    markers_phns = present_phns + absent_phns

                    #phenotypes = [layer_phn, cre_phn, cluster_phn] + markers_phns
                    phenotypes = [layer_phn, cluster_phn] + markers_phns
                    yield label, phenotypes


class TasicNeuron(NeuronEBM):
    owlClass = ilxtr.NeuronTasic2015
    shortname = 'Tasic2015'


def main(stop=None):
    setattr(Phenotype._predicates,
            'hasComputedMolecularPhenotype',
            ilxtr.hasComputedMolecularPhenotype)  # adds the custom phenotype
    metadata = {"cre": cre_df}
    ttl_test_path = '/mnt/c/Users/allen/Desktop/Neuron/Tasic2015/ttl_export'
    config = Config("tasic-2015", ttl_export_dir=Path(ttl_test_path))
    tb = TasicBagger(data=df_types, **metadata)
    ind = 0
    for label, bag in tb.bags:
        if stop:
            if ind == stop:
                break
        TasicNeuron(*bag, label=label, override=True)
        ind += 1
    config.write()
    config.write_python()


if __name__ == '__main__':
    main(stop=2)
