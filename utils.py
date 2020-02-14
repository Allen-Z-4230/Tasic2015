# tree stuff
import anytree
import anytree.util
from anytree.importer import DictImporter
from anytree.exporter import DotExporter


def name_mapping(major_class, major_dissectoin, sub_class):  # New Names for hiearchical groups
    custom_name = major_class + major_dissectoin + sub_class + 'cell'
    return custom_name


def parse_binary(mydict, bin_str, leaf_name):
    """Parses binary node position codes into python dictionary

    mydict: dict, root dictionary to parse through
    bin_str: str, string of binary characters
    leaf_name: str, assigned name of the leaf node
    """

    bls = [int(x) for x in list(bin_str)]
    root = mydict
    for ind, node_val in enumerate(bls):
        if 'children' in root:
            pass
        else:
            root.update({'children': []})

        if node_val == 0:  # left child
            if len(root['children']) == 0:
                root['children'].append({})
                root = root['children'][0]
            else:
                root = root['children'][0]

        if node_val == 1:  # right child
            if len(root['children']) == 1:
                root['children'].append({})
                root = root['children'][1]
            else:
                root = root['children'][1]

        if ind == len(bls) - 1:  # check leaf node (last iteration)
            root['label'] = leaf_name


def get_node(tree, pos):
    """gets node based on position
       pos: int, breadth-first node position

    """
    return anytree.search.findall(tree, filter_=lambda node: node.pos == pos)[0]


def cluster_converter(row, tree=None):
    """Function for pandas apply. Converts cluster_id into its position on
    the dendrogram based on least common acestor.

    row: DataFrame, single row of data
    tree: AnyTree, tree containing the dendrogram
    """
    if row['coretype'] == 'Core':
        target_node = anytree.search.findall_by_attr(tree, row['primary'],
                                                     name='cluster_id')[0]

    elif row['coretype'] == 'Transition':  # transition clusters
        n1 = anytree.search.findall_by_attr(tree, row['primary'],
                                            name='cluster_id')[0]
        n2 = anytree.search.findall_by_attr(tree, row['secondary'],
                                            name='cluster_id')[0]
        target_node = anytree.util.commonancestors(n1, n2)[-1]

    return target_node.pos


def gene_merge(row, df=None, index='markers_present'):
    """Function for pandas apply. Converts cluster_id into its position on
    the dendrogram based on least common acestor.

    row: DataFrame, single row of data
    df: DataFrame, dataframe that store the gene data.
    index: 'markers_present' or 'gene_absent'
    """
    if row['coretype'] == 'Core':
        cell = df[df['cluster_id'] == row['primary']][index]
        if cell.isnull().all():
            gene_set = set()
        else:
            gene_set = set(cell.str.split(",").tolist()[0])

    elif row['coretype'] == 'Transition':
        cell1 = df[df['cluster_id'] == row['primary']][index]
        cell2 = df[df['cluster_id'] == row['secondary']][index]
        if cell1.isnull().all():
            gene_set1 = set()
        else:
            gene_set1 = set(cell1.str.split(",").tolist()[0])
        if cell2.isnull().all():
            gene_set2 = set()
        else:
            gene_set2 = set(cell2.str.split(",").tolist()[0])
        if index == 'markers_present':
            gene_set = set()
        elif index == 'genes_absent':
            gene_set = gene_set1.union(gene_set2)

    return gene_set


def diff(l1, l2, symmetric=True):
    """Returns differences of two iterables"""
    if symmetric:
        return list(set(l1) ^ set(l2))
    else:
        return list(set(l1) - set(l2))
