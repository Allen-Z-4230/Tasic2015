# Tasic 2015 Neuron Type Ontology

## Project-level:

- Generation of Cell names
- names for higher subclasses
- penotypic & Genotypic -- Check threshold

--------------------------------------------------------------------------------

- Start building from individual neurons
- Naming scheme: Tasic 15, cluster [cluster_number]
- Function that changes names flexibly
- Templates & References

  - Bolser Lewis
  - OntTerm triplesimple
  - nistd/core.py
  - curation.py = triplesExport

## Code-level:

- cre line code should now be fixed, build it w/ NEG phenotypes
- output unmapped gene list for debugging -

Definitive phenotype:

```
- Cre line [NCBI Taxon] [Allen Brain Institute]
- V1
- Layer of dissection [UBERON]
```

Inferred gene markers phenotype

```
- Cluster ID [ILX]
- Genes present [NCBIGene]
- Genes Absent [NCBIGene]
```

## For the Future:

- Continuous Phenotypes (gene counts, etc.)
- Phenotype Inheritance in subclasses?

--------------------------------------------------------------------------------

## Interpretations & Assumptions:

- Layer VI was separated into two sub-compartments a and b, these are assumed to be computationally significant, though we are not certain what metric was used to divide the layer.

- Cells classified with a "transition" core type were ambiguous.

  - Markers present was discarded for these cells
  - Markers absent was the union of absent markers in their categories

## Useful Links

1. The neurons branch<br>
  <https://github.com/SciCrunch/NIF-Ontology/tree/neurons/ttl><br>
  <https://github.com/SciCrunch/NIF-Ontology/blob/neurons/ttl/neuron-development.ttl><br>
  <https://github.com/SciCrunch/NIF-Ontology/blob/neurons/ttl/bridge/neuron-bridge.ttl><br>
  <https://github.com/SciCrunch/NIF-Ontology/blob/neurons/ttl/phenotype-core.ttl><br>
  <https://github.com/SciCrunch/NIF-Ontology/blob/neurons/ttl/phenotypes.ttl><br>
  <https://github.com/SciCrunch/NIF-Ontology/blob/neurons/ttl/generated/neurons/common-usage-types.ttl><br>
  <https://github.com/SciCrunch/NIF-Ontology/blob/neurons/ttl/generated/neurons/allen-cell-types.ttl><br>
  <https://github.com/SciCrunch/NIF-Ontology/blob/neurons/ttl/generated/neurons/huang-2017.ttl><br>
  <https://github.com/SciCrunch/NIF-Ontology/blob/neurons/ttl/generated/neurons/markram-2015.ttl>

2. Mappings between hbp-cell and the new interlex identifiers.<br>
  <https://github.com/SciCrunch/NIF-Ontology/blob/neurons/ttl/generated/NIF-Neuron-HBP-cell-import.ttl><br>
  <https://github.com/tgbugs/pyontutils/blob/master/hbp_cell_conv.csv>

3. Python code that generates the new neuron ttl files<br>
  [neurondm/lang.py](https://github.com/tgbugs/pyontutils/blob/master/neurondm/neurondm/lang.py)<br>
  [neurondm/core.py](https://github.com/tgbugs/pyontutils/blob/master/neurondm/neurondm/core.py)<br>
  [neurondm/build.py](https://github.com/tgbugs/pyontutils/blob/master/neurondm/neurondm/build.py)<br>
  [neurondm/models](https://github.com/tgbugs/pyontutils/tree/master/neurondm/neurondm/models)<br>
  [nlxeol/lift_neuron_triples.py](https://github.com/tgbugs/nlxeol/blob/master/lift_neuron_triples.py)<br>
  [nifstd/resources/26451489 table 1.csv](https://github.com/tgbugs/pyontutils/blob/master/nifstd/resources/26451489%20table%201.csv)

4. Papers<br>
  Markram 2015 [PMID:26451489](https://www.ncbi.nlm.nih.gov/pubmed/26451489)<br>
  Huang 2017 [PMID:27053207](https://www.ncbi.nlm.nih.gov/pubmed/27053207)<br>
  Tasic 2015:

  - <http://casestudies.brain-map.org/celltax>
  - <https://www.nature.com/articles/nn.4216>

--------------------------------------------------------------------------------

Meeting (10/10)

1. Error when importing genes from Huang2017 - fixed
2. Some cre-lines & genes are not able to be mapped:

  - Neg Cre Lines: complement of (hasProteinExpression, someDriveLine)

3. somehow the layers are not serializing - fixed

  - NVM, they are serializing correctly

4. Have a good way to document the source of the data (csv, xlsx, etc)

--------------------------------------------------------------------------------
