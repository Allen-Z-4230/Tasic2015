# Project-Level Agenda:

- Document sources of data

## Code To-Dos:

- Cre line phenotypes have no JAX numbers
- Fix cre lines (Neg & AND phenotypes)

  - Correct by actually modifying the data (correct across all cells)
  - New function for combination cre lines

- join all unmapped genes into a list for debugging

- Implement threading/multiprocessing OR run them on the server

  - setup ssh

- Issue: SubClass Violation when compiling Layers

- Read the codebase when possible

## Project Journal

February 12, 2019:

Project rebooted after a long hiatus at the end of fall quarter. Working on:

1. Making cre-line conversion code
2. Exporting missing genes
3. Documenting progress

Meeting | February 13, 2019:

Questions:

Conceptual:

1. Why use cre lines at all if we're only trying to get a complete picture of cell types in the visual cortex? Why is it necessary to use transgenic lines to characterize what cell types (and markers) are expressed in each layer?

2. Directionality of inference. Do we do (cluster x and marker y) -> (usually found in layer a of cre line b) or (layer a of cre-line b) -> (cell usually exhibit maker y)?

Code:

1. Do we need the cre-line name for the phenotypes or just the stock number / ID?

  - Stock number, include leading zero, interpret as string.
  - names as labels, abbreviation as a synonym

2. New namespace for cre_lines? -- no, use pre-build dictionary

3. No source number for specific cre-lines

  a. Jackson -> JAX

  b. ctgr -> T2A, manualy entry of stock number & source

4. gitignore

5. rename Data to resources

6. code reorg & cleanup

# Resources

## [OWL Reference](https://www.w3.org/TR/owl-ref)

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
  Markram 2015: [PMID:26451489](https://www.ncbi.nlm.nih.gov/pubmed/26451489)<br>
  Huang 2017: [PMID:27053207](https://www.ncbi.nlm.nih.gov/pubmed/27053207)<br>
  Tasic 2015:

  - <http://casestudies.brain-map.org/celltax>
  - <https://www.nature.com/articles/nn.4216>

5. Codebase:

  - Bolser Lewis
  - OntTerm triplesimple
  - nistd/core.py
  - curation.py = triplesExport
