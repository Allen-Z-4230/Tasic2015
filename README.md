Interpretations:
- Cells classified with a "transition" coretype were ambiguous.
    - Markers present was discarded for these cells
    - Markers absent was the union of absent markers in their categories

Assumptions

### Useful Links
1. The neurons branch  
https://github.com/SciCrunch/NIF-Ontology/tree/neurons/ttl  
https://github.com/SciCrunch/NIF-Ontology/blob/neurons/ttl/neuron-development.ttl  
https://github.com/SciCrunch/NIF-Ontology/blob/neurons/ttl/bridge/neuron-bridge.ttl  
https://github.com/SciCrunch/NIF-Ontology/blob/neurons/ttl/phenotype-core.ttl  
https://github.com/SciCrunch/NIF-Ontology/blob/neurons/ttl/phenotypes.ttl  
https://github.com/SciCrunch/NIF-Ontology/blob/neurons/ttl/generated/neurons/common-usage-types.ttl  
https://github.com/SciCrunch/NIF-Ontology/blob/neurons/ttl/generated/neurons/allen-cell-types.ttl  
https://github.com/SciCrunch/NIF-Ontology/blob/neurons/ttl/generated/neurons/huang-2017.ttl  
https://github.com/SciCrunch/NIF-Ontology/blob/neurons/ttl/generated/neurons/markram-2015.ttl  

2. Mappings between hbp-cell and the new interlex identifiers.  
https://github.com/SciCrunch/NIF-Ontology/blob/neurons/ttl/generated/NIF-Neuron-HBP-cell-import.ttl  
https://github.com/tgbugs/pyontutils/blob/master/hbp_cell_conv.csv  

3. Python code that generates the new neuron ttl files  
[neurondm/lang.py](https://github.com/tgbugs/pyontutils/blob/master/neurondm/neurondm/lang.py)  
[neurondm/core.py](https://github.com/tgbugs/pyontutils/blob/master/neurondm/neurondm/core.py)  
[neurondm/build.py](https://github.com/tgbugs/pyontutils/blob/master/neurondm/neurondm/build.py)  
[neurondm/models](https://github.com/tgbugs/pyontutils/tree/master/neurondm/neurondm/models)  
[nlxeol/lift_neuron_triples.py](https://github.com/tgbugs/nlxeol/blob/master/lift_neuron_triples.py)  
[nifstd/resources/26451489 table 1.csv](https://github.com/tgbugs/pyontutils/blob/master/nifstd/resources/26451489%20table%201.csv)  

4. Papers  
Markram 2015 [PMID:26451489](https://www.ncbi.nlm.nih.gov/pubmed/26451489)  
Huang 2017 [PMID:27053207](https://www.ncbi.nlm.nih.gov/pubmed/27053207)  