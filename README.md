# Tasic 2015 Neuron Type Ontology

**Primary Goal**: To convert transcriptomic cell types in the Tasic 2015 Paper into a neuron type ontology. We are adopting a new framework used by NIF which models the neurons as discrete classes with phenotypes as properties.

Tasic 2015 paper:

- 49 transcriptomic cell types, including 23 GABAergic, 19 glutamatergic and 7 non-neuronal types in the mouse **visual cortex**
- some of our transcriptomic cell types displayed specific and differential electrophysiological and axon projection properties
- **Single Cell RNA-seq**, which offers single-cell resolutions of gene expression information, allowing for cell-type related analyses

  - isolating one individual cells from a population of cells (tissue sample).
  - isolate the RNA, convert to cDNA with reverse transcriptase.
  - Amplification with polymerase chain reaction (PCR)
  - Sequencing with Next-Generation Sequencing (NGS) or similar techniques (in this case, authors used the SMARTer PCR Protocol)
  - Transcriptome (genomic features) extracted from template matching or _de novo_ assembly

- Used transgenic mice expressive **Cre-Recombinase (cre-lines)**, which can regulate the expression of certain genes at specific locations and times. This can aid studying differential expression as well as genetic influences in development.

  - Cre-drivers: usually named from the gene(s) at their promotor sites

    - [NIH Cre-Driver Network](www.credrivermice.org/database)
    - [Jackson Laboratory Cre Database](https://www.jax.org/research-and-faculty/resources/cre-repository/characterized-cre-lines-jax-cre-resource)
    - [Allen Institute Transgenic Database](https://connectivity.brain-map.org/transgenic)

  - Cre-reporters: molecular markers that can verify that certain genes have been expressed.

- Clustering methods:

  - Used the intersection of **PCA** and iterative weighted gene coexpression network analysis (**WGCNA**)

    - PCA: iteratively finds orthogonal components of the feature space which are linear combinations of genes.
    - WGCNA: A form of hierarchical clustering which constructs an adjacency matrix as well as proximity measure between them to perform branch cutting. Modules resulting from this network branch cut are then characterized by their eigengene.

  - Validation of cluster membership using random forest.

- Clustering uncovered common broad neuronal categories (matched using cre-layer enrichment, previously known marker genes, and expression profiles) such as excitatory (glutamatergic) vs inhibitory (GABAergic).

- In total, 49 transcriptomic types have been identified, each with different expression profiles and marker genes. These are useful for us because they can be cross-referenced with cre-line specific layers. For any given transgenic mice, this would tell us what markers tend to be expressed in each layer.

- In summary, we are capturing the following phenotypes:

**Definitive phenotype**:

```
- Cre line [NCBI Taxon] [Allen Brain Institute]
- V1
- Layer of dissection [UBERON]
```

**Inferred gene markers phenotype**

```
- Cluster ID [ILX]
- Genes present [NCBIGene]
- Genes Absent [NCBIGene]
```

## For the Future:

- Continuous Phenotypes (gene counts, expression profiles, etc.)
- Phenotype Inheritance in subclasses?

--------------------------------------------------------------------------------

## Interpretations & Assumptions:

- Layer VI was separated into two sub-compartments a and b, these are assumed to be computationally significant, though we are not certain what metric was used to divide the layer.

- Cells classified with a "transition" core type were ambiguous.

  - Markers present was discarded for these cells
  - Markers absent was the union of absent markers in their categories

- "NEG" (knock-out) cre-lines are modeled as complement of (hasProteinExpr, someDriveLine)

- Combination Cre Lines & Mismatches:

  - PvalbD-Slc32a1: Pvalb-T2A-Dre AND Slc32a1-IRES-Cre

  - PvalbF-Gad2: Pvalb-2A-Flpo AND Gad2-IRES-Cre

  - Tac2: Tac2-IRES2-Cre

  - Nkx2-1: Nkx2-1-CreERT2

  - Ctgf: Ctgf-2A-dgCre

  --------------------------------------------------------------------------------

## Sources of Data
