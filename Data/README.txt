Description of files contained in this data download:
1. 	genes_counts.csv: File containing read count values obtained from the RSEM algorithm for each gene (row) for each cell (column)
2.	genes_rpkm.csv: File containing the RPKM values obtained from the RSEM algorithm for each gene (row) for each cell (column)
3.	ercc_counts.csv: File containing read count values obtained from the RSEM algorithm for each external ERCC spike-in control (row) for each cell (column)
4.  cell_metadata.csv: File containing information about each cell profiled, including its nomenclature, Cre line of origin, dissection, date of collection and sequencing, and read mapping statistics
5.	cluster_metadata.csv: File containing information about each data-driven cluster, including its label, the corresponding label in Tasic et al. (Nat. Neuro, 2106), the primary cell class membership, and marker genes (including genes with widespread expression in the cluster, sparse expression in the cluster, and no expression in the cluster). 
6. 	cell_classification.csv: File containing information about the cluster membership of each cell, including whether the cell is a "core" (unambiguously assigned to a single cluster) or "transition" (shares membership between two clusters) cell, as well as its membership score (from 0-10) for each cluster (labeled f01 to f49). 

To generate the count and RPKM data, 100bp single-end reads were aligned using RSEM to the mm10 mouse genome build with the RefSeq annotation downloaded on 11 June 2013.
Raw fastq files are available at Gene Expression Omnibus, accession ID GSE71585