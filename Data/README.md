#  `Data/`

## Processed mRNA and protein abundances

The processed mRNA and protein abundances used in our model fits are contained in the TSV file `mean_se_rna_protein.tsv`.

## Phylogenetic trees

This directory contains three phylogenetic trees: 

- `vertlife_mcc_dna_only_node_date.tre` (Primary analysis, Newick format, obtained from [VertLife MCC](https://vertlife.org/data/mammals/))
- `vertlife_mcc_dna_only_fossil_bd.tre` (Newick format, obtained from [VertLife MCC](https://vertlife.org/data/mammals/))
- `tree_11sp_noGpig.nex` (Nexus format, kindly provided by Dr. Gunter Wagner)

We note that `tree_11sp_noGpig.nex` places the species horse in the order Carnivora, which does not have much support in the literature. However, we found that using any of these trees generally results in the same conclusions. 

## Other data

Other directories contain `README.md` explaining data sources and processing. 