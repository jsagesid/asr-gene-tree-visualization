**Gene Tree Visualization with Ancestral Protein Properties**

R scripts for visualizing gene trees with ancestral sequence reconstructions and protein property evolution.

These tools integrate phylogenetic information, ancestral sequence reconstruction, and biochemical property calculations to examine how protein features change through evolutionary time.

**Overview**

This repository contains scripts for visualizing evolutionary changes in protein properties across a time-calibrated gene tree.

The workflow combines:

- Time-calibrated phylogenetic trees (e.g., generated with MCMCtree)

- Ancestral sequence reconstruction (CodeML / PAML)

- Protein property calculations for both extant and ancestral sequences

- Trait data visualization for species at the tips of the tree

Protein properties are calculated from the aligned sequences and mapped onto the phylogeny. The scripts then quantify and visualize changes in protein properties between ancestral nodes and descendant branches.

Example protein properties include:

- Net charge

- Hydrophobicity

The visualization can optionally incorporate phenotypic or ecological trait data, such as milk fat percentage.

**Features**

The scripts allow users to:

- Visualize a time-calibrated gene tree

- Map ancestral protein sequences onto internal nodes

- Calculate protein properties for extant and reconstructed sequences

- Quantify changes in protein properties between nodes

- Highlight branches with significant biochemical shifts

- Overlay trait data for extant taxa

This allows researchers to explore how protein evolution relates to phenotypic or ecological traits.

**Input Files**

The pipeline requires the following inputs.

1. Phylogenetic Tree

A time-calibrated tree file.

- Supported formats:

- .nwk

- .txt

Example sources:

- MCMCtree

- BEAST

- other phylogenetic inference software

2. Multiple Sequence Alignment

A protein multiple sequence alignment (MSA) in FASTA format.

- Requirements:

- Amino acid sequences

- Sequence names must match those used in the tree

Example:

>Homo_sapiens
MTEYKLVVVGAGGVGKSALTIQLIQ
>Phoca_vitulina
MTEYKLVVVGAGGVGKSALTIQLIQ

3. Protein Property Table

A table containing calculated protein properties for each sequence (e.g. net charge and hydrophobicity).

- This table can be generated using the protein property calculation script included in this repository.

4. Trait Data (Optional)

Optional trait data for extant species.

- Trait data can be used to visualize ecological or physiological variation alongside protein evolution.

**Protein Property Calculation**

This repository also includes a script for calculating protein biochemical properties directly from MSAs.

Currently implemented properties include:

- Net charge

- Hydrophobicity

These values are calculated for:

- extant sequences

- reconstructed ancestral sequences

The resulting table can be used directly as input for the visualization script.

**Output**

The scripts generate visualizations of a gene tree annotated with evolutionary protein properties, including:

- ancestral sequences mapped to internal nodes

- changes in protein properties along branches

- optional trait values displayed at tip taxa

This allows users to identify major evolutionary shifts in protein chemistry across the phylogeny.

**Example Applications**

This approach is useful for investigating:

- molecular adaptations

- functional shifts in protein evolution

- correlations between protein properties and ecological traits

- evolutionary changes in proteins, metabolic enzymes, or other functional genes
