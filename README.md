# Proteomic Analysis of Polyester Degradation in Carnivorous Plants

This repository contains data analysis notebooks and scripts used during the Master's Thesis *"Proteomic Analysis of Polyester Degradation in Carnivorous Plants"*.

The study investigates the digestive proteome of carnivorous pitcher plants (*Nepenthes* and *Sarracenia*) in response to different feeding conditions, focusing on identifying proteins involved in the hydrolysis of synthetic polyesters PET and PBAT. The repository includes R notebooks and scripts used for data processing, differential expression analysis, and homology search.

## Contents

- **`notebooks/`**: R notebooks that document the analysis workflow of FragPipe results:
    - `thesis_analysis_HELA.rmd`: Analysis of HeLa standards.
    - `thesis_analysis_nepenthes.rmd`: Analysis of *Nepenthes* samples.
    - `thesis_analysis_sarracenia.rmd`: Analysis of *Sarracenia* samples.

- **`scripts/`**: Scripts used for specific tasks, including:
    - `extract_sequences.py`: Provides functions to extract sequences from a FASTA or FASTQ file based on a list of sequence IDs.
    - `plot_mzml.py`: Generates 2D overview plots of mass spectrometry data from mzML files.
    - `reformat_fasta_header.ipynb`: Annotates FASTA file headers based on a supplied metadata file.
    - `write_fp_manifest.R`: Reads all mzML files in a directory and generates a manifest for FragPipe.

- **`README.md`**: This file, providing an overview of the repository and instructions.

## Requirements

To reproduce the analyses, you will need the following software and R packages:

- **R version 4.2.1**
- [**RStudio**](https://posit.co/download/rstudio-desktop/)
- **Required R packages**:
    - `here`
    - `tidyverse`
    - `magrittr`
    - `limma`
    - `ComplexHeatmap`
    - `Biobase`
    - `pheatmap`
    - `EnhancedVolcano`
    - `ggfortify`
    - `cowplot`
- **BLAST**

## Visualization and Docking Tools

The following tools were used for visualization and docking studies in this research:

- [**PyMol**](https://github.com/schrodinger/pymol-open-source): Used for molecular visualization. 
- **PyMol Packages**:
    - [**DockingPie**](https://github.com/paiardin/DockingPie): Employed for docking simulations.
- [**PLIP**](https://plip-tool.biotec.tu-dresden.de/plip-web/plip/index): Utilized for analyzing protein-ligand interactions.

**Note**: The specific steps and scripts related to the use of these tools are not included in this repository. The focus of this repository is on the data processing, differential expression analysis, and other R-based analyses performed during the study.

## License

This repository is licensed under the GNU General Public License v3.0. You may use, distribute, and modify this code under the terms of the GPL-3.0 license. For more details, see the [LICENSE](./LICENSE) file.

## Author
Sebastian Gritsch 201810[at]fhwn.ac.at
