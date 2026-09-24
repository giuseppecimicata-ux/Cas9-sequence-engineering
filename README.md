# Cas9 Sequence Engineering

Python tools for sequence selection, processing, and construct design in Cas9 protein engineering workflows.

## Overview

This repository contains Python workflows developed during a Cas9 protein-engineering project for the selection, processing, and preparation of designed protein variants.

The scripts connect computational variant-selection outputs with experimentally usable sequence designs by identifying selected candidates, retrieving and comparing protein sequences, reconstructing domain context, and generating corresponding DNA sequences.

The project was exploratory and did not result in a publication. The repository therefore focuses on the computational sequence-processing and construct-design tools developed during the work.

## Workflow

The computational workflow connects variant-selection outputs to construct-ready sequences:

```text
Variant-selection matrix
        ↓
Candidate identification
        ↓
Selected variant IDs
        ↓
Protein-sequence extraction
        ↓
Comparison with reference sequence
        ↓
Sequence processing and domain reconstruction
        ↓
DNA construct generation
```

## Repository contents

### `analyze_variant_selection.py`

Analyzes variant-selection matrices, visualizes the selection landscape, identifies selected designs, and exports their sequence identifiers for downstream processing.

### `extract_and_compare_variants.py`

Retrieves selected protein variants from FASTA files, includes the corresponding reference sequence, and calculates sequence differences relative to the reference.

### `prepare_cas9_variant_sequences.py`

Processes selected Cas9 variant sequences for downstream construct design, including sequence cleanup and optional restoration of required flanking regions.

### `generate_cas9_variant_dna.py`

Generates DNA sequences for designed protein variants while preserving the original codons at unchanged positions and assigning predefined codons at engineered residues.

## Input data

Experimental and project-specific input data are not included in this repository.

Users can provide paths to their own input files in the configuration section of each script. Expected input formats are documented within the individual scripts.

## Requirements

The workflows use standard Python scientific and sequence-processing libraries, including:

* pandas
* NumPy
* Matplotlib
* openpyxl

Exact dependencies are listed in `requirements.txt`.

## Notes

These scripts were developed as research tools during an exploratory Cas9 protein-engineering project. They have been refactored for readability and reuse while preserving the original analysis and sequence-design logic.

## Author

**Giuseppe Cimicata, PhD**
University of California, San Francisco

Research areas: protein engineering, structural biology, quantitative biophysics, genome-editing systems, and computational protein design.
