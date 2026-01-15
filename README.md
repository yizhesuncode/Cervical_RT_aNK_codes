Cervical_RT_aNK_codes

This repository contains all R and Python scripts used to investigate the treatment-specific roles of tumor-associated adaptive NK (aNK) cells in cervical cancer, with a particular focus on radiotherapy-associated immune modulation.

🔬 Project overview

This project aims to characterize:

Treatment-specific functional adaptations of tumor-infiltrating adaptive NK cells

Radiotherapy (RT)–associated transcriptional reprogramming

Immune memory–like features of aNK cells

Gene regulatory networks underlying treatment response

Analyses are based on single-cell RNA-seq data and in silico perturbation modeling.

📂 Repository structure
Cervical_RT_aNK_codes/
│
├── All_R_code_for_cervical_RT.R      # Main R analysis pipeline
├── add.flag.R                       # Custom visualization utilities
├── Celloracle_GRN_Perturbation.py   # GRN perturbation (CellOracle)
├── README.md
├── .gitignore
└── *.Rproj

🧬 Main analyses
R-based analyses

Data preprocessing & QC

Cell clustering and annotation

Differential expression analysis

Adaptive NK cell subset characterization

Visualization (UMAP, heatmaps, dotplots, etc.)

Python-based analyses

Gene regulatory network construction

In silico TF perturbation modeling

CellOracle-based vector field simulations

Developmental trajectory analysis

🛠 Requirements
R packages

Seurat

ggplot2

dplyr

ComplexHeatmap

SingleR

patchwork

(see script headers for full list)

Python packages

celloracle

scanpy

numpy

pandas

matplotlib

▶ How to run
R pipeline
source("All_R_code_for_cervical_RT.R")

GRN perturbation (Python)
python Celloracle_GRN_Perturbation.py


⚠ Please update file paths according to your local environment.

📊 Output

The pipeline generates:

UMAP visualizations

Heatmaps of marker genes

Differential gene expression tables

Perturbation score distributions

Regulatory network plots

📖 Citation

If you use this code, please cite:

Sun Y. et al.
Treatment-specific regulation of tumor-associated adaptive NK cells in cervical cancer
Manuscript in preparation.

👤 Author

Yizhe Sun
Karolinska Institutet
Email: (optional)

⚠ Notes

Raw sequencing data are not included

Large intermediate files are excluded via .gitignore

Scripts are provided for research and reproducibility purposes

🤝 Contact

For questions, suggestions, or collaboration, feel free to:

Open an issue

Contact the author directly