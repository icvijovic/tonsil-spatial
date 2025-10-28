# tonsil-spatial

Code associated with "Dynamics of local B cell migration during affinity maturation in the human tonsil"

snakemake_workflow/: Snakemake pipeline to process VDJ sequences, performing lineage and germline calling. Output can be found on Zenodo (details below) and is used as input to the figure generation files.

Figure_Code_Main_Text.ipynb: Code to analyze processed data and produce figures in main text.
Figure_Code_SI.ipynb: Code to analyze processed data and produce figures in supplementary notes.
data_analysis_utils.py: File containing helper functions used by the figure code generation files.

ZENODO:
Preprocessed_Data/: Contains some intermediate data with longer analysis times, mostly for spatial distance and lineage threshold analysis.
Lineage_Seqs/: Lineage sequences and trees used in phylogenetic analysis.
Clustering_Analysis/: Additional lineage clustering data used in lineage threshold analysis.
All other files are processed files from the snakemake pipeline (run on data from https://zenodo.org/record/7326538).
