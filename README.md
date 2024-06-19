# GTEx analysis

To reproduce analyses, first download the data from GTEx website and GENCODE or from the drive folder

Then create the conda environment with:
```bash
conda env create -f gtex.yml
```

The entire analysis is then reproduced by running scripts in the `analysis` folder:
```bash
conda activate gtex
python 0_load_data_subset.py
Rscript 1_yarn_normalize.R
python 2_specificity_analysis.py
```

Given the simplicity of the pipeline we don't use tools such as snakemake here.

