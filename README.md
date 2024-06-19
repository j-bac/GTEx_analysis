# GTEx analysis

To reproduce analyses, first download the data from GTEx website and GENCODE or from the drive folder

Then create the conda environment with:
```bash
conda env create -f gtex.yml
```

The analysis is then reproduced by running two notebooks and one script in the `analysis` folder. If preferred, notebooks can be executed as scripts using nbconvert. The full pipeline is then performed by:
```bash
conda activate gtex
jupyter nbconvert --execute 0_load_data_subset.ipynb
Rscript 1_yarn_normalize.R
jupyter nbconvert --execute 2_specificity_analysis.ipynb
```
