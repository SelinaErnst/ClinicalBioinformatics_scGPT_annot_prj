# scGPT

## Scripts
|Script | Comments|
|-----|---------|
|00_import.ipynb| imports required to run annotation tutorial |
|01_ms_dataset.ipynb| ms dataset used in annotation tutorial |
|03_annotation.py| script following the annotation tutorial |
|04_job_efficiency.ipynb| comparison between completed job runs |
|05_prediction_on_fine_tuned_model.py| annotation tutorial using already fine-tuned model |
|06_prediction_plot.ipynb| resulting plots (umap and confusion matrix) |
|07_prediction_on_human_model.py| annotation without fine-tuning |
|08_compare_predictions.ipynb| comparison of confusion matrices |

## Issues

### Environment setup: Dependency conflict

For the most part of the analysis I've been using `scgpt_env` which however was created by me trying to solve dependency issues as they came up and is therefore hard to replicate. Here I'm trying to figure out how I can get a working conda environment using environment.yml files.

The issue is that `scgpt 0.2.4` (latest public release) declares: `scvi-tools<1.0,>=0.16.0`

However, the older `scvi-tools 0.20.3` is incompatible with newer AnnData versions because it relies on an internal AnnData API (SparseDataset) that changed in newer releases. 

There are therefore two options:

- Use the older dependency stack with `scvi-tools 0.20.3` and an older AnnData version: `environment_old.yml`
    ```
    scvi-tools 0.20.3
    anndata 0.9.2
    ```
    This combination works with scGPT 0.2.4. However, both the older AnnData 0.9.2 and the newer AnnData 0.11.4 can have problems reading newer `.h5ad` files. For example:
    ```python
    IORegistryError: No read method registered for IOSpec(
        encoding_type='null',
        encoding_version='0.1.0'
    )
    ```
- Use newer AnndData and and scvi-tools: `environment.yml` 
    ```
    anndata 0.11.4
    scvi-tools 1.3.3
    ```
    Although scGPT 0.2.4 declares scvi-tools<1.0, I tested that it works with scvi-tools 1.3.3. Because pip's dependency resolver rejects this combination based on scGPT's metadata, scGPT need to be installed separately without resolving its dependencies:
    ```bash
    pip install --no-deps scgpt==0.2.4
    ```

