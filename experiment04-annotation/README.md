## Scripts
|Script | Comments|
|-----|---------|
|00_data_prep.ipynb| preparation of datasets|
|01_annot_wo_fine_tuning.py| annotation script without initial fine tuning (FAIL)|
|02_annot_with_fine_tuning.py| fine-tuning on reference and prediction of annotation on query dataset |
|03_annot_using_fine_tuned.py| prediction of annotation on query dataset using already fine-tuned model|
|04_prediction_results.ipynb| plotting results of fine-tune and predict (see 02) |
|05_ped_on_fine_tuned_results.ipynb| plotting results of prediction using already fine-tuned model (see 03)|
|06_scvi_pools_dataset.py| subsetting sarcoma dataset (pool7 and pool8) |
|07_subset_heroes_dataset.ipynb| HVGs of sarcoma dataset (pool7 and pool8) for model |