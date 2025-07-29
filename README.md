# Propagate go-terms using Python
**Propagate_go_terms.py** is a very simple script to propagate GO terms from GAF annotations or [DeepFRI](https://github.com/flatironinstitute/DeepFRI) csv files using the Gene Ontology hierarchy.

Propagation means adding parent go-terms (related, but usually less specific) to your query terms. 
deepFRI predictions also have scores - in this implementation the parent terms receive the **highes score** of their children. 

The main propagation engine is the get_all_parents() method from [goatools](https://github.com/tanghaibao/goatools). The current implementation uses only **is_a** [relationship](https://github.com/tanghaibao/goatools/blob/main/notebooks/parents_and_ancestors.ipynb).

## Main script supports three modes:

--mode gaf: propagate from a GAF file (optionally filtered by your protein list).

--mode pred: propagate from DeepFRI CSV predictions.

--mode both: run both in one invocation.

## Key flags:

--gaf: path to GOA GAF file.

--proteins: optional one-per-line protein ID list for GAF filtering.

--pred_pattern: glob pattern for DeepFRI CSVs.

--score_cutoff: float, threshold for DeepFRI scores.

--id_from_filename: when set, overrides the default naming from csv "protein" column by parsing the filename.

--obo: path to GO OBO file.

--ic: optional CSV with information content (go_term,IC).

--output_gaf / --output_pred: output CSV paths.

## Installation
```
git clone https://github.com/FilipSchymik/propagate_go_terms.git
cd propagate_go_terms

# create conda env
conda env create -f environment.yml
conda activate propagate_go_terms
```

## Usage
1. GAF only, all proteins:
```
python propagate_go_terms.py \
  --mode gaf \
  --gaf data/goa_uniprot_all.gaf \
  --obo data/go.obo \
  --output_gaf results/gaf_propagated.csv
```
2. DeepFRI only, load IDs from filenames:
```
python propagate_go_terms.py \
  --mode pred \
  --pred_pattern "predictions/*predictions.csv" \
  --id_from_filename \
  --obo data/go.obo \
  --output_pred results/deepfri_propagated.csv
```
3. Both modes with IC integration:
```
python propagate_go_terms.py \
  --mode both \
  --proteins data/af_best_4mln.txt \
  --gaf data/goa_uniprot_all.gaf \
  --pred_pattern "predictions/*_predictions.csv" \
  --obo data/go.obo \
  --ic data/IC_swissprot.csv \
  --output_gaf results/gaf.csv \
  --output_pred results/pred.csv
```
