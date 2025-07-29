import pandas as pd
import glob
import time
import sys
from pathlib import Path
from goatools.obo_parser import GODag
from tqdm.auto import tqdm


def load_protein_ids(path: str) -> set:
    """
    Load a one-column text file of protein IDs (one per line).
    """
    with open(path) as fh:
        return {line.strip() for line in fh if line.strip()}


def filter_gaf_by_proteins(
    gaf_path: str,
    proteins: set = None,
    chunksize: int = 2_000_000,
) -> pd.DataFrame:
    """
    Read a GAF file in chunks. If proteins is provided, filter to those proteins;
    otherwise include all entries. Always show a progress bar and timing.
    Returns DataFrame with columns: ['DB_Object_ID','GO_ID','Aspect'].
    """
    start = time.time()
    cols = ['DB','DB_Object_ID','GO_ID','Aspect']
    usecols = [0, 1, 4, 8]
    reader = pd.read_csv(
        gaf_path,
        sep='\t', comment='!', header=None,
        usecols=usecols,
        names=cols,
        dtype={'DB_Object_ID': 'category', 'GO_ID': 'category', 'Aspect': 'category'},
        chunksize=chunksize
    )
    collected = []
    total_rows = 0
    for chunk in tqdm(reader, desc="Filtering GAF"):
        if proteins is not None:
            chunk = chunk[chunk['DB_Object_ID'].isin(proteins)]
        if not chunk.empty:
            collected.append(chunk[['DB_Object_ID','GO_ID','Aspect']])
            total_rows += len(chunk)
    df = pd.concat(collected, ignore_index=True) if collected else pd.DataFrame(columns=['DB_Object_ID','GO_ID','Aspect'])
    print(f"Filtered GAF: {total_rows} rows in {time.time()-start:.2f}s")
    return df


def load_deepfri_predictions(
    pattern: str,
    score_cutoff: float = None,
    id_from_filename: bool = False
) -> pd.DataFrame:
    """
    Load all CSVs matching `pattern` for deepFRI predictions.
    Skips lines starting with '#'. Keeps original term name and aspect.
    If id_from_filename=True, overrides the 'protein' column using the filename.
    Optionally filter by score_cutoff. Always show timing.
    Returns DataFrame with ['protein','GO_ID','score','orig_name','aspect'].
    """
    start = time.time()
    files = glob.glob(pattern)
    dfs = []
    for f in files:
        df = pd.read_csv(f, comment='#')
        expected = ['Protein','GO_term/EC_number','Score','GO_term/EC_number name']
        missing = [c for c in expected if c not in df.columns]
        if missing:
            raise KeyError(f"Missing columns in {f}: {missing}")
        df = df.rename(columns={
            'Protein': 'protein',
            'GO_term/EC_number': 'GO_ID',
            'Score': 'score',
            'GO_term/EC_number name': 'orig_name'
        })[['protein','GO_ID','score','orig_name']]
        if df.empty:
            continue
        df['aspect'] = Path(f).stem.split('_')[-2]
        if id_from_filename:
            actual_prot = Path(f).stem.split('_')[0]
            df['protein'] = actual_prot
        if score_cutoff is not None:
            df = df[df['score'] >= score_cutoff]
        if not df.empty:
            dfs.append(df)
    preds = pd.concat(dfs, ignore_index=True) if dfs else pd.DataFrame(
        columns=['protein','GO_ID','score','orig_name','aspect']
    )
    print(f"Loaded deepFRI: {len(preds)} rows in {time.time()-start:.2f}s")
    return preds


def build_ancestor_map(
    godag: GODag,
    all_terms: set,
) -> pd.DataFrame:
    """
    Build ancestor mapping DataFrame: columns ['GO_ID','ancestor','aspect','term_name'].
    Always show progress bar and timing.
    """
    start = time.time()
    records = []
    ns_map = {'biological_process':'BP','molecular_function':'MF','cellular_component':'CC'}
    for go in tqdm(all_terms, desc="Building ancestor map"):
        node = godag.get(go)
        if node is None:
            continue
        ancestors = node.get_all_parents() | {go}
        for anc in ancestors:
            anc_node = godag[anc]
            aspect = ns_map.get(anc_node.namespace)
            records.append((go, anc, aspect, anc_node.name))
    df = pd.DataFrame.from_records(records, columns=['GO_ID','ancestor','aspect','term_name'])
    print(f"Built ancestor map: {len(records)} entries in {time.time()-start:.2f}s")
    return df


def propagate_gaf(
    filtered_gaf: pd.DataFrame,
    ancestor_map: pd.DataFrame,
    ic_df: pd.DataFrame = None
) -> pd.DataFrame:
    """
    Propagate GAF terms: merge with ancestor_map, dedupe, add optional IC and term_name.
    Always show timing.
    Returns DataFrame with ['protein','go_term','go_name','aspect'] plus 'IC' if provided.
    """
    start = time.time()
    df = filtered_gaf.rename(columns={'DB_Object_ID':'protein','GO_ID':'go_term'})
    merged = df.merge(ancestor_map, left_on='go_term', right_on='GO_ID', how='inner')
    out = merged[['protein','ancestor','aspect','term_name']]
    out = out.rename(columns={'ancestor':'go_term','term_name':'go_name'})
    out = out.drop_duplicates(['protein','go_term','aspect'])
    if ic_df is not None:
        out = out.merge(ic_df.rename(columns={'go_term':'go_term','IC':'IC'}), on='go_term', how='left')
    print(f"Propagated GAF: {len(out)} rows in {time.time()-start:.2f}s")
    return out


def propagate_deepfri(
    preds: pd.DataFrame,
    ancestor_map: pd.DataFrame,
    ic_df: pd.DataFrame = None
) -> pd.DataFrame:
    """
    Propagate deepFRI predictions: merge with ancestor_map, aggregate by max score,
    add ancestor term_name and optional IC.
    Always show timing.
    Returns DataFrame with ['protein','go_term','go_name','aspect','score'] plus 'IC' if provided.
    """
    start = time.time()
    merged = pd.merge(preds, ancestor_map, on='GO_ID', how='inner', suffixes=('_pred',''))
    merged = merged.rename(columns={'ancestor':'go_term','term_name':'go_name','score':'child_score'})
    merged = merged[['protein','go_term','go_name','aspect','child_score']]
    agg = merged.groupby(['protein','go_term','go_name','aspect'], as_index=False)
    agg = agg.child_score.max().rename(columns={'child_score':'score'})
    if ic_df is not None:
        agg = agg.merge(ic_df.rename(columns={'go_term':'go_term','IC':'IC'}), on='go_term', how='left')
    print(f"Propagated deepFRI: {len(agg)} rows in {time.time()-start:.2f}s")
    return agg


if __name__ == '__main__':
    import argparse

    parser = argparse.ArgumentParser("Propagate GO terms for GAF or deepFRI")
    parser.add_argument('--mode', choices=['gaf','pred','both'], default='both', help='Which data to propagate')
    parser.add_argument('--proteins', help='Optional protein list (one per line)')
    parser.add_argument('--gaf', help='GAF file path (required for gaf or both)')
    parser.add_argument('--pred_pattern', help='Glob for deepFRI CSVs (required for pred or both)')
    parser.add_argument('--score_cutoff', type=float, help='Optional score cutoff')
    parser.add_argument('--id_from_filename', action='store_true', help='Use protein ID from filename instead of protein column')
    parser.add_argument('--obo', required=True, help='GO OBO file path')
    parser.add_argument('--ic', help='Optional IC CSV with [go_term,IC]')
    parser.add_argument('--output_gaf', help='Output CSV for propagated GAF')
    parser.add_argument('--output_pred', help='Output CSV for propagated deepFRI')
    args = parser.parse_args()

    # Validate OBO
    if not Path(args.obo).is_file():
        parser.error(f"OBO file not found: {args.obo}")

    overall_start = time.time()

    # GAF mode
    if args.mode in ('gaf','both'):
        if not args.gaf or not args.output_gaf:
            parser.error('GAF mode needs --gaf and --output_gaf')
        proteins = load_protein_ids(args.proteins) if args.proteins else None
        filtered_gaf = filter_gaf_by_proteins(args.gaf, proteins)
        gaf_terms = set(filtered_gaf['GO_ID'])
    else:
        gaf_terms = set()

    # deepFRI mode
    if args.mode in ('pred','both'):
        if not args.pred_pattern or not args.output_pred:
            parser.error('Pred mode needs --pred_pattern and --output_pred')
        preds = load_deepfri_predictions(args.pred_pattern, args.score_cutoff, args.id_from_filename)
        pred_terms = set(preds['GO_ID'])
    else:
        pred_terms = set()

    # Build ancestor map
    godag = GODag(args.obo, load_obsolete=True)
    ancestor_map = build_ancestor_map(godag, gaf_terms.union(pred_terms))

    # Load IC
    ic_df = pd.read_csv(args.ic, dtype={'go_term':'category','IC':float}) if args.ic else None

    # Propagate and write, then summary
    stats = {}
    if args.mode in ('gaf','both'):
        gogaf = propagate_gaf(filtered_gaf, ancestor_map, ic_df)
        gogaf.to_csv(args.output_gaf, index=False)
        if ic_df is not None:
            u = gogaf['go_term'].nunique(); m = gogaf['IC'].isna().sum(); stats['GAF']=(u,m,m/u*100)
    if args.mode in ('pred','both'):
        gopred = propagate_deepfri(preds, ancestor_map, ic_df)
        gopred.to_csv(args.output_pred, index=False)
        if ic_df is not None:
            u = gopred['go_term'].nunique(); m = gopred['IC'].isna().sum(); stats['deepFRI']=(u,m,m/u*100)

    if ic_df is not None:
        print("\nIC availability summary:")
        for src,(u,m,p) in stats.items():
            print(f"  {src}: {u} unique GO terms, {m} missing IC ({p:.2f}% missing)")

    print(f"Total execution time: {time.time()-overall_start:.2f}s")
