""" Detect cryptic mutation clusters in wastewater sequencing data """

import warnings
import argparse
import sys, os
import pandas as pd
import pickle
from pathlib import Path
from tqdm import tqdm

from outbreak_tools import crumbs
from outbreak_data import outbreak_data as od

from gisaid_authentication import check_authentication

warnings.simplefilter(action='ignore', category=FutureWarning)

# Date window for clinical data query
START_DATE = "2020-01-01"
END_DATE = "2025-12-31"
FREYJA_BARCODES = "reference/freyja/usher_barcodes.feather"

parser = argparse.ArgumentParser(
    description="Detect cryptic mutation clusters in wastewater sequencing data"
)
parser.add_argument(
    "--covar_dir", help="Directory containing coVar (linked mutations) output", type=str
)
parser.add_argument(
    "--metadata", help="Metadata file", type=str
)
parser.add_argument(
    "--freyja_barcodes", help="Freyja barcodes file", type=str, default=FREYJA_BARCODES
)
parser.add_argument(
    "--output", help="Output file name", default="covar_clinical_detections.tsv"
)

# Hide print statements from API calls
class HiddenPrints:
    def __enter__(self):
        self._original_stdout = sys.stdout
        sys.stdout = open(os.devnull, 'w')

    def __exit__(self, exc_type, exc_val, exc_tb):
        sys.stdout.close()
        sys.stdout = self._original_stdout

def main():
    # Authenticate with GISAID credentials
    try:
        check_authentication()
    except:
        raise Exception("Error authenticating with GISAID credentials\n"
                        "Please run `python gisaid_authentication.py` "
                        "to authenticate and try again")

    args = parser.parse_args()

    aggregate_covariants = parse_covariants(args.covar_dir)

    # Add metadata to aggregate covariants
    aggregate_covariants = join_metadata(aggregate_covariants, args.metadata)

    # Query clinical data
    cryptic_variants = query_clinical_data(
        aggregate_covariants, FREYJA_BARCODES, START_DATE, END_DATE
    )

    # Save output
    cryptic_variants.to_csv(args.output, sep="\t", index=False)



def parse_aa_muts(muts):
    """Prepare query for aa mutations"""
    output = []
    if type(muts) != str:
        return []
    for m in muts.split(' '):
        if m == "Unknown" or m == "NA": # Frameshift indels, most likely due to sequencing errors
            return []
        if m.split(":")[1][0] == m.split(":")[1][-1]: # Omit synonymous mutations
            continue
        if ":INS" in m: # Omit insertions
            continue
        if ":DEL" in m:
            if '/' not in m: # Workaround for single aa deletion query bug (e.g. S:DEL144 -> S:DEL144/144)
                output.append(f'{m}/{m.split(":DEL")[1]}')
            else:
                output.append(m)
        else:
            output.append(m)
    return list(set(output))


def parse_covariants(covar_dir):
    """Parse covar output files, aggregate into one dataframe"""

    agg_covariants = pd.DataFrame()
    for file in os.listdir(covar_dir):
        df = pd.read_csv(f'{covar_dir}/{file}', sep="\t")

        df["query"] = df["aa_mutations"].apply(parse_aa_muts)
        df = df[df["query"].apply(len) > 0]
        df["sample_id"] = file
        agg_covariants = pd.concat([agg_covariants, df])

    return agg_covariants

def join_metadata(aggregate_covariants, metadata_file):
    """Add metadata to aggregate covariants dataframe"""

    metadata = pd.read_csv(metadata_file)
    aggregate_covariants['sample_id'] = aggregate_covariants['sample_id'].str.split('.trimmed').str[0]

    # Join metadata
    aggregate_covariants = aggregate_covariants.merge(
        metadata[metadata.columns],
        on="sample_id",
        how="left",
    )

    return aggregate_covariants


def query_clinical_data(aggregate_covariants, freyja_barcodes, START_DATE, END_DATE, cache_dir=".cache"):
    """Query outbreak.info API for clinical detections of mutation clusters"""
    
    # Create cache directory
    Path(cache_dir).mkdir(exist_ok=True)
    cache_file = Path(cache_dir) / "clinical_data_cache.pkl"
    
    # Load existing cache
    if cache_file.exists():
        with open(cache_file, 'rb') as f:
            cache = pickle.load(f)
    else:
        cache = {}
    
    lineage_key = crumbs.get_alias_key()
    barcode_muts = pd.read_feather(freyja_barcodes).columns
    
    # Pre-filter and deduplicate
    unique_clusters = aggregate_covariants.drop_duplicates(subset=['query'])
    
    # Filter out barcode matches upfront
    unique_clusters = unique_clusters[
        ~unique_clusters['nt_mutations'].apply(
            lambda x: all([m in barcode_muts for m in str(x).split(" ")])
        )
    ]
    
    # Query only uncached clusters
    for idx, row in tqdm(unique_clusters.iterrows(), 
                         total=len(unique_clusters),
                         desc="Querying clinical data"):
        cluster = row["query"]
        cluster_key = str(cluster)
        
        if cluster_key in cache:
            continue
            
        try:
            with HiddenPrints():
                mut_data = od.lineage_cl_prevalence(
                    ".",
                    descendants=True,
                    mutations=cluster,
                    datemin=START_DATE,
                    datemax=END_DATE,
                    lineage_key=lineage_key,
                )
        except NameError as e:
            print(f"Error querying outbreak.info for cluster {cluster}: {e}")
            cache[cluster_key] = 0
            continue

        if mut_data is not None:
            cache[cluster_key] = int(mut_data["lineage_count"].sum())
        else:
            cache[cluster_key] = 0
    
    # Save cache
    with open(cache_file, 'wb') as f:
        pickle.dump(cache, f)
    
    # Map results back
    aggregate_covariants["num_clinical_detections"] = aggregate_covariants["query"].apply(
        lambda x: cache.get(str(x))
    )
    
    aggregate_covariants = aggregate_covariants.dropna(subset=["num_clinical_detections"])
    
    return aggregate_covariants

if __name__ == '__main__':
    main()