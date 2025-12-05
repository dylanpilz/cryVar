""" Detect cryptic SARS-CoV-2 mutation clusters in wastewater sequencing data """

import warnings
import argparse
import sys, os
import pandas as pd
import pickle
from datetime import datetime
from pathlib import Path
from tqdm import tqdm
from multiprocessing import Pool, cpu_count

from outbreak_tools import crumbs
from outbreak_data import outbreak_data as od

from gisaid_authentication import check_authentication

warnings.simplefilter(action='ignore', category=FutureWarning)

# Date window for clinical data query
START_DATE = "2020-01-01"
END_DATE = datetime.now().strftime("%Y-%m-%d")

parser = argparse.ArgumentParser(
    description="Detect cryptic mutation clusters in wastewater sequencing data"
)
parser.add_argument(
    "--covar_dir", help="Directory containing coVar (linked mutations) output", type=str, required=True, default="results/covar"
)
parser.add_argument(
    "--metadata", help="Metadata file", type=str, required=True
)
parser.add_argument(
    "--output_dir", help="Output directory", default="results/detect_cryptic"
)
parser.add_argument(
    "--processes", help="Number of parallel processes to use", type=int, default=None
)
parser.add_argument(
    "--max_clinical_detections", help="Maximum number of clinical detections to consider a cryptic variant", type=int, default=10
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
    covariants_with_clinical_data = query_clinical_data(
        aggregate_covariants, START_DATE, END_DATE, 
        processes=args.processes
    )

    # Save output
    covariants_with_clinical_data.to_csv(f'{args.output_dir}/covar_clinical_detections.tsv', sep="\t", index=False)


    # Filter for cryptic variants
    cryptic_variants = covariants_with_clinical_data[covariants_with_clinical_data["num_clinical_detections"] <= args.max_clinical_detections]
    cryptic_variants = cryptic_variants[cryptic_variants["query"].apply(len) >= 2]

    # select for clusters that appear at least 2 times in the wastewater data
    counts = cryptic_variants.groupby("nt_mutations").size().reset_index(name="count")
    counts = counts[counts["count"] >= 2]
    cryptic_variants = cryptic_variants[cryptic_variants["nt_mutations"].isin(counts["nt_mutations"])]
    cryptic_variants.to_csv(f'{args.output_dir}/cryptic_variants.tsv', sep="\t", index=False)


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
        df["sample_id"] = file.split('.covar.tsv')[0]
        agg_covariants = pd.concat([agg_covariants, df])

    return agg_covariants

def join_metadata(aggregate_covariants, metadata_file):
    """Add metadata to aggregate covariants dataframe"""

    metadata = pd.read_csv(metadata_file)
    # Join metadata
    aggregate_covariants = aggregate_covariants.merge(
        metadata[metadata.columns],
        on="sample_id",
        how="left",
    )

    return aggregate_covariants


def query_single_cluster(args_tuple):
    """Worker function to query a single cluster"""
    cluster, START_DATE, END_DATE, lineage_key = args_tuple
    cluster_key = str(cluster)
    
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
    except Exception as e:
        print(f"Error querying outbreak.info for cluster {cluster}: {e}")
        return (cluster_key, 0)

    if mut_data is not None:
        return (cluster_key, int(mut_data["lineage_count"].sum()))
    else:
        return (cluster_key, 0)


def query_clinical_data(aggregate_covariants, START_DATE, END_DATE, 
                        cache_dir=".cache", processes=None):
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
    
    # Pre-filter and deduplicate
    unique_clusters = aggregate_covariants.drop_duplicates(subset=['query'])
    
    # Filter out already cached clusters
    uncached_clusters = []
    for idx, row in unique_clusters.iterrows():
        cluster = row["query"]
        cluster_key = str(cluster)
        if cluster_key not in cache:
            uncached_clusters.append(cluster)
    
    print(uncached_clusters)
    # Query uncached clusters in parallel
    if uncached_clusters:
        if processes is None:
            processes = cpu_count()
        
        # Prepare arguments for worker function
        query_args = [(cluster, START_DATE, END_DATE, lineage_key) 
                     for cluster in uncached_clusters]
        
        # Query in parallel
        with Pool(processes=processes) as pool:
            results = list(tqdm(
                pool.imap(query_single_cluster, query_args),
                total=len(query_args),
                desc="Querying clinical data"
            ))
        
        # Update cache with results
        for cluster_key, count in results:
            cache[cluster_key] = count
    
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