import pandas as pd
import pyarrow as pa
import pyarrow.parquet as pq
import hashlib
import os

########################
# Utility Functions
########################

def hash_id(sample_id: str) -> str:
    """Hashes sample identifiers for anonymization."""
    return hashlib.sha256(sample_id.encode()).hexdigest()[:16]

def load_1000_genomes(vcf_path: str) -> pd.DataFrame:
    """Parse VCF genotype files (basic parser, no cyvcf2/vcfpy)."""
    rows = []
    with open(vcf_path, "r") as f:
        for line in f:
            if line.startswith("##"):
                continue  # skip metadata
            if line.startswith("#CHROM"):
                header = line.strip().split("\t")
                samples = header[9:]  # sample IDs
                continue
            parts = line.strip().split("\t")
            chrom, pos, var_id, ref, alt = parts[0], parts[1], parts[2], parts[3], parts[4]
            genotypes = parts[9:]
            for sample, gt in zip(samples, genotypes):
                gt_clean = gt.split(":")[0]  # only GT field (e.g., "0/1")
                rows.append({
                    "Sample_ID": hash_id(sample),
                    "Data_Type": "Genotype",
                    "Gene/Protein": var_id if var_id != "." else chrom,
                    "Expression/Variant": gt_clean,
                    "Metadata": {"pos": pos, "ref": ref, "alt": alt}
                })
    return pd.DataFrame(rows)


def load_gtex(expression_path: str) -> pd.DataFrame:
    """Load GTEx gene expression matrices."""
    df = pd.read_csv(expression_path, sep="\t", index_col=0)
    melted = df.melt(var_name="Sample_ID", value_name="Expression", ignore_index=False)
    melted.reset_index(inplace=True)
    melted.rename(columns={"index": "Gene"}, inplace=True)
    melted["Sample_ID"] = melted["Sample_ID"].apply(hash_id)
    melted["Data_Type"] = "GeneExpression"
    melted["Metadata"] = None
    melted.rename(columns={"Gene": "Gene/Protein", "Expression": "Expression/Variant"}, inplace=True)
    return melted


def load_pride(proteomics_path: str) -> pd.DataFrame:
    """Load PRIDE proteomics datasets."""
    df = pd.read_csv(proteomics_path)
    df["Sample_ID"] = df["Sample_ID"].apply(hash_id)
    df["Data_Type"] = "ProteinExpression"
    df.rename(columns={"Protein": "Gene/Protein", "Intensity": "Expression/Variant"}, inplace=True)

    # Add Metadata column if missing
    if "Metadata" not in df.columns:
        df["Metadata"] = None

    return df[["Sample_ID", "Data_Type", "Gene/Protein", "Expression/Variant", "Metadata"]]


########################
# Harmonization & Export
########################
def harmonize_and_export(datasets: list, output_file: str = "genome_phoenix.parquet"):
    """Merge all datasets into one standardized Parquet table."""
    merged = pd.concat(datasets, ignore_index=True)

    # Ensure consistent datatypes
    merged["Expression/Variant"] = merged["Expression/Variant"].astype(str)
    merged["Gene/Protein"] = merged["Gene/Protein"].astype(str)

    table = pa.Table.from_pandas(merged, preserve_index=False)
    pq.write_table(table, output_file)
    print(f"[+] Harmonized dataset written to {output_file}")



########################
# Query Function
########################

def query_data(parquet_file: str, filters: dict) -> pd.DataFrame:
    """Filter harmonized dataset by metadata fields."""
    table = pq.read_table(parquet_file).to_pandas()
    query_df = table.copy()
    for key, value in filters.items():
        query_df = query_df[query_df["Metadata"].apply(lambda m: m.get(key) == value if isinstance(m, dict) else False)]
    return query_df


########################
# Example Workflow
########################

if __name__ == "__main__":
    # Example paths (replace with actual datasets)
    g1 = load_1000_genomes("data/example.vcf")
    g2 = load_gtex("data/expression.tsv")
    g3 = load_pride("data/proteomics.csv")

    harmonize_and_export([g1, g2, g3])

    # Example query
    result = query_data("genome_phoenix.parquet", {"tissue": "liver", "population": "AFR"})
    print(result.head())
