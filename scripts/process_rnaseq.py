import pandas as pd
import json

# --------- USER CONFIGURATION ---------
RNA_FILE = "data/raw_data/rnaseq/tcga_acc_rnaseq_rsem_log2.cct"
GEO_META_FILE = "geo_response_metadata.csv"         # Optional: responder labels
WSI_JSON_PATH = r"study_patient_mapping.json"                 # ✅ JSON file as described
IMMUNE_MARKERS = ["CD274", "CTLA4", "PDCD1", "LAG3", "TIGIT", "CD8A", "FOXP3", "GZMB"]
# --------------------------------------

# --------- 1. LOAD RNA-SEQ DATA ---------
print("Loading RNA-seq data...")
rna_df = pd.read_csv(RNA_FILE, sep="\t", index_col=0)

# Extract case IDs from sample barcodes
def get_case_id(sample_id):
    return "-".join(sample_id.split("-")[:3])

rna_df.columns = [get_case_id(col) for col in rna_df.columns]
rna_df = rna_df.loc[~rna_df.index.duplicated(keep='first')]  # Remove duplicate genes

# --------- 2. FILTER FOR IMMUNE MARKERS ---------
print("Filtering immune-related genes...")
immune_df = rna_df[rna_df.index.isin(IMMUNE_MARKERS)]

# Transpose so rows = samples
immune_df = immune_df.T
immune_df.index.name = "case_id"

# --------- 3. OPTIONAL: ADD RESPONSE LABELS ---------
try:
    print("Adding responder/non-responder labels...")
    geo_meta = pd.read_csv(GEO_META_FILE)
    geo_meta["case_id"] = geo_meta["case_id"].astype(str)
    immune_df = immune_df.merge(geo_meta[["case_id", "response"]], on="case_id", how="left")
except FileNotFoundError:
    print("WARNING: GEO response label file not found - skipping this step.")

# --------- 4. LOAD & MERGE WSI METADATA ---------
try:
    print("Merging with WSI metadata from JSON...")
    with open(WSI_JSON_PATH, "r") as f:
        wsi_json = json.load(f)

    acc_cases = wsi_json.get("acc", {})
    wsi_df = pd.DataFrame.from_dict(acc_cases, orient="index")
    wsi_df.index.name = "case_id"

    # Normalize case_ids so dots and dashes match
    immune_df.index = immune_df.index.str.replace(r"\.", "-", regex=True)
    wsi_df.index = wsi_df.index.str.replace(r"\.", "-", regex=True)

    immune_df = immune_df.merge(wsi_df, on="case_id", how="left")

    # Fill missing values in WSI columns with placeholder
    immune_df["num_patches"] = immune_df["num_patches"].fillna("No WSI Data")
    immune_df["file"] = immune_df["file"].fillna("No WSI Data")

except FileNotFoundError:
    print("WARNING: WSI JSON file not found - skipping this step.")


# --------- 5. EXPORT ---------
output_file = "data/processed_data/rnaseq_immune_markers_with_metadata.csv"
immune_df.to_csv(output_file)
print(f"SUCCESS: Done. Output saved to: {output_file}")
