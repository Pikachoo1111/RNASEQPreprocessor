import pandas as pd
import json
import gzip

# --------- USER CONFIGURATION ---------
MIRNA_FILE = "data/raw_data/mirna/tcga_acc_mirna_rpm_log2.cct"
GEO_META_FILE = "geo_response_metadata.csv"         # Optional: responder labels
WSI_JSON_PATH = r"study_patient_mapping.json"                 # JSON file as described

# Immune-related miRNAs based on literature
IMMUNE_MIRNAS = [
    # T cell regulation and immune checkpoints
    "hsa-mir-155",      # T cell activation, immune response
    "hsa-mir-146a",     # Immune tolerance, T-reg function
    "hsa-mir-21",       # Immune suppression, M2 macrophages
    "hsa-mir-150",      # T cell development, NK cell function
    "hsa-mir-125b-1",   # T cell differentiation
    "hsa-mir-125b-2",   # T cell differentiation
    
    # Immune checkpoint regulation
    "hsa-mir-200c",     # PD-L1 regulation
    "hsa-mir-200a",     # PD-L1 regulation
    "hsa-mir-200b",     # PD-L1 regulation
    "hsa-mir-34a",      # PD-L1 regulation
    "hsa-mir-138",      # CTLA-4 regulation
    
    # Cytokine and immune signaling
    "hsa-mir-23a",      # IL-17 pathway
    "hsa-mir-23b",      # IL-17 pathway
    "hsa-mir-27a",      # Immune cell migration
    "hsa-mir-27b",      # Immune cell migration
    "hsa-mir-29a",      # Immune cell function
    "hsa-mir-29b-1",    # Immune cell function
    "hsa-mir-29b-2",    # Immune cell function
    "hsa-mir-29c",      # Immune cell function
    
    # Tumor microenvironment
    "hsa-mir-10a",      # Immune cell infiltration
    "hsa-mir-10b",      # Immune cell infiltration
    "hsa-mir-126",      # Angiogenesis, immune cell recruitment
    "hsa-mir-181a-1",   # T cell sensitivity
    "hsa-mir-181a-2",   # T cell sensitivity
    "hsa-mir-181b-1",   # T cell sensitivity
    "hsa-mir-181b-2",   # T cell sensitivity
    
    # Interferon response
    "hsa-mir-100",      # Interferon signaling
    "hsa-mir-101-1",    # Interferon response
    "hsa-mir-101-2",    # Interferon response
    
    # Let-7 family (immune regulation)
    "hsa-let-7a-1",     # Immune cell development
    "hsa-let-7a-2",     # Immune cell development
    "hsa-let-7a-3",     # Immune cell development
    "hsa-let-7b",       # Immune cell development
    "hsa-let-7c",       # Immune cell development
    "hsa-let-7d",       # Immune cell development
    "hsa-let-7e",       # Immune cell development
    "hsa-let-7f-1",     # Immune cell development
    "hsa-let-7f-2",     # Immune cell development
    
    # Additional immune-related miRNAs
    "hsa-mir-17",       # T cell proliferation
    "hsa-mir-18a",      # T cell proliferation
    "hsa-mir-19a",      # T cell proliferation
    "hsa-mir-19b-1",    # T cell proliferation
    "hsa-mir-19b-2",    # T cell proliferation
    "hsa-mir-20a",      # T cell proliferation
    "hsa-mir-92a-1",    # T cell proliferation
    "hsa-mir-92a-2",    # T cell proliferation
    
    # Macrophage polarization
    "hsa-mir-124-1",    # M1/M2 macrophage polarization
    "hsa-mir-124-2",    # M1/M2 macrophage polarization
    "hsa-mir-124-3",    # M1/M2 macrophage polarization
    "hsa-mir-223",      # Macrophage function
    
    # NK cell function
    "hsa-mir-30a",      # NK cell cytotoxicity
    "hsa-mir-30b",      # NK cell cytotoxicity
    "hsa-mir-30c-1",    # NK cell cytotoxicity
    "hsa-mir-30c-2",    # NK cell cytotoxicity
    "hsa-mir-30d",      # NK cell cytotoxicity
    "hsa-mir-30e",      # NK cell cytotoxicity
]
# --------------------------------------

# --------- 1. LOAD MIRNA DATA ---------
print("Loading miRNA data...")
try:
    with gzip.open(MIRNA_FILE, 'rt') as f:
        mirna_df = pd.read_csv(f, sep="\t", index_col=0)
except Exception as e:
    print(f"Error reading compressed file: {e}")
    print("Trying to read as regular file...")
    mirna_df = pd.read_csv(MIRNA_FILE, sep="\t", index_col=0)

print(f"Loaded miRNA data: {mirna_df.shape[0]} miRNAs, {mirna_df.shape[1]} samples")

# Extract case IDs from sample barcodes
def get_case_id(sample_id):
    return "-".join(sample_id.split(".")[:3]).replace(".", "-")

mirna_df.columns = [get_case_id(col) for col in mirna_df.columns]
mirna_df = mirna_df.loc[~mirna_df.index.duplicated(keep='first')]  # Remove duplicate miRNAs

# --------- 2. FILTER FOR IMMUNE-RELATED MIRNAS ---------
print("Filtering immune-related miRNAs...")
print(f"Looking for {len(IMMUNE_MIRNAS)} immune-related miRNAs...")

# Find available immune miRNAs in the dataset
available_immune_mirnas = [mirna for mirna in IMMUNE_MIRNAS if mirna in mirna_df.index]
print(f"Found {len(available_immune_mirnas)} immune miRNAs in dataset:")
for mirna in available_immune_mirnas:
    print(f"  - {mirna}")

if not available_immune_mirnas:
    print("WARNING: No immune miRNAs found with exact names. Checking for partial matches...")
    # Try partial matching
    available_immune_mirnas = []
    for mirna in IMMUNE_MIRNAS:
        partial_matches = [idx for idx in mirna_df.index if mirna.replace("hsa-", "") in idx]
        if partial_matches:
            available_immune_mirnas.extend(partial_matches)
    
    available_immune_mirnas = list(set(available_immune_mirnas))  # Remove duplicates
    print(f"Found {len(available_immune_mirnas)} immune miRNAs with partial matching:")
    for mirna in available_immune_mirnas[:10]:  # Show first 10
        print(f"  - {mirna}")
    if len(available_immune_mirnas) > 10:
        print(f"  ... and {len(available_immune_mirnas) - 10} more")

if available_immune_mirnas:
    immune_df = mirna_df.loc[available_immune_mirnas]
else:
    print("WARNING: No immune miRNAs found. Using top 50 most variable miRNAs instead...")
    # Calculate variance for each miRNA and select top 50
    mirna_var = mirna_df.var(axis=1).sort_values(ascending=False)
    top_mirnas = mirna_var.head(50).index.tolist()
    immune_df = mirna_df.loc[top_mirnas]
    print(f"Selected top {len(top_mirnas)} most variable miRNAs")

# Transpose so rows = samples
immune_df = immune_df.T
immune_df.index.name = "case_id"

print(f"Final immune miRNA matrix: {immune_df.shape[0]} samples, {immune_df.shape[1]} miRNAs")

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
    if "num_patches" in immune_df.columns:
        immune_df["num_patches"] = immune_df["num_patches"].fillna("No WSI Data")
    if "file" in immune_df.columns:
        immune_df["file"] = immune_df["file"].fillna("No WSI Data")

except FileNotFoundError:
    print("WARNING: WSI JSON file not found - skipping this step.")

# --------- 5. EXPORT ---------
output_file = "data/processed_data/mirna_immune_markers_with_metadata.csv"
immune_df.to_csv(output_file)
print(f"SUCCESS: Done. Output saved to: {output_file}")
print(f"Final dataset shape: {immune_df.shape}")

# --------- 6. SUMMARY STATISTICS ---------
print("\nSummary Statistics:")
print(f"Number of samples: {immune_df.shape[0]}")
print(f"Number of immune miRNAs: {immune_df.shape[1] - (2 if 'num_patches' in immune_df.columns else 0)}")

if available_immune_mirnas:
    print(f"\nTop 5 most variable immune miRNAs:")
    mirna_cols = [col for col in immune_df.columns if col not in ['num_patches', 'file', 'response']]
    if mirna_cols:
        mirna_variance = immune_df[mirna_cols].var().sort_values(ascending=False)
        for i, (mirna, var) in enumerate(mirna_variance.head(5).items()):
            print(f"  {i+1}. {mirna}: variance = {var:.4f}")
