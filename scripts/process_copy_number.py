import pandas as pd
import json
import gzip

# --------- USER CONFIGURATION ---------
CN_GENE_FILE = "data/raw_data/copy_number/tcga_acc_copy_number_gene_level.cct.gz"
CN_FOCAL_FILE = "data/raw_data/copy_number/tcga_acc_copy_number_focal_level.cct"
GEO_META_FILE = "geo_response_metadata.csv"         # Optional: responder labels
WSI_JSON_PATH = r"study_patient_mapping.json"                 # JSON file as described

# Immune-related genes for copy number analysis
IMMUNE_GENES = [
    # Immune checkpoint genes
    "CD274",    # PD-L1 (9p24.1)
    "PDCD1",    # PD-1 (2q37.3)
    "CTLA4",    # CTLA-4 (2q33.2)
    "LAG3",     # LAG-3 (12p13.32)
    "TIGIT",    # TIGIT (3q13.31)
    "HAVCR2",   # TIM-3 (5q33.3)
    "BTLA",     # BTLA (3q13.2)
    "ICOS",     # ICOS (2q33.3)
    "ICOSLG",   # ICOS-L (21q22.3)
    
    # T cell function genes
    "CD8A",     # CD8A (2p11.2)
    "CD8B",     # CD8B (2p11.2)
    "CD4",      # CD4 (12p13.31)
    "FOXP3",    # FOXP3 (Xp11.23)
    "GZMB",     # Granzyme B (14q12)
    "PRF1",     # Perforin 1 (10q22.1)
    "IL2",      # IL-2 (4q27)
    "IFNG",     # IFN-gamma (12q15)
    
    # MHC genes (6p21.3)
    "HLA-A",    "HLA-B",    "HLA-C",
    "HLA-DRA",  "HLA-DRB1", "HLA-DQA1", "HLA-DQB1",
    "B2M",      # Beta-2 microglobulin (15q21.1)
    "TAP1",     # TAP1 (6p21.32)
    "TAP2",     # TAP2 (6p21.32)
    
    # Cytokines and chemokines
    "TNF",      # TNF-alpha (6p21.33)
    "IL1B",     # IL-1 beta (2q14.1)
    "IL6",      # IL-6 (7p15.3)
    "IL10",     # IL-10 (1q32.1)
    "TGFB1",    # TGF-beta1 (19q13.2)
    "CCL2",     # MCP-1 (17q12)
    "CCL5",     # RANTES (17q12)
    "CXCL9",    # MIG (4q21.1)
    "CXCL10",   # IP-10 (4q21.1)
    
    # Immune signaling pathways
    "STAT1",    # STAT1 (2q32.2)
    "STAT3",    # STAT3 (17q21.2)
    "STAT4",    # STAT4 (2q32.2-q32.3)
    "IRF1",     # IRF1 (5q31.1)
    "IRF3",     # IRF3 (19q13.33)
    "IRF7",     # IRF7 (11p15.5)
    "NFKB1",    # NF-kB1 (4q24)
    "RELA",     # RelA (11q13.1)
    
    # Complement genes
    "C1QA",     # C1QA (1p36.12)
    "C1QB",     # C1QB (1p36.12)
    "C1QC",     # C1QC (1p36.12)
    "C3",       # C3 (19p13.3)
    "C5",       # C5 (9q33.2)
    
    # NK cell genes
    "NCAM1",    # CD56 (11q23.2)
    "KLRK1",    # NKG2D (12p13.2)
    "KIR2DL1",  # KIR genes (19q13.4)
    "KIR2DL3",  "KIR3DL1",
    
    # B cell genes
    "CD19",     # CD19 (16p11.2)
    "MS4A1",    # CD20 (11q12.2)
    "CD79A",    # CD79A (19q13.2)
    "CD79B",    # CD79B (17q23.3)
    
    # Macrophage genes
    "CD68",     # CD68 (17p13.1)
    "CD163",    # CD163 (12p13.31)
    "CD86",     # CD86 (3q13.33)
    "CD80",     # CD80 (3q13.33)
]

# Important immune-related chromosomal regions
IMMUNE_REGIONS = [
    "1p36",     # Complement genes
    "2q32",     # STAT genes
    "2q33",     # CTLA4, ICOS
    "3q13",     # TIGIT, BTLA, CD86, CD80
    "4q21",     # Chemokine cluster
    "4q24",     # NFKB1
    "5q31",     # IRF1, cytokine cluster
    "6p21",     # MHC region
    "9p24",     # PD-L1 region
    "11q13",    # RELA
    "12p13",    # LAG3, CD4
    "14q12",    # GZMB
    "17q12",    # Chemokine cluster
    "19q13",    # Multiple immune genes
]
# --------------------------------------

def load_copy_number_data(file_path):
    """Load copy number data and handle different formats"""
    try:
        if file_path.endswith('.gz'):
            with gzip.open(file_path, 'rt') as f:
                df = pd.read_csv(f, sep="\t", index_col=0)
        else:
            df = pd.read_csv(file_path, sep="\t", index_col=0)
        print(f"Loaded copy number data from {file_path}: {df.shape[0]} features, {df.shape[1]} samples")
        return df
    except Exception as e:
        print(f"Error loading {file_path}: {e}")
        return None

def get_case_id(sample_id):
    """Extract case ID from sample barcode"""
    return "-".join(sample_id.split(".")[:3]).replace(".", "-")

def find_immune_features(df, feature_list):
    """Find immune-related features in the dataset"""
    available_features = []
    
    # Direct matching
    for feature in feature_list:
        if feature in df.index:
            available_features.append(feature)
    
    # Partial matching for gene names
    for feature in feature_list:
        partial_matches = [idx for idx in df.index if feature.upper() in idx.upper()]
        available_features.extend(partial_matches)
    
    # Remove duplicates while preserving order
    available_features = list(dict.fromkeys(available_features))
    
    return available_features

# --------- 1. LOAD COPY NUMBER DATA ---------
print("Loading copy number data...")

# Load gene-level copy number data
cn_gene_df = load_copy_number_data(CN_GENE_FILE)

# Load focal-level copy number data
cn_focal_df = load_copy_number_data(CN_FOCAL_FILE)

# Use gene-level data as primary, focal as supplementary
if cn_gene_df is not None:
    cn_df = cn_gene_df
    print("Using gene-level copy number data")
elif cn_focal_df is not None:
    cn_df = cn_focal_df
    print("Using focal-level copy number data")
else:
    print("❌ No copy number data could be loaded!")
    exit(1)

# Extract case IDs from sample barcodes
cn_df.columns = [get_case_id(col) for col in cn_df.columns]
cn_df = cn_df.loc[~cn_df.index.duplicated(keep='first')]  # Remove duplicates

# --------- 2. FILTER FOR IMMUNE-RELATED GENES/REGIONS ---------
print("Filtering immune-related copy number alterations...")
print(f"Looking for immune-related features from {len(IMMUNE_GENES)} candidates...")

available_immune_features = find_immune_features(cn_df, IMMUNE_GENES)

print(f"Found {len(available_immune_features)} immune-related features in dataset:")
for feature in available_immune_features[:20]:  # Show first 20
    print(f"  - {feature}")
if len(available_immune_features) > 20:
    print(f"  ... and {len(available_immune_features) - 20} more")

# If using focal data, also look for immune regions
if cn_focal_df is not None and available_immune_features:
    print("\nLooking for immune-related chromosomal regions in focal data...")
    focal_immune_regions = []
    for region in IMMUNE_REGIONS:
        matching_regions = [idx for idx in cn_focal_df.index if region in idx]
        focal_immune_regions.extend(matching_regions)
    
    if focal_immune_regions:
        print(f"Found {len(focal_immune_regions)} immune-related focal regions:")
        for region in focal_immune_regions[:10]:  # Show first 10
            print(f"  - {region}")

if available_immune_features:
    immune_cn_df = cn_df.loc[available_immune_features]
else:
    print("WARNING: No immune features found. Using most variable features...")
    # Select top 100 most variable features
    feature_variance = cn_df.var(axis=1).sort_values(ascending=False)
    top_features = feature_variance.head(100).index.tolist()
    immune_cn_df = cn_df.loc[top_features]
    print(f"Selected top {len(top_features)} most variable features")

# Transpose so rows = samples
immune_cn_df = immune_cn_df.T
immune_cn_df.index.name = "case_id"

print(f"Final immune copy number matrix: {immune_cn_df.shape[0]} samples, {immune_cn_df.shape[1]} features")

# --------- 3. CALCULATE COPY NUMBER SUMMARY STATISTICS ---------
print("Calculating copy number summary statistics...")

# Calculate overall genomic instability metrics
immune_cn_df['total_alterations'] = (immune_cn_df.abs() > 0.3).sum(axis=1)  # Count significant alterations
immune_cn_df['amplification_burden'] = (immune_cn_df > 0.3).sum(axis=1)     # Count amplifications
immune_cn_df['deletion_burden'] = (immune_cn_df < -0.3).sum(axis=1)        # Count deletions
immune_cn_df['cn_instability_score'] = immune_cn_df.iloc[:, :-3].std(axis=1)  # Copy number variability

# Identify samples with high genomic instability
q75_instability = immune_cn_df['cn_instability_score'].quantile(0.75)
immune_cn_df['high_instability'] = (immune_cn_df['cn_instability_score'] > q75_instability).astype(int)

print(f"Copy number alteration summary:")
print(f"  Mean alterations per sample: {immune_cn_df['total_alterations'].mean():.1f}")
print(f"  Mean amplifications per sample: {immune_cn_df['amplification_burden'].mean():.1f}")
print(f"  Mean deletions per sample: {immune_cn_df['deletion_burden'].mean():.1f}")
print(f"  High instability samples: {immune_cn_df['high_instability'].sum()}")

# --------- 4. OPTIONAL: ADD RESPONSE LABELS ---------
try:
    print("Adding responder/non-responder labels...")
    geo_meta = pd.read_csv(GEO_META_FILE)
    geo_meta["case_id"] = geo_meta["case_id"].astype(str)
    immune_cn_df = immune_cn_df.merge(geo_meta[["case_id", "response"]], on="case_id", how="left")
except FileNotFoundError:
    print("WARNING: GEO response label file not found - skipping this step.")

# --------- 5. LOAD & MERGE WSI METADATA ---------
try:
    print("Merging with WSI metadata from JSON...")
    with open(WSI_JSON_PATH, "r") as f:
        wsi_json = json.load(f)

    acc_cases = wsi_json.get("acc", {})
    wsi_df = pd.DataFrame.from_dict(acc_cases, orient="index")
    wsi_df.index.name = "case_id"

    # Normalize case_ids so dots and dashes match
    immune_cn_df.index = immune_cn_df.index.str.replace(r"\.", "-", regex=True)
    wsi_df.index = wsi_df.index.str.replace(r"\.", "-", regex=True)

    immune_cn_df = immune_cn_df.merge(wsi_df, on="case_id", how="left")

    # Fill missing values in WSI columns with placeholder
    if "num_patches" in immune_cn_df.columns:
        immune_cn_df["num_patches"] = immune_cn_df["num_patches"].fillna("No WSI Data")
    if "file" in immune_cn_df.columns:
        immune_cn_df["file"] = immune_cn_df["file"].fillna("No WSI Data")

except FileNotFoundError:
    print("WARNING: WSI JSON file not found - skipping this step.")

# --------- 6. EXPORT ---------
output_file = "data/processed_data/copy_number_immune_features_with_metadata.csv"
immune_cn_df.to_csv(output_file)
print(f"SUCCESS: Done. Output saved to: {output_file}")
print(f"Final dataset shape: {immune_cn_df.shape}")

# --------- 7. SUMMARY STATISTICS ---------
print("\nSummary Statistics:")
print(f"Number of samples: {immune_cn_df.shape[0]}")

# Count different types of columns
cn_cols = [col for col in immune_cn_df.columns if col not in ['total_alterations', 'amplification_burden', 'deletion_burden', 'cn_instability_score', 'high_instability', 'num_patches', 'file', 'response']]
print(f"Number of copy number features: {len(cn_cols)}")

print(f"\nGenomic Instability Metrics:")
print(f"  Total alterations: {immune_cn_df['total_alterations'].mean():.1f} ± {immune_cn_df['total_alterations'].std():.1f}")
print(f"  Amplification burden: {immune_cn_df['amplification_burden'].mean():.1f} ± {immune_cn_df['amplification_burden'].std():.1f}")
print(f"  Deletion burden: {immune_cn_df['deletion_burden'].mean():.1f} ± {immune_cn_df['deletion_burden'].std():.1f}")
print(f"  CN instability score: {immune_cn_df['cn_instability_score'].mean():.3f} ± {immune_cn_df['cn_instability_score'].std():.3f}")
print(f"  High instability samples: {immune_cn_df['high_instability'].sum()} ({immune_cn_df['high_instability'].mean()*100:.1f}%)")

if cn_cols:
    print(f"\nTop 5 most altered immune features:")
    cn_alteration_freq = (immune_cn_df[cn_cols].abs() > 0.3).mean().sort_values(ascending=False)
    for i, (feature, freq) in enumerate(cn_alteration_freq.head(5).items()):
        print(f"  {i+1}. {feature}: altered in {freq*100:.1f}% of samples")
