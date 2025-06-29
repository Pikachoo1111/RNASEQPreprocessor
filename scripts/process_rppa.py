import pandas as pd
import json

# --------- USER CONFIGURATION ---------
RPPA_GENE_FILE = "data/raw_data/rppa/tcga_acc_rppa_gene_level.cct"
RPPA_ANALYTE_FILE = "data/raw_data/rppa/tcga_acc_rppa_analyte_level.cct"
GEO_META_FILE = "geo_response_metadata.csv"         # Optional: responder labels
WSI_JSON_PATH = r"study_patient_mapping.json"                 # JSON file as described

# Immune-related proteins and pathways for RPPA analysis
IMMUNE_PROTEINS = [
    # Immune checkpoint proteins
    "PD-L1", "PDL1", "CD274",
    "PD-1", "PD1", "PDCD1",
    "CTLA4", "CTLA-4",
    "LAG3", "LAG-3",
    "TIM3", "HAVCR2",
    "TIGIT",
    
    # T cell signaling
    "CD3", "CD3E", "CD3G",
    "CD4", "CD8", "CD8A",
    "TCR", "TCRB",
    "CD28", "CD80", "CD86",
    "ICOS", "ICOSLG",
    
    # Cytokines and growth factors
    "TNF", "TNFA", "TNF-α",
    "IL2", "IL-2",
    "IL6", "IL-6",
    "IL10", "IL-10",
    "IFNG", "IFN-γ", "IFNGR1",
    "TGFB", "TGFB1", "TGF-β",
    "IL1B", "IL-1β",
    "IL17", "IL17A",
    
    # Transcription factors
    "STAT1", "STAT1_pY701",
    "STAT3", "STAT3_pY705", "STAT3_pS727",
    "STAT4", "STAT5", "STAT6",
    "FOXP3", "FOXO1", "FOXO3A",
    "NFKB", "NFKB1", "RELA", "RELB",
    "IRF1", "IRF3", "IRF4", "IRF7",
    "CREB", "CREB_pS133",
    "JUN", "JUNB", "JUND",
    "FOS", "FOSB",
    
    # Kinases and signaling
    "AKT", "AKT_pS473", "AKT_pT308",
    "MTOR", "MTOR_pS2448",
    "P70S6K1", "P70S6K1_pT389",
    "S6", "S6_pS235_S236", "S6_pS240_S244",
    "ERK1", "ERK2", "MAPK1", "MAPK3",
    "P38", "P38_pT180_Y182",
    "JNK", "JNK_pT183_Y185",
    "SRC", "SRC_pY416", "SRC_pY527",
    "LYN", "LYN_pY397",
    "SYK", "SYK_pY525_Y526",
    "ZAP70", "ZAP70_pY319",
    
    # Cell cycle and apoptosis
    "P53", "P53_pS15",
    "MDM2", "MDM2_pS166",
    "P21", "CDKN1A",
    "P27", "CDKN1B",
    "RB", "RB_pS807_S811",
    "CYCLIN_D1", "CCND1",
    "CYCLIN_E1", "CCNE1",
    "BAD", "BAD_pS112",
    "BAX", "BCL2", "BCLXL",
    "CASP3", "CASP7", "CASP9",
    "PARP", "PARP_cleaved",
    
    # Metabolism
    "AMPKA", "AMPKA_pT172",
    "ACC", "ACC_pS79",
    "FASN", "ACLY",
    "HIF1A", "HIF1A_pS797",
    "LDHA", "PKM2",
    
    # Angiogenesis
    "VEGFR2", "KDR", "VEGFR2_pY1175",
    "VEGFA", "VEGFB", "VEGFC",
    "PDGFR", "PDGFRB", "PDGFRB_pY751",
    "EGFR", "EGFR_pY1068", "EGFR_pY1173",
    "HER2", "ERBB2", "HER2_pY1248",
    "HER3", "ERBB3", "HER3_pY1289",
    
    # DNA damage response
    "ATM", "ATM_pS1981",
    "ATR", "ATR_pT1989",
    "CHEK1", "CHEK1_pS345",
    "CHEK2", "CHEK2_pT68",
    "BRCA1", "BRCA1_pS1524",
    "H2AX", "H2AX_pS139",
    "53BP1", "TP53BP1",
    
    # Autophagy
    "LC3A", "LC3B", "MAP1LC3A", "MAP1LC3B",
    "BECN1", "BECLIN1",
    "ATG5", "ATG7", "ATG12",
    "ULK1", "ULK1_pS757",
    
    # Complement and innate immunity
    "C3", "C5", "C1QA", "C1QB", "C1QC",
    "TLR4", "TLR2", "TLR9",
    "MYD88", "IRAK1", "IRAK4",
]
# --------------------------------------

def load_rppa_data(file_path):
    """Load RPPA data and handle different formats"""
    try:
        df = pd.read_csv(file_path, sep="\t", index_col=0)
        print(f"Loaded RPPA data from {file_path}: {df.shape[0]} proteins, {df.shape[1]} samples")
        return df
    except Exception as e:
        print(f"Error loading {file_path}: {e}")
        return None

def get_case_id(sample_id):
    """Extract case ID from sample barcode"""
    return "-".join(sample_id.split(".")[:3]).replace(".", "-")

def find_immune_proteins(df, protein_list):
    """Find immune-related proteins in the dataset"""
    available_proteins = []
    
    # Direct matching
    for protein in protein_list:
        if protein in df.index:
            available_proteins.append(protein)
    
    # Partial matching for protein names
    for protein in protein_list:
        partial_matches = [idx for idx in df.index if protein.upper() in idx.upper()]
        available_proteins.extend(partial_matches)
    
    # Remove duplicates while preserving order
    available_proteins = list(dict.fromkeys(available_proteins))
    
    return available_proteins

# --------- 1. LOAD RPPA DATA ---------
print("Loading RPPA data...")

# Load gene-level RPPA data
rppa_gene_df = load_rppa_data(RPPA_GENE_FILE)

# Load analyte-level RPPA data (includes phosphorylation states)
rppa_analyte_df = load_rppa_data(RPPA_ANALYTE_FILE)

# Choose which dataset to use (analyte-level has more detailed protein states)
if rppa_analyte_df is not None:
    rppa_df = rppa_analyte_df
    print("Using analyte-level RPPA data (includes phosphorylation states)")
elif rppa_gene_df is not None:
    rppa_df = rppa_gene_df
    print("Using gene-level RPPA data")
else:
    print("❌ No RPPA data could be loaded!")
    exit(1)

# Extract case IDs from sample barcodes
rppa_df.columns = [get_case_id(col) for col in rppa_df.columns]
rppa_df = rppa_df.loc[~rppa_df.index.duplicated(keep='first')]  # Remove duplicate proteins

# --------- 2. FILTER FOR IMMUNE-RELATED PROTEINS ---------
print("Filtering immune-related proteins...")
print(f"Looking for immune-related proteins from {len(IMMUNE_PROTEINS)} candidates...")

available_immune_proteins = find_immune_proteins(rppa_df, IMMUNE_PROTEINS)

print(f"Found {len(available_immune_proteins)} immune-related proteins in dataset:")
for protein in available_immune_proteins[:20]:  # Show first 20
    print(f"  - {protein}")
if len(available_immune_proteins) > 20:
    print(f"  ... and {len(available_immune_proteins) - 20} more")

if available_immune_proteins:
    immune_rppa_df = rppa_df.loc[available_immune_proteins]
else:
    print("WARNING: No immune proteins found. Using all available proteins...")
    immune_rppa_df = rppa_df

# Transpose so rows = samples
immune_rppa_df = immune_rppa_df.T
immune_rppa_df.index.name = "case_id"

print(f"Final immune protein matrix: {immune_rppa_df.shape[0]} samples, {immune_rppa_df.shape[1]} proteins")

# --------- 3. CALCULATE PROTEIN PATHWAY SCORES ---------
print("Calculating immune pathway scores...")

# Define pathway groups
pathway_groups = {
    'checkpoint_proteins': ['PD-L1', 'PDL1', 'CD274', 'PD-1', 'PD1', 'PDCD1', 'CTLA4', 'CTLA-4', 'LAG3', 'TIM3', 'TIGIT'],
    'tcell_signaling': ['CD3', 'CD4', 'CD8', 'CD28', 'CD80', 'CD86', 'ICOS'],
    'cytokine_signaling': ['TNF', 'IL2', 'IL6', 'IL10', 'IFNG', 'TGFB1'],
    'stat_signaling': ['STAT1', 'STAT3', 'STAT4', 'STAT5', 'STAT6'],
    'nfkb_signaling': ['NFKB', 'NFKB1', 'RELA', 'RELB'],
    'pi3k_akt_pathway': ['AKT', 'MTOR', 'P70S6K1', 'S6'],
    'mapk_pathway': ['ERK1', 'ERK2', 'P38', 'JNK'],
}

# Calculate pathway scores as mean of available proteins in each pathway
for pathway_name, pathway_proteins in pathway_groups.items():
    # Find available proteins for this pathway
    pathway_cols = []
    for protein in pathway_proteins:
        matching_cols = [col for col in immune_rppa_df.columns if protein.upper() in col.upper()]
        pathway_cols.extend(matching_cols)
    
    pathway_cols = list(set(pathway_cols))  # Remove duplicates
    
    if pathway_cols:
        immune_rppa_df[f'{pathway_name}_score'] = immune_rppa_df[pathway_cols].mean(axis=1)
        print(f"  {pathway_name}: {len(pathway_cols)} proteins")
    else:
        print(f"  {pathway_name}: No proteins found")

# --------- 4. OPTIONAL: ADD RESPONSE LABELS ---------
try:
    print("Adding responder/non-responder labels...")
    geo_meta = pd.read_csv(GEO_META_FILE)
    geo_meta["case_id"] = geo_meta["case_id"].astype(str)
    immune_rppa_df = immune_rppa_df.merge(geo_meta[["case_id", "response"]], on="case_id", how="left")
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
    immune_rppa_df.index = immune_rppa_df.index.str.replace(r"\.", "-", regex=True)
    wsi_df.index = wsi_df.index.str.replace(r"\.", "-", regex=True)

    immune_rppa_df = immune_rppa_df.merge(wsi_df, on="case_id", how="left")

    # Fill missing values in WSI columns with placeholder
    if "num_patches" in immune_rppa_df.columns:
        immune_rppa_df["num_patches"] = immune_rppa_df["num_patches"].fillna("No WSI Data")
    if "file" in immune_rppa_df.columns:
        immune_rppa_df["file"] = immune_rppa_df["file"].fillna("No WSI Data")

except FileNotFoundError:
    print("WARNING: WSI JSON file not found - skipping this step.")

# --------- 6. EXPORT ---------
output_file = "data/processed_data/rppa_immune_proteins_with_metadata.csv"
immune_rppa_df.to_csv(output_file)
print(f"SUCCESS: Done. Output saved to: {output_file}")
print(f"Final dataset shape: {immune_rppa_df.shape}")

# --------- 7. SUMMARY STATISTICS ---------
print("\nSummary Statistics:")
print(f"Number of samples: {immune_rppa_df.shape[0]}")

# Count different types of columns
protein_cols = [col for col in immune_rppa_df.columns if not col.endswith('_score') and col not in ['num_patches', 'file', 'response']]
pathway_cols = [col for col in immune_rppa_df.columns if col.endswith('_score')]

print(f"Number of immune proteins: {len(protein_cols)}")
print(f"Number of pathway scores: {len(pathway_cols)}")

if pathway_cols:
    print(f"\nPathway Scores:")
    for pathway in pathway_cols:
        score_mean = immune_rppa_df[pathway].mean()
        score_std = immune_rppa_df[pathway].std()
        print(f"  {pathway}: {score_mean:.3f} ± {score_std:.3f}")

if protein_cols:
    print(f"\nTop 5 most variable proteins:")
    protein_variance = immune_rppa_df[protein_cols].var().sort_values(ascending=False)
    for i, (protein, var) in enumerate(protein_variance.head(5).items()):
        print(f"  {i+1}. {protein}: variance = {var:.4f}")
