import pandas as pd
import numpy as np
import json
import subprocess
import sys
import os
from pathlib import Path

# --------- USER CONFIGURATION ---------
SCRIPTS_DIR = "scripts"
PROCESSED_DATA_DIR = "data/processed_data"
WSI_JSON_PATH = r"study_patient_mapping.json"

# TPM Normalization settings
APPLY_TPM_NORMALIZATION = True  # Set to False to skip TPM normalization
TPM_NORMALIZE_DATA_TYPES = ['rnaseq', 'mirna']  # Data types to apply TPM normalization to

# Processing scripts to run
PROCESSING_SCRIPTS = [
    "process_rnaseq.py",
    "process_mirna.py", 
    "process_methylation.py",
    "process_rppa.py",
    "process_copy_number.py"
]

# Expected output files from each script
EXPECTED_OUTPUTS = {
    "process_rnaseq.py": "rnaseq_immune_markers_with_metadata.csv",
    "process_mirna.py": "mirna_immune_markers_with_metadata.csv",
    "process_methylation.py": "methylation_immune_sites_with_metadata.csv",
    "process_rppa.py": "rppa_immune_proteins_with_metadata.csv",
    "process_copy_number.py": "copy_number_immune_features_with_metadata.csv"
}
# --------------------------------------

def run_processing_script(script_name):
    """Run a processing script and return success status"""
    script_path = os.path.join(SCRIPTS_DIR, script_name)
    
    if not os.path.exists(script_path):
        print(f"ERROR: Script not found: {script_path}")
        return False

    print(f"\nRunning {script_name}...")
    try:
        result = subprocess.run([sys.executable, script_path],
                              capture_output=True, text=True, timeout=300)

        if result.returncode == 0:
            print(f"SUCCESS: {script_name} completed successfully")
            print("Output:", result.stdout[-200:] if len(result.stdout) > 200 else result.stdout)
            return True
        else:
            print(f"ERROR: {script_name} failed with return code {result.returncode}")
            print("Error:", result.stderr)
            return False

    except subprocess.TimeoutExpired:
        print(f"TIMEOUT: {script_name} timed out after 5 minutes")
        return False
    except Exception as e:
        print(f"ERROR: Error running {script_name}: {e}")
        return False

def tpm_normalize(df, data_type='rnaseq', gene_length_col=None):
    """
    Apply TPM (Transcripts Per Million) normalization to expression data

    Parameters:
    - df: DataFrame with samples as rows and genes/features as columns
    - data_type: Type of data being normalized (for logging)
    - gene_length_col: Column name containing gene lengths (if available)

    Returns:
    - DataFrame with TPM-normalized values
    """
    print(f"  Applying TPM normalization to {data_type} data...")

    # Identify numeric columns (exclude metadata columns)
    metadata_cols = ['num_patches', 'file', 'response']
    numeric_cols = [col for col in df.columns if col not in metadata_cols]

    if len(numeric_cols) == 0:
        print(f"  WARNING: No numeric columns found for TPM normalization in {data_type}")
        return df

    # Create a copy to avoid modifying original data
    df_normalized = df.copy()

    # Extract numeric data for normalization
    numeric_data = df[numeric_cols].copy()

    # Handle negative values by shifting to positive range if needed
    min_val = numeric_data.min().min()
    if min_val < 0:
        print(f"  WARNING: Negative values detected in {data_type}. Shifting by {abs(min_val) + 1}")
        numeric_data = numeric_data + abs(min_val) + 1

    # If gene lengths are not provided, assume uniform length (standard for processed data)
    if gene_length_col is None or gene_length_col not in df.columns:
        # Standard TPM calculation without gene length correction
        # Step 1: Calculate RPK (Reads Per Kilobase) - assuming uniform gene length
        rpk_data = numeric_data  # No length normalization

        # Step 2: Calculate scaling factor (sum of RPK per sample / 1,000,000)
        scaling_factors = rpk_data.sum(axis=1) / 1e6

        # Step 3: Divide RPK by scaling factor to get TPM
        tpm_data = rpk_data.div(scaling_factors, axis=0)
    else:
        # TPM calculation with gene length correction
        gene_lengths = df[gene_length_col]

        # Step 1: Calculate RPK (Reads Per Kilobase)
        rpk_data = numeric_data.div(gene_lengths / 1000, axis=1)

        # Step 2: Calculate scaling factor
        scaling_factors = rpk_data.sum(axis=1) / 1e6

        # Step 3: Calculate TPM
        tpm_data = rpk_data.div(scaling_factors, axis=0)

    # Replace numeric columns with TPM-normalized values
    df_normalized[numeric_cols] = tpm_data

    # Verify TPM normalization (each sample should sum to ~1M)
    sample_sums = tpm_data.sum(axis=1)
    mean_sum = sample_sums.mean()
    print(f"  TPM normalization complete. Mean sample sum: {mean_sum:.0f} (target: 1,000,000)")

    return df_normalized

def load_processed_data(apply_tpm_normalization=True):
    """Load all processed data files with optional TPM normalization"""
    data_dict = {}

    for script, output_file in EXPECTED_OUTPUTS.items():
        file_path = os.path.join(PROCESSED_DATA_DIR, output_file)

        if os.path.exists(file_path):
            try:
                df = pd.read_csv(file_path, index_col=0)
                data_type = script.replace("process_", "").replace(".py", "")

                # Apply TPM normalization if requested
                if apply_tpm_normalization:
                    # Apply TPM normalization to specified data types
                    if data_type in TPM_NORMALIZE_DATA_TYPES:
                        df = tpm_normalize(df, data_type)
                    else:
                        print(f"  Skipping TPM normalization for {data_type} (not in configured types: {TPM_NORMALIZE_DATA_TYPES})")

                data_dict[data_type] = df
                print(f"SUCCESS: Loaded {data_type} data: {df.shape}")
            except Exception as e:
                print(f"ERROR: Error loading {output_file}: {e}")
        else:
            print(f"WARNING: File not found: {file_path}")

    return data_dict

def integrate_multiomics_data(data_dict):
    """Integrate multi-omics data by sample ID"""
    print("\nIntegrating multi-omics data...")

    if not data_dict:
        print("ERROR: No data to integrate!")
        return None
    
    # Find common samples across all data types
    all_samples = None
    for data_type, df in data_dict.items():
        samples = set(df.index)
        if all_samples is None:
            all_samples = samples
        else:
            all_samples = all_samples.intersection(samples)
        print(f"  {data_type}: {len(samples)} samples")
    
    print(f"  Common samples across all data types: {len(all_samples)}")
    
    if len(all_samples) == 0:
        print("ERROR: No common samples found across data types!")
        return None
    
    # Create integrated dataset
    integrated_data = {}
    metadata_cols = ['num_patches', 'file', 'response']
    
    for data_type, df in data_dict.items():
        # Filter to common samples
        df_filtered = df.loc[list(all_samples)]
        
        # Separate data columns from metadata
        data_cols = [col for col in df_filtered.columns if col not in metadata_cols]
        meta_cols = [col for col in df_filtered.columns if col in metadata_cols]
        
        # Store data with prefixed column names
        for col in data_cols:
            integrated_data[f"{data_type}_{col}"] = df_filtered[col]
        
        # Store metadata (will be the same across data types)
        if data_type == list(data_dict.keys())[0]:  # Only store metadata once
            for col in meta_cols:
                if col in df_filtered.columns:
                    integrated_data[col] = df_filtered[col]
    
    # Create integrated DataFrame
    integrated_df = pd.DataFrame(integrated_data, index=list(all_samples))
    integrated_df.index.name = "case_id"
    
    print(f"SUCCESS: Integrated dataset created: {integrated_df.shape}")
    return integrated_df

def calculate_multiomics_scores(integrated_df):
    """Calculate multi-omics immune scores"""
    print("\nCalculating multi-omics immune scores...")
    
    # Define data type prefixes
    data_types = ['rnaseq', 'mirna', 'methylation', 'rppa', 'copy_number']
    
    # Calculate per-data-type immune scores
    for data_type in data_types:
        data_cols = [col for col in integrated_df.columns if col.startswith(f"{data_type}_")]
        
        if data_cols:
            # Calculate mean z-score across features for each sample
            data_subset = integrated_df[data_cols]
            
            # Standardize features (z-score)
            data_standardized = (data_subset - data_subset.mean()) / data_subset.std()
            
            # Calculate immune score as mean of standardized values
            integrated_df[f"{data_type}_immune_score"] = data_standardized.mean(axis=1)
            
            print(f"  {data_type} immune score: {integrated_df[f'{data_type}_immune_score'].mean():.3f} ± {integrated_df[f'{data_type}_immune_score'].std():.3f}")
    
    # Calculate overall multi-omics immune score
    score_cols = [col for col in integrated_df.columns if col.endswith('_immune_score')]
    if score_cols:
        integrated_df['multiomics_immune_score'] = integrated_df[score_cols].mean(axis=1)
        print(f"  Multi-omics immune score: {integrated_df['multiomics_immune_score'].mean():.3f} ± {integrated_df['multiomics_immune_score'].std():.3f}")
    
    # Classify samples based on immune score
    if 'multiomics_immune_score' in integrated_df.columns:
        q33 = integrated_df['multiomics_immune_score'].quantile(0.33)
        q67 = integrated_df['multiomics_immune_score'].quantile(0.67)
        
        integrated_df['immune_classification'] = 'Intermediate'
        integrated_df.loc[integrated_df['multiomics_immune_score'] <= q33, 'immune_classification'] = 'Cold'
        integrated_df.loc[integrated_df['multiomics_immune_score'] >= q67, 'immune_classification'] = 'Hot'
        
        print(f"  Immune classification:")
        print(f"    Cold: {(integrated_df['immune_classification'] == 'Cold').sum()}")
        print(f"    Intermediate: {(integrated_df['immune_classification'] == 'Intermediate').sum()}")
        print(f"    Hot: {(integrated_df['immune_classification'] == 'Hot').sum()}")
    
    return integrated_df

def generate_summary_report(integrated_df, data_dict):
    """Generate a summary report of the integration"""
    print("\nGenerating summary report...")

    report = []
    report.append("# TCGA-ACC Multi-Omics Integration Report")
    report.append(f"Generated on: {pd.Timestamp.now().strftime('%Y-%m-%d %H:%M:%S')}")
    report.append("")

    # Processing configuration
    report.append("## Processing Configuration")
    report.append(f"- **TPM Normalization**: {'Enabled' if APPLY_TPM_NORMALIZATION else 'Disabled'}")
    if APPLY_TPM_NORMALIZATION:
        report.append(f"- **TPM Normalized Data Types**: {', '.join(TPM_NORMALIZE_DATA_TYPES)}")
    report.append("")
    
    # Data summary
    report.append("## Data Summary")
    report.append(f"- **Total samples in integrated dataset**: {integrated_df.shape[0]}")
    report.append(f"- **Total features**: {integrated_df.shape[1]}")
    report.append("")
    
    # Per data type summary
    report.append("## Data Type Summary")
    for data_type, df in data_dict.items():
        original_features = df.shape[1] - 3  # Subtract metadata columns
        integrated_features = len([col for col in integrated_df.columns if col.startswith(f"{data_type}_")])
        report.append(f"- **{data_type.upper()}**: {original_features} original features -> {integrated_features} integrated features")
    report.append("")
    
    # Immune scores summary
    if 'multiomics_immune_score' in integrated_df.columns:
        report.append("## Multi-Omics Immune Scores")
        score_cols = [col for col in integrated_df.columns if col.endswith('_immune_score')]
        for col in score_cols:
            mean_score = integrated_df[col].mean()
            std_score = integrated_df[col].std()
            report.append(f"- **{col}**: {mean_score:.3f} ± {std_score:.3f}")
        report.append("")
        
        # Classification summary
        report.append("## Immune Classification")
        for classification in ['Cold', 'Intermediate', 'Hot']:
            count = (integrated_df['immune_classification'] == classification).sum()
            pct = count / len(integrated_df) * 100
            report.append(f"- **{classification}**: {count} samples ({pct:.1f}%)")
        report.append("")
    
    # WSI data availability
    if 'num_patches' in integrated_df.columns:
        wsi_available = (integrated_df['num_patches'] != 'No WSI Data').sum()
        report.append("## WSI Data Availability")
        report.append(f"- **Samples with WSI data**: {wsi_available}/{len(integrated_df)} ({wsi_available/len(integrated_df)*100:.1f}%)")
        report.append("")
    
    # Save report
    report_text = "\n".join(report)
    report_file = os.path.join(PROCESSED_DATA_DIR, "multiomics_integration_report.md")
    with open(report_file, 'w') as f:
        f.write(report_text)
    
    print(f"SUCCESS: Summary report saved to: {report_file}")
    return report_text

def main():
    """Main integration pipeline"""
    print("Starting TCGA-ACC Multi-Omics Integration Pipeline")
    print("=" * 60)

    # Create output directory if it doesn't exist
    os.makedirs(PROCESSED_DATA_DIR, exist_ok=True)

    # Step 1: Run all processing scripts
    print("\nStep 1: Running individual processing scripts")
    success_count = 0
    
    for script in PROCESSING_SCRIPTS:
        if run_processing_script(script):
            success_count += 1
    
    print(f"\nSUCCESS: Completed {success_count}/{len(PROCESSING_SCRIPTS)} processing scripts")

    if success_count == 0:
        print("ERROR: No processing scripts completed successfully. Exiting.")
        return

    # Step 2: Load processed data with optional TPM normalization
    print(f"\nStep 2: Loading processed data (TPM normalization: {'enabled' if APPLY_TPM_NORMALIZATION else 'disabled'})")
    data_dict = load_processed_data(apply_tpm_normalization=APPLY_TPM_NORMALIZATION)

    if not data_dict:
        print("ERROR: No processed data found. Exiting.")
        return

    # Step 3: Integrate multi-omics data
    print("\nStep 3: Integrating multi-omics data")
    integrated_df = integrate_multiomics_data(data_dict)

    if integrated_df is None:
        print("ERROR: Integration failed. Exiting.")
        return

    # Step 4: Calculate multi-omics scores
    print("\nStep 4: Calculating multi-omics immune scores")
    integrated_df = calculate_multiomics_scores(integrated_df)

    # Step 5: Save integrated dataset
    print("\nStep 5: Saving integrated dataset")
    output_file = os.path.join(PROCESSED_DATA_DIR, "tcga_acc_multiomics_integrated.csv")
    integrated_df.to_csv(output_file)
    print(f"SUCCESS: Integrated dataset saved to: {output_file}")

    # Step 6: Generate summary report
    print("\nStep 6: Generating summary report")
    generate_summary_report(integrated_df, data_dict)

    print("\nSUCCESS: Multi-omics integration pipeline completed successfully!")
    print("=" * 60)
    print(f"Final integrated dataset: {integrated_df.shape[0]} samples × {integrated_df.shape[1]} features")

if __name__ == "__main__":
    main()
