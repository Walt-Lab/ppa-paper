#!/usr/bin/env python3

import os
import pandas as pd
import sys
import matplotlib.pyplot as plt
from pathlib import Path

sys.path.append(str(Path().absolute().parent / "src"))

from olink_fractionation import analyze_fractionation
from raw_data_preprocessing import (
    clean_up_raw_data,
    plot_protein_fractionation,
    calculate_fractionation_scores,
)

# Get data and output directories relative to notebook
data_dir = Path().absolute().parent / "data"
output_dir = Path().absolute().parent / "outputs"

# Define paths relative to those directories
assay_list_path = data_dir / "231220_ht_panel_assay_list.xlsx"
brain_rna_seq_raw_path = data_dir / "240411_brain_rna_seq_raw.csv"
output_directory = output_dir / "ht_output"
plate_layout_path = data_dir / "231204_Walt_Olink_HT_Plate.xlsx"
raw_data = data_dir / "240214_Walt_Olink_HT_Raw.parquet"
uniprot_fasta_database = data_dir / "uniprot_fasta_database.gz"
gtex_path = data_dir / "GTEx_Analysis_v10_RNASeQCv2.4.2_gene_median_tpm.gct.gz"
lod_path = data_dir / "Explore HT_Fixed LOD_2024-06-19.csv"

fraction_list_path = "fraction_list.txt"
fractionation_scores_path = "fractionation_scores.txt"

# Create a tidy dataframe from the raw data file.
tidy_data = clean_up_raw_data(raw_data, plate_layout_path)

pattern1 = {
    "high_fraction_options" : [
        ["8", "9", "10"]
    ],
    "low_fraction_options" : [
        ["11", "12", "13", "14", "15"]
    ]
}
pattern2 = {
    "high_fraction_options" : [
        ["8", "9", "10", "14", "15"]
    ],
    "low_fraction_options" : [
        ["11", "12", "13"]
    ]
}
pattern3 = {
    "high_fraction_options" : [
        ["11", "12", "13", "14", "15"]
    ],
    "low_fraction_options" : [
        ["6", "7", "8", "9", "10"]
    ]
}

for pattern in [pattern1, pattern2, pattern3]:
    fraction_calcs = {}
    fraction_list = []
    fraction_dict = {}

    for high_fraction in pattern["high_fraction_options"]:
        for low_fraction in pattern["low_fraction_options"]:
            print(f"{high_fraction} vs {low_fraction}")
            fraction_list.append((high_fraction, low_fraction))
            results = analyze_fractionation(
                tidy_data, 
                high_fraction,
                low_fraction,
                sample_health="healthy",
                mean_median_individual="individual_median"
            )
            fraction_dict[f"{high_fraction} vs {low_fraction}"] = results
            
            # Calculate scores for ALL proteins at once
            scores = calculate_fractionation_scores(results, tidy_data, high_fraction, low_fraction)
            
            # Create dictionary mapping protein IDs to scores
            protein_score_dict = dict(zip(results, scores))
            fraction_calcs[f"{high_fraction} vs {low_fraction}"] = protein_score_dict


    # # Make file of the high/low fraction combos
    # with open(fraction_list_path, 'w') as f:
    #     for item in fraction_list:
    #         f.write(f"{item}\n")

    # # Make file of the fractionation scores for each high/low fraction combo
    # with open(fractionation_scores_path, 'w') as f:
    #     for key, value in fraction_calcs.items():
    #         f.write(f"{key} {value}\n\n")

    # Make graphs for each target, in folders named after the high/low fraction pairs
    outermost_folder_name = "pattern1" if pattern == pattern1 else "pattern2" if pattern == pattern2 else "pattern3"
    os.makedirs(outermost_folder_name, exist_ok=True)

    for key, value_list in fraction_dict.items():
        folder_name = key  
        safe_folder_name = folder_name.replace('[', '').replace(']', '').replace("'", '')
        full_path = os.path.join(outermost_folder_name, safe_folder_name)
        os.makedirs(full_path, exist_ok=True)
        
        # Create plots for ALL proteins
        for i, value in enumerate(value_list):
            filename = f"{full_path}/{value}.png"
            ax = plot_protein_fractionation(tidy_data, value)
            
            fig = ax.get_figure()
            fig.savefig(filename, dpi=300, bbox_inches='tight')
            plt.close(fig)
        
        # Create a separate "top5" folder for the top 5 proteins by score
        top5_path = os.path.join(full_path, "top5")
        os.makedirs(top5_path, exist_ok=True)
        
        # Get top 5 proteins by score
        protein_scores = fraction_calcs[key]
        top5_proteins = sorted(protein_scores.items(), key=lambda x: x[1], reverse=True)[:5]
        
        # Create plots for top 5
        for protein_id, score in top5_proteins:
            filename = f"{top5_path}/{protein_id}_score_{score:.2f}.png"
            ax = plot_protein_fractionation(tidy_data, protein_id)
            
            fig = ax.get_figure()
            fig.savefig(filename, dpi=300, bbox_inches='tight')
            plt.close(fig)

# Should give an output like:
# pattern1/
#   8, 9, 10 vs 11, 12, 13, 14, 15/
#     A0JNW5.png
#     A1L4H1.png
#     ... (all proteins)
#     top5/
#       P12345_score_2.34.png
#       Q67890_score_2.15.png
#       ... (only top 5)
# pattern2/
#   ...