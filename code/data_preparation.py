from pyprojroot.here import here
import pandas as pd
import numpy as np
import squidpy as sq


# Set up root directory
base_path = pyprojroot.find_root(pyprojroot.has_dir(".git"))


# Load data 
vs_data = sq.read.visium(path = str(base_path) + '/data/raw',
                             counts_file = 'CytAssist_FFPE_Human_Skin_Melanoma_filtered_feature_bc_matrix.h5',
                             library_id = 'CytAssist_FFPE_Human_Skin_Melanoma',
                             load_images = True)
vs_data.var_names_make_unique()
IF_matrix_raw = pd.read_csv(here('data/raw/spatial/barcode_fluorescence_intensity.csv'))
tissue_position = pd.read_csv(here('data/raw/spatial/tissue_positions.csv'))


# Gene data processing
gene_counts = pd.DataFrame(vs_data.X.toarray(), columns=vs_data.var.index.tolist())
gene_counts = gene_counts.reset_index(drop=True)
gene_counts.index = vs_data.obs.index.tolist()
gene_matrix = vs_data.obs.copy()
gene_matrix.insert(0, 'barcode', gene_matrix.index.tolist())
gene_matrix = pd.concat([gene_matrix, gene_counts], axis=1)
# Save the processed gene data
# vs_data.write_h5ad(here('data/processed/gene_matrix.h5ad')) 


# IF data processing
IF_matrix = pd.merge(IF_matrix_raw, tissue_position, on=['barcode', 'in_tissue'], how='left')
IF_matrix = IF_matrix[IF_matrix['in_tissue'] == 1]
IF_matrix = IF_matrix[['barcode', 'channel3_mean', 'array_row', 'array_col']]
IF_matrix = IF_matrix.rename(columns={'channel3_mean': 'PCNA_IF'})
# Save the processed IF data
# IF_matrix.to_csv(here('data/processed/IF_matrix.csv'), index=True)


# Combine gene and IF data
protein_gene_matrix = pd.merge(gene_matrix, IF_matrix, on=['barcode', 'array_row', 'array_col'], how='left')
protein_gene_matrix['PCNA_IF'] = protein_gene_matrix['PCNA_IF'].fillna(0)
# Save the combined data
# protein_gene_matrix.to_csv(here('data/processed/protein_gene_matrix.csv'), index=True)
