import sys
import pyprojroot
from pyprojroot.here import here
import pandas as pd
import numpy as np
from scipy.spatial import distance
import matplotlib.pyplot as plt



# Set up root directory
base_path = pyprojroot.find_root(pyprojroot.has_dir(".git"))

# Import the calculate_RipleyK and RipleyK_statistic from the scripts 
sys.path.append(str(base_path.joinpath('code')))
from Ripleys_K_function import calculate_RipleyK
from Ripleys_K_based_statistics import RipleyK_statistic
sys.path.remove(str(base_path.joinpath('code')))



############################
# 1 - Ripley's K function
############################

# calculate_RipleyK(x, y, protein_weights, gene_matrix, r_min=1, r_max=None, r_step=1, 
#                     area=None, weighted=True, n_genes=None, return_df=True, gene_names=None)


# Example of how to use the function using demo data
## Load demo data
protein_gene_matrix = pd.read_csv(here('data/processed/protein_gene_matrix.csv'))

## Prepare the input of the function
x = protein_gene_matrix['array_row'].values
y = protein_gene_matrix['array_col'].values
protein_weights = protein_gene_matrix["PCNA_IF"].values
gene_matrix = protein_gene_matrix.iloc[:, 4:-1].values

## Calculate Ripley's K function
k_weighted_result = calculate_RipleyK(x = x, y = y, 
                               protein_weights=protein_weights, gene_matrix=gene_matrix, 
                               weighted=True, return_df=True, n_genes=10)

print("\nWeighted Ripley's K function results:\n")
print(k_weighted_result)

k_unweight_result = calculate_RipleyK(x = x, y = y, 
                               protein_weights=protein_weights, gene_matrix=gene_matrix, 
                               weighted=False, return_df=True, n_genes=10)

print("\nUnweighted Ripley's K function results:\n")
print(k_unweight_result)


## Save the demo results
k_weighted_result.to_csv(here('output/weighted_Ripleys_K_demo_results.csv'), index=False)
k_unweight_result.to_csv(here('output/unweighted_Ripleys_K_demo_results.csv'), index=False)




###################################
# 2 - Ripley's K based statistic
###################################

# RipleyK_statistic(weighted_RipleyK, unweighted_RipleyK, radius_low=1, radius_up=None)


# Example of how to use the function using demo data
deduction_df, statistic_df = RipleyK_statistic(k_weighted_result, k_unweight_result, radius_low=1, radius_up=20)
print('\nRipley\'s K based statistic:\n')
print(statistic_df)

## Plot the difference between weighted and unweighted Ripley's K functions
plt.figure(figsize=(10, 6))
for i in range(10):
    plt.plot(deduction_df["radius"].values, deduction_df.iloc[:,i+1].values, label=deduction_df.columns[i+1], marker='o', markersize=1.5, alpha=0.7)
plt.xlabel('Radius (Spot)')
plt.ylabel('Difference between Weighted and Unweighted Ripley\'s K Functions')
plt.title('Unweighted Ripley\'s K Function in Melanoma (Demo) Dataset')
plt.legend(loc='center left', bbox_to_anchor=(1, 0.5))
plt.show()