import pyprojroot
import pandas as pd
import numpy as np
from scipy.spatial import distance
import warnings


# Set up root directory
base_path = pyprojroot.find_root(pyprojroot.has_dir(".git"))


# Define the function
def calculate_RipleyK(x, y, protein_weights, gene_matrix, r_min=1, r_max=None, r_step=1, 
                        area=None, weighted=True, n_genes=None, return_df=True, gene_names=None):
    """
    Calculate the Ripley's K function for spatial point pattern analysis.

    Parameters:
    - x, y (numpy.array): Coordinates of points.
    - protein_weights (numpy.array): Array of protein weights corresponding to each point.
    - gene_matrix (numpy.array): Matrix of gene expression data; rows correspond to points, columns to genes.
    - r_min, r_max, r_step (float, optional): Min, max, and step for the radii to be used in the analysis. Defaults r_min to 1, r_max to the smaller of the maximum values of x and y coordinates, and r_step to 1, if None.
    - area (float, optional): Area over which points are distributed. Defaults to bounding rectangle if None.
    - weighted (bool): If True, calculated mark (both gene and protein )weighted Ripley's K, else calculated unweighted Ripley's K.
    - n_genes (int, optional): Number of genes to analyze. Defaults to the number of columns in gene_matrix if None.
    - return_df (bool): If True, returns results as a pandas DataFrame, otherwise as a numpy array.
    - gene_names (list of str, optional): Names of genes; used for DataFrame column labeling.

    Returns:
    - Numpy array of Ripley's K values or a Pandas DataFrame with coordinates and gene names for each gene and target protein across specified radii.
    """

    # Prepare point coordinates and calculate pairwise distances
    points = np.vstack((x, y)).T
    dist_matrix = distance.cdist(points, points)

    # Unweighted: convert non-zero weights to 1 (binary presence)
    if not weighted:
        protein_weights[np.where(protein_weights != 0)] = 1
        gene_matrix[np.where(gene_matrix != 0)] = 1

    # Set default values for radii if not specified
    if r_max is None:
        r_max = max(x) if max(x) <= max(y) else max(y)

    # Define the radii to analyze
    radii = np.arange(r_min, r_max + 1, r_step)

    # Set default area if not specified
    if area is None:
        area = max(x) * max(y)

    # Determine number of genes if not specified
    if n_genes is None:
        n_genes = gene_matrix.shape[1]

    # Calculate the K values across all radii for each gene
    k_matrix = []
    for gene in range(n_genes):
        gene_weights = gene_matrix[:, gene]
        weight_products = np.outer(gene_weights, protein_weights)

        gene_weighted_sum = weight_products.sum(axis=1)

        k_values = []
        for radius in radii:
            within_radius = dist_matrix <= radius
            gene_weighted_sum_mask = (weight_products * within_radius).sum(axis=1)
            with warnings.catch_warnings():
                warnings.simplefilter("ignore", category=RuntimeWarning)
                k = np.nansum((gene_weighted_sum_mask / gene_weighted_sum)) * area
            k_values.append(k)

        k_matrix.append(k_values)

    # Return the results either as a DataFrame or as a raw matrix
    if return_df:
        if gene_names is None:
            print("\nGene names are not provided. Using default names.\n")
            gene_names = [f'Gene_{i+1}' for i in range(n_genes)]
        k_matrix_df = pd.DataFrame(k_matrix, index=gene_names).T
        k_matrix_df.insert(0, 'radius', radii)
        return k_matrix_df
    else:
        return k_matrix

