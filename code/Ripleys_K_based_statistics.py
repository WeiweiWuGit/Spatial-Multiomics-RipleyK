import sys
import pyprojroot
import pandas as pd
import numpy as np


# Set up root directory
base_path = pyprojroot.find_root(pyprojroot.has_dir(".git"))

# Import the calculate_RipleyK from the script before the new function
sys.path.append(str(base_path.joinpath('code')))
from Ripleys_K_function import calculate_RipleyK
sys.path.remove(str(base_path.joinpath('code')))


def RipleyK_statistic(weighted_RipleyK, unweighted_RipleyK, radius_low=1, radius_up=None):
    """
    Calculate the difference between weighted and unweighted Ripley's K functions and provide Ripley's K based statistic.

    Parameters:
    - weighted_RipleyK (pd.DataFrame): DataFrame containing the mark-weighted Ripley's K results.
    - unweighted_RipleyK (pd.DataFrame): DataFrame containing the unweighted Ripley's K results.
    - radius_low (int): Lower bound of the radius range for Ripley's K based statistic. Defaults to 1.
    - radius_up (int): Upper bound of the radius range for Ripley's K based statistic. Defaults to maximum radius.

    Returns:
    - tuple: Contains two DataFrames, deduction_df and statistic_df.
    """
    # Calculate differences between weighted and unweighted Ripley's K
    deduction_df = weighted_RipleyK.set_index('radius').subtract(unweighted_RipleyK.set_index('radius')).reset_index()
    
    # Set default value for radius_up if not specified
    if radius_up is not None:
        radius_up = deduction_df['radius'].max()

    # Select rows within the specified radius range
    select_df = deduction_df[(deduction_df['radius'] >= radius_low) & (deduction_df['radius'] <= radius_up)]

    # Calculate statistical summary
    statistic_df = pd.DataFrame(select_df.drop('radius', axis=1).mean(axis=0)).T

    return deduction_df, statistic_df


