#! /usr/bin/env python

'''
PCA analysis two star files.
'''

import pandas as pd
import numpy as np
from common.starparse import StarFile
from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA
import sys

def compare_dataframes_pca(df1, df2):
    # Get numeric columns
    numeric_cols = df1.select_dtypes(include=['float64', 'int64']).columns
    
    # Combine dataframes for scaling
    df1_numeric = df1[numeric_cols]
    df2_numeric = df2[numeric_cols]
    combined = pd.concat([df1_numeric, df2_numeric])
    
    # Scale the data
    scaler = StandardScaler()
    scaled_data = scaler.fit_transform(combined)
    
    # Perform PCA
    pca = PCA()
    pca_result = pca.fit_transform(scaled_data)
    
    # Calculate loadings
    loadings = pd.DataFrame(
        pca.components_.T * np.sqrt(pca.explained_variance_),
        columns=[f'PC{i+1}' for i in range(len(numeric_cols))],
        index=numeric_cols
    )
    
    # Split back into original groups
    df1_pca = pca_result[:len(df1)]
    df2_pca = pca_result[len(df1):]
    
    # Calculate mean difference along principal components
    mean_diff = np.abs(df1_pca.mean(axis=0) - df2_pca.mean(axis=0))
    
    # Get most significant features
    significant_features = pd.DataFrame({
        'Feature': numeric_cols,
        'Contribution': np.abs(loadings.iloc[:, 0])
    }).sort_values('Contribution', ascending=False)
    
    return {
        'explained_variance_ratio': pca.explained_variance_ratio_,
        'loadings': loadings,
        'significant_features': significant_features,
        'mean_difference': mean_diff
    }




starfile = StarFile( sys.argv[1])
df1 = starfile.particles_df

starfile = StarFile( sys.argv[2])
df2 = starfile.particles_df

# Usage
results = compare_dataframes_pca(df1, df2)
print("Most significant features:")
print(results['significant_features'])
print("\nExplained variance ratio:")
print(results['explained_variance_ratio'])