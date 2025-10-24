#!/usr/bin/env python
# coding: utf-8

# In[ ]:


import pandas as pd
import numpy as np
import matplotlib
from matplotlib import pylab as plt
import re
import os 
from sklearn.compose import ColumnTransformer
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import StandardScaler, OneHotEncoder, OrdinalEncoder, MinMaxScaler
from sklearn.model_selection import train_test_split


# read & merge all count files 
data_path = '/oscar/data/larschan/etaner1/DATA1030_ML/count_matrix_120'
file_list = [f for f in os.listdir(data_path) if f.endswith('.tsv')]
# print(file_list)

# read in metadata 
meta = pd.read_csv('/oscar/data/larschan/etaner1/DATA1030_ML/tsv/experiment_report_2025_120.tsv', sep='\t', skiprows=1)

# inspect columns 
print(meta.columns.tolist())
# a lot of these don't matter or don't add anything 
# NOTE!! sex is missing and needs to be appended 

# Keep only relevant columns
# Keep only columns of interest
meta = meta[[
    'Accession', 'Biosample accession', 'Biosample summary',
    'Biosample age', 'Biological replicate', 'Technical replicate'
]]


# In[ ]:


# parse & later append 'sex' from Biosample summary
meta['Sex'] = meta['Biosample summary'].str.extract(r'(?i)\b(female|male)\b', expand=False)

# parse & clean "Age" into categories (for now, might change later?)
def clean_age(age_string):
    if pd.isna(age_string):
        return None
    text = str(age_string).lower()
    if "90" in text:
        return "90+"
    elif "80" in text:
        return "80s"
    elif "70" in text:
        return "70s"
    else:
        num = re.findall(r'\d+', text)
        return num[0] if num else None

meta['Age_Category'] = meta['Biosample age'].apply(clean_age)

# remove any white space & ensure consistent formatting
meta['Sex'] = meta['Sex'].str.lower()
meta = meta.drop_duplicates(subset='Accession')

print(meta[['Accession', 'Sex', 'Age_Category', 'Biosample accession']].head())


# In[ ]:


# read one file as a dataframe
def read_count_file(filename):
    df = pd.read_csv(os.path.join(data_path, filename), sep='\t', usecols=['gene_id', 'expected_count'])
    df = df.set_index('gene_id')
    sample_id = filename.replace('.tsv', '')
    df.columns = [sample_id]
    return df


# In[ ]:


# merge via concatenation (not iterative .merge)!!
dfs = [read_count_file(f) for f in file_list]
merged = pd.concat(dfs, axis=1)

print(f"Merged expression matrix shape: {merged.shape}")  # genes × samples


# In[ ]:


# aligning metadata w exp matrix 
expr_T = merged.T

# need to reset index for merging by sample accession (Accession)
expr_T.index.name = 'Accession'
expr_T = expr_T.reset_index()

# merge metadata with expression data
df_full = expr_T.merge(meta, on='Accession', how='inner')
df_full = df_full.set_index('Accession')

print(f"Combined dataset shape: {df_full.shape}")  # rows = samples, columns = genes + metadata


# In[ ]:


# splitting by group 
X = df_full.drop(columns=['Sex'])
y = df_full['Sex']
groups = df_full['Biosample accession']


# In[ ]:


# Group-aware split (80/10/10)
gss_outer = GroupShuffleSplit(n_splits=1, test_size=0.1, random_state=42)
train_val_idx, test_idx = next(gss_outer.split(X, y, groups=groups))

X_trainval, X_test = X.iloc[train_val_idx], X.iloc[test_idx]
y_trainval, y_test = y.iloc[train_val_idx], y.iloc[test_idx]
groups_trainval = groups.iloc[train_val_idx]


gss_inner = GroupShuffleSplit(n_splits=1, test_size=0.111, random_state=42)
train_idx, val_idx = next(gss_inner.split(X_trainval, y_trainval, groups=groups_trainval))

X_train, X_val = X_trainval.iloc[train_idx], X_trainval.iloc[val_idx]
y_train, y_val = y_trainval.iloc[train_idx], y_trainval.iloc[val_idx]

print(f"Train: {X_train.shape[0]} | Val: {X_val.shape[0]} | Test: {X_test.shape[0]}")


# In[ ]:


# preprocessing
# (Fit on Train Only)
numeric_features = [col for col in X_train.columns if col.startswith('ENSG') or col.isnumeric()]
categorical_features = ['Age_Category']


# In[ ]:


# column transformation: scale numeric, one-hot encode categorical
preprocessor = ColumnTransformer(
    transformers=[
        ('num', StandardScaler(), numeric_features),
        ('cat', OneHotEncoder(handle_unknown='ignore'), categorical_features)
    ])

# build preprocessing pipeline
pipeline = Pipeline(steps=[('preprocessor', preprocessor)])


# In[ ]:


# fit ONLY on training data
pipeline.fit(X_train)

# transform train, val, test sets
X_train_processed = pipeline.transform(X_train)
X_val_processed = pipeline.transform(X_val)
X_test_processed = pipeline.transform(X_test)


# to confirm correct splitting 
print(f"Train matrix transformed: {X_train_processed.shape}")
print(f"Validation matrix transformed: {X_val_processed.shape}")
print(f"Test matrix transformed: {X_test_processed.shape}")

