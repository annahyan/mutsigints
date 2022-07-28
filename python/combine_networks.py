#!/usr/bin/env python
# coding: utf-8
# you can run it like
# python3 common_subnets_types.py input_dir [output_dir]
# In[2]:


import os
import collections
import pandas as pd
import numpy as np
import copy as cp
import itertools
import hashlib
import joblib
import inspect
import sys
import glob

input_dir = sys.argv[1]
if (len(sys.argv) > 2):
    output_dir = sys.argv[2] + "/"
else:
    output_dir = input_dir + "/"

# In[4]:

# get_ipython().run_line_magic('load_ext', 'rpy2.ipython')


# ### Reading in signature interactions

# In[5]:


all_matrices = dict()
tissues = list()


# In[5]:


for filename in os.listdir(input_dir):
    
    if not filename.endswith(".txt"):
        continue  

    tissue = filename.split(".")[0]
    tissues.append(tissue)
    
    all_matrices[tissue] = pd.read_table(input_dir + "/" + filename)


# In[6]:


tissues = set(tissues)


# ### Functions

# In[7]:


def intersect_cor_matrices(mat1, mat2):
    
    mat1 = mat1.copy()
    mat2 = mat2.copy()
    
    np.fill_diagonal(mat1.values, 0)
    np.fill_diagonal(mat2.values, 0)
    
    mat1 = mat1.loc[ ( mat1 != 0 ).any(axis = 0), : ]
    mat1 = mat1.loc[ :, ( mat1 != 0 ).any(axis = 0) ]
    
    mat2 = mat2.loc[ ( mat2 != 0 ).any(axis = 0), : ]
    mat2 = mat2.loc[ :, ( mat2 != 0 ).any(axis = 0) ]
    
    common_sigs = pd.Series(mat1.columns.intersection(mat2.columns))
    
    mat1 = mat1.loc[ common_sigs, common_sigs ] 
    mat2 = mat2.loc[ common_sigs, common_sigs ]
    
    out = ((mat1 == mat2) & (mat1 != 0)) * mat1

    out = out.loc[ ( out != 0 ).any(axis = 0), : ]
    out = out.loc[ :, ( out != 0 ).any(axis = 0) ]
        
    return(out)


# In[8]:


def break_group(comb_int_obj, tissue_set, matrix_combo):
    
    tissue_subset = tissue_set - comb_int_obj.tissues
    
    tissue_subset = list(tissue_subset) 
    
    for tissue_i in range(len(tissue_subset)):
        tissue = tissue_subset[tissue_i]
        
        mat_int = intersect_cor_matrices(comb_int_obj.sig_adj, all_matrices[tissue])
        
        summed_tissues = cp.deepcopy(comb_int_obj.tissues)
        summed_tissues.add(tissue)
        
#         print(summed_tissues)
#         print(type(summed_tissues))
        sig_clust_el = SignatureCluster(summed_tissues, mat_int)
        
        mat_int_hash = joblib.hash(mat_int)
    
        if mat_int_hash in matrix_combo: 
            try:
                matrix_combo[mat_int_hash] = matrix_combo[mat_int_hash].merge_with(sig_clust_el)
            except AttributeError:
                print("matrix combo:")
                print(matrix_combo[mat_int_hash])
                print("\n###\nsig_clust_el:")
                print(sig_clust_el)
                print("count = %d" % count )
                break
        else:
            matrix_combo[mat_int_hash] = cp.deepcopy(sig_clust_el)
            
            break_group(matrix_combo[mat_int_hash], set(tissue_subset[tissue_i:]), matrix_combo)


# In[9]:


def rescue_code(function):
    import inspect
    get_ipython().set_next_input("".join(inspect.getsourcelines(function)[0]))


# In[10]:


# rescue_code(compare_cor_matrices)


# ### Class definition

# In[11]:


class SignatureCluster:
   
    def __init__(self, tissues, sig_adj):
        self.tissues = set(tissues)
        self.sig_adj = sig_adj
        
    def __repr__(self):
        print("Tissues:")
        print(self.tissues)
        print("Signature interaction matrix:")
        print(self.sig_adj)
        return("")

    def __str__(self):
        print("Tissues:")
        print(self.tissues)
        print("Signature interaction matrix:")
        print(self.sig_adj)
        return("")
        
    def merge_with(self, other):
        
        hash_self = joblib.hash(self.sig_adj)
        hash_other = joblib.hash(other.sig_adj)
        
        if (hash_self != hash_other):
            raise Exception("Signature adjacency matrices between two objects has to be the same.")
            
        self.tissues.update(other.tissues)
        
        return(self)
    
        
    def __deepcopy__(self, memo):
        copy = type(self)(set(), None)
        memo[id(self)] = copy
        copy.tissues = cp.deepcopy(self.tissues, memo)
        copy.sig_adj = cp.deepcopy(self.sig_adj, memo)
        return(copy)


# ### Real code

# Generation of all pairwise combinations of tissues

# In[40]:


tissue_combinations = itertools.combinations(tissues, 2)


# Creating the SingatureCluster elements which are the largest common segments

# In[41]:


pair_interactions = list()

for i in tissue_combinations:
    
    tissue_pairs = i
    sig_intersection = intersect_cor_matrices(all_matrices[tissue_pairs[0]], all_matrices[tissue_pairs[1]])
    
    sig_cluster_el = SignatureCluster(tissue_pairs, sig_intersection)
    
    pair_interactions.append(sig_cluster_el)
    


# In[42]:


pair_interactions[1]


# Creating matrix_dict and adding networks of only pairwise interactions to the dictionary.

# In[54]:


matrix_dict = {}
count = 0
for pair_int in pair_interactions:
    
    count += 1
    
    sig_adj_hash = joblib.hash(pair_int.sig_adj)
    
    if sig_adj_hash in matrix_dict: 
        try:
            matrix_dict[sig_adj_hash] = matrix_dict[sig_adj_hash].merge_with(pair_int)
        except AttributeError:
            print("matrix dict:")
            print(matrix_dict[sig_adj_hash])
            print("\n###\npairint:")
            print(pair_int)
            print("count = %d" % count )
            break
    else:
        matrix_dict[sig_adj_hash] = cp.deepcopy(pair_int)
    


# In[56]:


len(matrix_dict)


# In[16]:


# for ii in matrix_dict:
#     print(matrix_dict[ii].sig_adj.shape)


# Collecting larger matrices after only pairwise intersections

# In[50]:


pack_of_combs = list()

for key in matrix_dict:
    if matrix_dict[key].sig_adj.shape[0] > 2:
        pack_of_combs.append(matrix_dict[key])


# In[55]:


len(pack_of_combs)


# Breaking down larger matrices with iterative intersections

# In[57]:


for sig_comb in pack_of_combs:
    break_group(sig_comb, tissues, matrix_dict)


# Checking if the identified tissues indeed contain the singature interactions

# In[58]:


for key in matrix_dict:
    
    obj = matrix_dict[key]
    
    for tissue in obj.tissues:
        mat_intersect = intersect_cor_matrices(obj.sig_adj, all_matrices[tissue])
        
        if not mat_intersect.equals(obj.sig_adj):
            print ("Tissue: ", tissue, "doesn't contain")
            print("Signatures:", obj.sig_adj)


# Checking if these signatures appear only among the tissues identified

# In[77]:


for key in matrix_dict:
    
    obj = matrix_dict[key]
        
    other_tissues = list(set(tissues) - obj.tissues)
    for tissue in other_tissues:
        mat_intersect = intersect_cor_matrices(obj.sig_adj, all_matrices[tissue])
        
        if mat_intersect.equals(obj.sig_adj):
            print(key)
            print ("Tissue: ", tissue, "contains")
            print("Signatures:", obj.sig_adj)    
    


# Writing tissues and signatures into text files

# In[20]:

output_dir = output_dir + '/matrix_dict_out/'
if not os.path.exists(output_dir):
    os.makedirs(output_dir)
else:
    files = glob.glob(output_dir + '*.txt')
    for f in files:
        os.remove(f)

for key in matrix_dict:
    
    el_tissues = matrix_dict[key].tissues
    el_sigs = matrix_dict[key].sig_adj
    
    if el_sigs.shape[0] == 0:
        continue
    
    f = open(output_dir + key + "_tissues.txt", "w")
    for ts in el_tissues:
        f.write(ts + '\n')
    f.close()
    
    el_sigs.to_csv(output_dir + key + "_sig_adj.txt")

