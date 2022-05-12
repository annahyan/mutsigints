# module load IPython/7.15.0-foss-2020a-Python-3.8.2
# module load Keras/2.3.1-foss-2020a-Python-3.8.2

import logging

# from scipy.sparse import csc_matrix
import pandas as pd
import numpy as np
from sklearn.decomposition import TruncatedSVD
from sklearn.model_selection import train_test_split

# from sklearn.datasets import make_classification
from sklearn.ensemble import ExtraTreesClassifier
from sklearn.tree import DecisionTreeClassifier
from sklearn.metrics import ConfusionMatrixDisplay, confusion_matrix
from sklearn.ensemble import RandomForestClassifier

import matplotlib.pyplot as plt

# from sklearn.linear_model import LinearRegression


### transposing all the files
# for file in *.htseq.fpkm.dat; do echo $file; transpose.awk $file > ${file%.htseq.fpkm.dat}".transpose.dat"; done

#### concating all the files
# for file in `ls *.transpose.dat`; do echo $file;
# tail -n +2 $file >> tcga_all_exps.htseq.fpkm.dat;
# done


# logging.info('Performing dimensionality reduction on modality 2 values...')
# embedder_mod2 = TruncatedSVD(n_components=50)
# mod2_pca = embedder_mod2.fit_transform(input_train_mod2.X)

# # split dimred back up
# X_train = mod1_pca[input_train.obs['group'] == 'train']
# X_test = mod1_pca[input_train.obs['group'] == 'test']
# y_train = mod2_pca

# assert len(X_train) + len(X_test) == len(mod1_pca)

# spp_table = pd.read_excel('1-s2.0-S1074761318301213-mmc2.xlsx')

## in R
## library(openxlsx)
## a = read.xlsx('1-s2.0-S1074761318301213-mmc2.xlsx')
## a['Immune.Subtype'][is.na(a['Immune.Subtype'])] = "CNA"
## write.table(a, file= "immune_subtypes.tsv", sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)

## cat RNA_seq/tcga_all_exps.htseq.fpkm.dat | cut -f1 | cut -c1-12 |  sort | uniq -c | awk '$1 == 1 {print $2}'> trimmed_unique_sample_names.dat

## awk 'FNR == NR {a[$1]++; next}{if (a[$1]) print $0;}' trimmed_unique_sample_names.dat immune_subtypes.tsv > immune_subtypes_matched.tsv


## awk 'FNR==NR{a[$1]++; next} {if (a[substr($1, 0, 12)]) print $0}' immune_subtypes_matched.tsv  RNA_seq/tcga_all_exps.htseq.fpkm.dat > tcga_all_exps.subset.dat

## awk 'FNR==NR {imm[$1] = $0; next} substr($1, 0, 12) in imm {print imm[substr($1, 0, 12)]}' immune_subtypes_matched.tsv tcga_all_exps.subset.dat > imm_subtypes_sorted.tsv


## cd ../; comm -12 <(sort PCAWG_RNA_seq/tophat_ensg_names.dat) <(sort RNA_seq/TCGA_ensg_names.dat) > common_tcga_pcawg_genes.dat


all_exps = pd.read_csv("tcga_all_exps.subset.dat", sep = "\t", index_col = 0, header = None) 
tcga_gene_names = pd.read_csv("RNA_seq/TCGA_ensg_names.dat", sep = "\t", header = None, index_col = 0)

all_exps.columns = tcga_gene_names.index

common_gene_names = pd.read_csv("common_tcga_pcawg_genes.dat", sep = "\t", header = None, index_col = 0)


all_exps_sub = all_exps[common_gene_names.index]


exps_np = all_exps_sub.to_numpy() 
 
# embedder_svd = TruncatedSVD(n_components=200) 
# mod_exps = embedder_svd.fit_transform(exps_np)


exps_log = np.log(exps_np + 1)

with open('npy_files/exps_log.npy', 'wb') as f:
    np.save(f, exps_log)

# mod_log_exps = embedder_svd.fit_transform(exps_log)



unique_samples = pd.read_csv('trimmed_unique_sample_names.dat', header = None)

spp_table = pd.read_csv('imm_subtypes_sorted.tsv', sep = '\t', header = None)
# spp_table = spp_table.set_index('TCGA.Participant.Barcode')

#### Setting up X and y ######

# X = mod_exps
# X_log = mod_log_exps
y = spp_table[2].to_numpy()

# X_train, X_test, y_train, y_test = train_test_split(
#     X, y, test_size=0.3, random_state=1121218
# )


# X_log_train, X_log_test, y_log_train, y_log_test = train_test_split(
#     X_log, y, test_size=0.3, random_state=1121218
# )


# ####

# etc = ExtraTreesClassifier()
# _ = etc.fit(X_train, y_train)
# y_pred = etc.predict(X_test)

# # Plot confusion matrix

# fig, ax = plt.subplots(figsize=(8, 5))
# cmp = ConfusionMatrixDisplay(
#     confusion_matrix(y_test, y_pred),
#     display_labels=['C1', 'C2', 'C3', 'C4', 'C5', 'C6', 'CNA'],
# )

# cmp.plot(ax=ax)

# plt.savefig('figures/extra_treesclassifier.png')
# # plt.show();


# etc = ExtraTreesClassifier()
# _ = etc.fit(X_log_train, y_log_train)
# y_log_pred = etc.predict(X_log_test)

# # Plot confusion matrix

# fig, ax = plt.subplots(figsize=(8, 5))
# cmp = ConfusionMatrixDisplay(
#     confusion_matrix(y_log_test, y_log_pred),
#     display_labels=['C1', 'C2', 'C3', 'C4', 'C5', 'C6', 'CNA'],
# )

# cmp.plot(ax=ax)

# plt.savefig('figures/log_extra_treesclassifier.png')
# # plt.show();



# #### DecsionTreeClassifier ####


# for max_d in np.arange(4, 12):
#     print(max_d)
#     clf = DecisionTreeClassifier(max_depth = max_d)
#     clf = clf.fit(X_train, y_train)
    
#     y_pred = clf.predict(X_test)
    
    
#     fig, ax = plt.subplots(figsize=(8, 5))
#     cmp = ConfusionMatrixDisplay(
#         confusion_matrix(y_test, y_pred),
#         display_labels=['C1', 'C2', 'C3', 'C4', 'C5', 'C6', 'CNA'],
#     )
    
#     cmp.plot(ax=ax)
    
#     plt.savefig('figures/decisiontreeclassifier_md_' + str(max_d) + '.png')





# #### Hyper-parameter tuning

# import numpy as np
from scipy.stats import randint
from sklearn.experimental import enable_halving_search_cv  # noqa
from sklearn.model_selection import HalvingRandomSearchCV
from sklearn.ensemble import RandomForestClassifier
# from sklearn.datasets import make_classification

rng = np.random.RandomState(0)

# X, y = make_classification(n_samples=700, random_state=rng)

clf = RandomForestClassifier(random_state=rng)

param_dist = {
    "max_depth": [7, None],
    "max_features": randint(1, 30),
    "min_samples_split": randint(2, 30),
    "bootstrap": [True, False],
    "criterion": ["gini", "entropy"],
}

rsh = HalvingRandomSearchCV(
    estimator=clf, param_distributions=param_dist, factor=2, random_state=rng
)

rsh.fit(X_log_train, y_log_train)
rsh.best_params_

# # rsh.decision_function(X_test)

# y_log_pred = rsh.predict(X_log_test)


# fig, ax = plt.subplots(figsize=(8, 5))
# cmp = ConfusionMatrixDisplay(
#     confusion_matrix(y_log_test, y_log_pred),
#     display_labels=['C1', 'C2', 'C3', 'C4', 'C5', 'C6', 'CNA'],
# )

# cmp.plot(ax=ax)

# plt.savefig('figures/log_random_forest_best_params.png')
# # plt.show();



# #### OVR classification ####


# from sklearn.multiclass import OneVsRestClassifier

# # Init/fit
# ovr = OneVsRestClassifier(estimator=Perceptron())
# _ = ovr.fit(X_log_train, y_log_train)
# print(len(ovr.estimators_))


#### Preparing the PCAWG dataset
# cd PCAWG_RNA_seq; awk 'BEGIN{OFS="\t"}FNR == NR {a[$1]++; next}{if (FNR==1) {print; next};  split($1, gn, "."); if (a[gn[1]]) {$1=gn[1]; print}}' ../common_tcga_pcawg_genes.dat tophat_star_fpkm_uq.v2_aliquot_gl.tsv > tophat_star_fpkm_uq.v2_aliquot_gl_filtered.tsv

pcawg_exps = pd.read_csv("PCAWG_RNA_seq/tophat_star_fpkm_uq.v2_aliquot_gl_filtered.tsv", sep = "\t", index_col = 0, header = 0)

pcawg_exps_np = pcawg_exps.to_numpy().transpose()


pcawg_log = np.log(pcawg_exps_np + 1)

with open('npy_files/pcawg_log.npy', 'wb') as f:
    np.save(f, pcawg_log)


comb_exps = np.concatenate((exps_log, pcawg_log))


### finding the best parameter for learning

out_dict = {}

for i in np.arange(20, 160, 20):

    embedder_svd_i = TruncatedSVD(n_components = i)
    
    mod_log_all = embedder_svd_i.fit_transform(comb_exps)

    tcga_set = mod_log_all[: exps_log.shape[0], :]
    
    X_log_train, X_log_test, y_log_train, y_log_test = train_test_split(
        tcga_set, y, test_size=0.3, random_state=1121218
    )

    rsh.fit(X_log_train, y_log_train)
    rsh.best_params_
    
    # rsh.decision_function(X_test)
    
    y_log_pred = rsh.predict(X_log_test)
    
    
    fig, ax = plt.subplots(figsize=(8, 5))
    cmp = ConfusionMatrixDisplay(
        confusion_matrix(y_log_test, y_log_pred),
        display_labels=['C1', 'C2', 'C3', 'C4', 'C5', 'C6', 'CNA'],
    )
    
    cmp.plot(ax=ax)

    out_dict[i] = confusion_matrix(y_log_test, y_log_pred)
    
    plt.savefig('figures/log_random_forest_best_params_concat_' + str(i) + '.png')
    # plt.show();




#######################

with open('npy_files/exps_log.npy', 'rb') as f:
    exps_log = np.load(f)

with open('npy_files/pcawg_log.npy', 'rb') as f:
    pcawg_log = np.load(f)

comb_exps = np.concatenate((exps_log, pcawg_log))

spp_table = pd.read_csv('imm_subtypes_sorted.tsv', sep = '\t', header = None)
# spp_table = spp_table.set_index('TCGA.Participant.Barcode')

#### Seting up X and y ######

# X = mod_exps
# X_log = mod_log_exps
y = spp_table[2].to_numpy()


#### Best with i = 60
i = 60

embedder_svd_i = TruncatedSVD(n_components = i)

mod_log_all = embedder_svd_i.fit_transform(comb_exps)

tcga_set = mod_log_all[: exps_log.shape[0], :]

X_log_train, X_log_test, y_log_train, y_log_test = train_test_split(
    tcga_set, y, test_size=0.3, random_state=1121218
)

rsh.fit(X_log_train, y_log_train)
rsh.best_params_


y_log_pred = rsh.predict(X_log_test)


fig, ax = plt.subplots(figsize=(8, 5))
cmp = ConfusionMatrixDisplay(
    confusion_matrix(y_log_test, y_log_pred),
    display_labels=['C1', 'C2', 'C3', 'C4', 'C5', 'C6', 'CNA'],
    )

cmp.plot(ax=ax)

out_dict[i] = confusion_matrix(y_log_test, y_log_pred)

plt.savefig('figures/log_random_forest_best_params_concat_complex' + str(i) + '.png')
# plt.show();



pcawg_set = mod_log_all[exps_log.shape[0] :, :]

pcawg_predicted_labels = rsh.predict(pcawg_set)

np.savetxt("npy_files/RF_" + str(i) + "_predicted_complex_labels.txt", pcawg_predicted_labels, fmt = "%s")

