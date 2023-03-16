# MLP & MLP+CNN models for Ki 

Harnessing Shannon entropy of molecular symbols in deep neural networks to enhance prediction accuracy
------------------------------------------------------------------------------------------------------
This repository holds the codes pertaining to Figs. 2a-b of the article 'Harnessing Shannon entropy-based descriptors in machine learning models to enhance the prediction accuracy of molecular properties'.

Description
-----------
Shannon entropy framework has been demonstrated as an efficient descriptor for regression-type machine learning problem using (i) MLP-based and (ii) MLP+CNN based-deep neural networks. In this specific case, we model inhibition constant (Ki) in the form of pCHEMBL/MW values of potential binding molecules to the target protein: coagulation factor F11a. The specific objectives of the codes are described in the Notes section below. The basic dataset has been provided in the repository in the form of .csv files.

Usage
-----
1. Download or make a clone of the repository
2. Make a new conda environment using the environment file 'mlp_dnn.yml'
3. Run the python files directly using a python IDE or from command line

Example: python MLP_only_pchembl_MW_shannon_partial_shannon_f11a_Ki.py

Notes
-----
The function files are KiNet_mlp.py and datasets_molpred_2D_1image_resnet.py under image_processing folder. Therefore, directly run the other python files apart from these ones.

The objectives and usage of the rest of the scripts are as follows: Please run the python scripts directly or using the command line 'python <script_name.py> from the terminal

(i) Image dataset download or data acquisition: Run the chembl_target_mol_featurizer_Ki_atom_labels.py file directly to build the image dataset which will be saved in the folder Ki_targets_atom_labels_wo_Milvexian_BMS. This script also extracts and saves a descriptor set from the CHEMBL website to Ki_dataset_all_filtered_features_atom_labels_wo_Milvexian_BMS.csv.

(ii) CNN_MLP_pchembl_MW_partial_shannon_f11a_Ki.py: This script models Ki values in the form of pCHEMBL/MW values of molecules for the target protein coagulation factor F11a with fractional Shannon entropy and MW as descriptors. The model is a hybrid one comprising of MLP and CNN based hybrid deep neural network. The model predicts Ki values of the molecules as per the test data set.

(iii) MLP_only_pchembl_MW_partial_shannon_f11a_Ki.py: This program predicts Ki values (on test dataset) in the form of pCHEMBL/MW values of molecules for the target protein coagulation factor F11a with fractional Shannon entropy and MW as descriptors.

(iv) MLP_only_pchembl_MW_shannon_f11a_Ki.py: This program predicts Ki values (on test dataset) in the form of pCHEMBL/MW values of molecules for the target protein coagulation factor F11a with Shannon entropy and MW as descriptors.

(v) MLP_only_pchembl_MW_shannon_partial_shannon_bond_freq_f11a_Ki.py: This program predicts Ki values (on test dataset) in the form of pCHEMBL/MW values of molecules for the target protein coagulation factor F11a with Shannon entropies (SMILES/SMARTS/InChiKey-based), fractional Shannon entropy, bond frequency and MW as descriptors.

(vi) MLP_only_pchembl_MW_shannon_partial_shannon_f11a_Ki.py: This program predicts Ki values (on test dataset) in the form of pCHEMBL/MW values of molecules for the target protein coagulation factor F11a with fractional Shannon entropy and MW as descriptors.

(vii) MLP_only_pchembl_MW_shannon_partial_shannon_ligand_BEI_f11a_Ki.py: This program predicts Ki values in the form of pCHEMBL/MW values of molecules for the target protein coagulation factor F11a with Shannon entropies, fractional Shannon entropy, ligand BEI and MW as descriptors on the test dataset.

(viii) Query_MLP_only_MILVEXIAN.py: This program predicts Ki value of a particular, unknown ligand molecule (here MILVEXIAN as example) to target protein coagulation factor F11 by querying an already developed model (through running script (vii) above) using the following descriptors- Shannon entropy (SMILES), ligand BEI and MW of the molecule.
