1. Download or make a clone of the repository
2. Make a new conda environment using the environment file 'mlp_dnn.yml'
3. Run the python files directly using a python IDE or from command line

Notes:

The function file is KiNet_mlp.py. Therefore, directly run the other python files apart from this one.

The objectives and usage of the rest of the scripts are as follows: Please run the python scripts directly or using the command line 'python <script_name.py> from the terminal

(i) MLP_only_train_test_hybrid_pchembl_MW_without_shannon_entropy.py: This script could run MLP-based deep neural network model for pchembl/MW prediction (of IC50 of Tissue Factor Pathway Inhibitor) without using any Shannon entropies as descriptor and only using the MW as the sole descriptor. The target is MW-normalized pchembl value or pchembl/MW.

(ii) MLP_only_train_test_hybrid_pchembl_MW_with_shannon_entropy.py: This script could run MLP-based deep neural network model for pchembl/MW prediction (of IC50 of Tissue Factor Pathway Inhibitor) using Shannon entropy (SMILES/SMARTS/InChiKey-based) as descriptor along with the MW as the other descriptor

(iii) MLP_only_train_test_hybrid_pchembl_MW_with_shannon_entropy_smiles_smarts_inchikey_partial_shannon_smiles.py: This script could run MLP only model for pchembl/MW prediction (of IC50 of Tissue Factor Pathway Inhibitor) using combination of Shannon entropies (SMILES Shannon, SMARTS Shannon, InChiKey Shannon and SMILES partial/ fractional Shannon) and MW as descriptors. The SMILES-based fractional Shannon could add more accuracy to the model prediction.

(iv) MLP_only_train_test_hybrid_ligand_eff_without_shannon_entropy.py: This script could run MLP-based deep neural network model for predicting BEI (Ligand Binding Efficiency Index) values of molecules (ligands of Tissue Factor Pathway Inhibitor) without using any Shannon entropies as descriptor and only using the MW as the sole descriptor

(v) MLP_only_train_test_hybrid_ligand_eff_with_shannon_entropy.py:This script could run an MLP only model for predicting BEI values of molecules (ligands of Tissue Factor Pathway Inhibitor) with using Shannon entropy as a descriptor along with using the MW as the other descriptor

(vi) MLP_only_train_test_hybrid_ligand_eff_with_shannon_entropy_smiles_smarts_inchikey_with_partial_shannon.py: This script could run MLP-based deep neural network model for predicting BEI values of molecules (ligands of Tissue Factor Pathway Inhibitor) using combination of Shannon entropies (SMILES Shannon, SMARTS Shannon, InChiKey Shannon and SMILES partial/ fractional Shannon) and MW as descriptors

(vii) MLP_only_train_test_hybrid_ligand_eff_with_morgan_fingerprint.py:This script could run MLP-based deep neural network model for predicting BEI values of molecules (ligands of Tissue Factor Pathway Inhibitor) using different combinations of Morgan Fingerprint (with or without), Shannon entropies (SMILES Shannon and SMILES partial/ fractional Shannon) and MW as descriptors

(viii) MLP_only_train_test_hybrid_pCHEMBL_with_shannon_entropy_and_other descriptors.py: This script could run MLP-based deep neural network model for pchembl/MW prediction (of IC50 of Tissue Factor Pathway Inhibitor) using Shannon entropy as descriptor along with using BEI and MW as other descriptors (Ligand Efficiency BEI prediction from script#(vii))
