# This program predicts Ki value of a particular ligand molecule to target protein coagulation factor F11 by querying an already developed model developed with the following descriptors which are Shannon entropies, fractional Shannon entropy, ligand BEI and MW 

# getting the necessary packages
from tensorflow.keras.models import load_model
from sklearn.preprocessing import MinMaxScaler

import numpy as np
import pandas as pd

import matplotlib.pyplot as plt

import re
import math

# K-mer tokenization (optional)
from SmilesPE.pretokenizer import kmer_tokenizer


def F11a_Ki_query( input_vector, scaling_factor ):

    # make predictions on the testing data
    print("[INFO] predicting Ki...")
    model = load_model("F11a_target_hybrid_model_Ki_atom_labels_wo_Milvexian_BMS")
    
   	#preds = model.predict([image_series])
    preds = model.predict(input_vector)
  
    scaling_factor = scaling_factor   
  
    Ki_pred = preds.flatten()*scaling_factor 
    print("predicted pCHEMBL/MW:", Ki_pred )
  
    return Ki_pred


### Estimation of Shannon entropy of SMILES as a feature

SMI_REGEX_PATTERN = r"""(\[[^\]]+]|Br?|Cl?|N|O|S|P|F|I|b|c|n|o|s|p|\(|\)|\.|=|#|-|\+|\\|\/|:|~|@|\?|>>?|\*|\$|\%[0-9]{2}|[0-9])"""
regex = re.compile(SMI_REGEX_PATTERN)

def shannon_entropy_smiles(mol_smiles):
    
    molecule = mol_smiles 
    tokens = regex.findall(molecule)
    
    
    # tokens = kmer_tokenizer(molecule, ngram=10)
    # print(tokens)
    
    ### Frequency of each token generated
    L = len(tokens)
    L_copy = L
    tokens_copy = tokens
    
    num_token = []
    
    
    for i in range(0,L_copy):
        
        token_search = tokens_copy[0]
        num_token_search = 0
        
        if len(tokens_copy) > 0:
            for j in range(0,L_copy):
                if token_search == tokens_copy[j]:
                    # print(token_search)
                    num_token_search += 1
            # print(tokens_copy)        
                    
            num_token.append(num_token_search)   
                
            while token_search in tokens_copy:
                    
                tokens_copy.remove(token_search)
                    
            L_copy = L_copy - num_token_search
            
            if L_copy == 0:
                break
        else:
            pass
        
    # print(num_token)
    
    ### Calculation of Shannon entropy
    total_tokens = sum(num_token)
    
    shannon = 0
    
    for k in range(0,len(num_token)):
        
        pi = num_token[k]/total_tokens
        
        # print(num_token[k])
        # print(math.log2(pi))
        
        shannon = shannon - pi * math.log2(pi)
    
    # shannon = math.exp(-shannon)    
        
    return shannon    


if __name__ == '__main__':
    

    # Importing necessary packages
    from rdkit import Chem

    from rdkit.Chem.Draw import MolToImage, MolDrawOptions
    opts = MolDrawOptions()
    
    
    ### Prediction of Ki of MILVEXIAN
    
    # Importing SMILES of the MILVEXIAN molecule
    smiles_string = 'C[C@@H]1CCC[C@H](n2cnc(-c3cc(Cl)ccc3-n3cc(Cl)nn3)cc2=O)c2cc(ccn2)-c2c(cnn2C(F)F)NC1=O'
    mol = Chem.MolFromSmiles(smiles_string)
    
    SMILES_Shannon = shannon_entropy_smiles(smiles_string)
    
    # img =  MolToImage(mol)
    # img.save("Ki_targets_query_img/MILVEXIAN_atom_labels.png")
    

    # Take the above molecule as input in the following array form: [MW, SMILES Shannon, ligands BEI]
    MW = 626.46
    ligands_BEI = 15.96
    input_vector = [MW, SMILES_Shannon, ligands_BEI ]
    
    # Changing the array to DataFrame & transpose the DataFrame to make it horizontal
    input_vector= pd.DataFrame(input_vector).T
    
    # loading the dataset
    df = pd.read_csv('Ki_dataset_all_filtered_features_atom_labels_wo_Milvexian_BMS_with_smiles.csv', encoding='cp1252') 
    df_smiles = df['smiles'].values
    
    
    # estimating the SMILES Shannon entropy of the queried molecule
    # smiles regex definition
    SMI_REGEX_PATTERN = r"""(\[[^\]]+]|Br?|Cl?|N|O|S|P|F|I|b|c|n|o|s|p|\(|\)|\.|=|#|-|\+|\\|\/|:|~|@|\?|>>?|\*|\$|\%[0-9]{2}|[0-9])"""
    regex = re.compile(SMI_REGEX_PATTERN)
    shannon_arr = []
    
    for p in range(0,len(df_smiles)):
    
        molecule = df_smiles[p]
        tokens = regex.findall(molecule)
        
        # Frequency of each token generated
        L = len(tokens)
        L_copy = L
        tokens_copy = tokens
        
        num_token = []
        
        
        for i in range(0,L_copy):
            
            token_search = tokens_copy[0]
            num_token_search = 0
            
            if len(tokens_copy) > 0:
                for j in range(0,L_copy):
                    if token_search == tokens_copy[j]:
                        # print(token_search)
                        num_token_search += 1
                # print(tokens_copy)        
                        
                num_token.append(num_token_search)   
                    
                while token_search in tokens_copy:
                        
                    tokens_copy.remove(token_search)
                        
                L_copy = L_copy - num_token_search
                
                if L_copy == 0:
                    break
            else:
                pass
            
        # print(num_token)
        
        # Calculation of Shannon entropy
        total_tokens = sum(num_token)
        
        shannon = 0
        
        for k in range(0,len(num_token)):
            
            pi = num_token[k]/total_tokens
            
            # print(num_token[k])
            # print(math.log2(pi))
            
            shannon = shannon - pi * math.log2(pi)
            
        print("shannon entropy: ", shannon)
        shannon_arr.append(shannon)
        
        # shannon_arr.append(math.exp(-shannon))
        
        # print(shannon)     
    df['shannon_smiles']= shannon_arr  

    
    # constructing a new df column containng only MW, shannon values, ligand BEI and targets/labels
    df_1 = pd.DataFrame( df['full_mwt'].values)
    df_2 = pd.DataFrame(df['shannon_smiles'].values)
    df_5 = pd.DataFrame(df['pChEMBL Value'].values)

    df_new = pd.concat([ df_1, df_2, df['Ligand Efficiency BEI'] ], axis = 1)
    
    ## np.concatenate with DataFrame constructor could append DataFrames with a particular DataFrames's header
    XData_Ki_mod = pd.DataFrame(np.concatenate([df_new.values, input_vector.values]), columns=df_new.columns.values)
    

    cs = MinMaxScaler()
    XData_Ki_mod_cs = cs.fit_transform(XData_Ki_mod)
    
    # The last row is the query_vector for MELVEXIAN
    query_vector = XData_Ki_mod_cs[-1,:]
    
    
    # Converting the query_vector back into DataFrame for feeding it into neural network
    input_data = pd.DataFrame(query_vector).T
    print("The query vector to input: ", input_data)
    

    target_scaling_factor = 0.023161293914592727  # maxPrice for pChEMBL/MW Value for Ki in the entire dataset
    # target_scaling_factor = 0.023
    
    # predicted pchembl/MW value for MILVEXIAN
    Ki_pred = F11a_Ki_query( input_data, target_scaling_factor )
    
    # Denormalized by multiplying with MW
    Ki_pred = Ki_pred * MW
    
    pCHEMBL_val =  Ki_pred 
    print("The Predicted pCHEMBL value is (for MILVEXIAN): ", pCHEMBL_val)
    actual_Ki_val = 10**(9 - pCHEMBL_val)
    print("The Predicted Ki in nM is: ", actual_Ki_val)
    
  
