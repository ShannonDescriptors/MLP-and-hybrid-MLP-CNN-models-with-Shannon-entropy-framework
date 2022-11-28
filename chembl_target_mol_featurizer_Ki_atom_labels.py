### This script will download/ retrieve all 2D images of molecules and their corresponding features (or properties) as available from CHEMBL database given an input file list


from chembl_webresource_client.utils import utils
from chembl_webresource_client.new_client import new_client
import requests
import urllib, json

import pandas as pd


## Example for one specific molecule only
# mol = 'CHEMBL4082616'
# url = ('https://www.ebi.ac.uk/chembl/api/data/molecule/{}.json'.format(mol))
# response = urllib.request.urlopen(url)
# data = json.loads(response.read())
# print (data)
# print("\n\n\n\n")
# print("The molecule_properties: ", data["molecule_properties"])


# Getting the list of molecules from the .csv file
IC50_data_list = pd.read_csv('F11a_targets_Ki_only_equal_values_edited_wo_Milvexian_BMS.csv', encoding='cp1252') 
mol_list = IC50_data_list['Molecule ChEMBL ID'].values

# Defining the dictionary elements
alogp = []
aromatic_rings = []
cx_logd = []
cx_logp = []
cx_most_apka = []
cx_most_bpka = []
full_mwt = []
hba = []
hba_lipinski = []
hbd = []
hbd_lipinski = []
heavy_atoms = []
molecular_species = []
mw_freebase = []
mw_monoisotopic = []
num_lipinski_ro5_violations = []
num_ro5_violations = []
psa = []
qed_weighted = []
ro3_pass = []
rtb = []

# loop over the entire molecule list to build the database
for i in range(len(mol_list)):
    
    mol = mol_list[i]
    url = ('https://www.ebi.ac.uk/chembl/api/data/molecule/{}.json'.format(mol))
    response = urllib.request.urlopen(url)
    data = json.loads(response.read())
    
    if data['molecule_properties']['alogp'] is None:
        alogp.append('')
    else:
        alogp.append(   (data['molecule_properties']['alogp']) )
        
        
    if data['molecule_properties']['aromatic_rings'] is None:
        
        aromatic_rings.append('')     
    else:      
        aromatic_rings.append(   float(data['molecule_properties']['aromatic_rings'])  )
    
    
    if data['molecule_properties']['cx_logd'] is None:
        
          cx_logd.append('')     
    else:      
          cx_logd.append(   float(data['molecule_properties']['cx_logd'])  )
        
     
    if data['molecule_properties']['cx_logp'] is None:
        
          cx_logp.append('')     
    else:      
          cx_logp.append(   float(data['molecule_properties']['cx_logp'])  )    
    

    if data['molecule_properties']['cx_most_apka'] is None: 
        cx_most_apka.append(   ''  )
    else:
        cx_most_apka.append( float(data['molecule_properties']['cx_most_apka']) )
    

    if data['molecule_properties']['cx_most_bpka'] is None: 
        cx_most_bpka.append(   ''  )
    else:
        cx_most_bpka.append( float(data['molecule_properties']['cx_most_bpka']) )
        
    
    if data['molecule_properties']['full_mwt'] is None: 
        full_mwt.append(   ''  )
    else:
        full_mwt.append( float(data['molecule_properties']['full_mwt']) )
        
    if data['molecule_properties']['hba'] is None: 
        hba.append(   ''  )
    else:
        hba.append( float(data['molecule_properties']['hba']) )   
        
    if data['molecule_properties']['hba_lipinski'] is None: 
        hba_lipinski.append(   ''  )
    else:
        hba_lipinski.append( float(data['molecule_properties']['hba_lipinski']) )       
        
    if data['molecule_properties']['hbd'] is None: 
        hbd.append(   ''  )
    else:
        hbd.append( float(data['molecule_properties']['hbd']) )       
        
    if data['molecule_properties']['hbd_lipinski'] is None: 
        hbd_lipinski.append(   ''  )
    else:
        hbd_lipinski.append( float(data['molecule_properties']['hbd_lipinski']) )  
        
        
    if data['molecule_properties']['heavy_atoms'] is None: 
        heavy_atoms.append(   ''  )
    else:
        heavy_atoms.append( float(data['molecule_properties']['heavy_atoms']) )     
        
    
    if data['molecule_properties']['molecular_species'] is None:
        molecular_species.append(  '' )    
        
    else:    
        if data['molecule_properties']['molecular_species'] == 'NEUTRAL':
            molecular_species.append(  0.0  ) 
        
        if data['molecule_properties']['molecular_species'] == 'BASE':
            molecular_species.append(  1.0  )   
        
        if data['molecule_properties']['molecular_species'] == 'ACID':
            molecular_species.append(  -1.0  )
            
        if (data['molecule_properties']['molecular_species'] is not None) and ( data['molecule_properties']['molecular_species'] != 'NEUTRAL') and (data['molecule_properties']['molecular_species'] != 'BASE') and (data['molecule_properties']['molecular_species'] != 'ACID'):
            molecular_species.append(  0.5 )    
            
        
    if data['molecule_properties']['mw_freebase'] is None: 
        mw_freebase.append(   ''  )
    else:
        mw_freebase.append( float(data['molecule_properties']['mw_freebase']) ) 
        
    if data['molecule_properties']['mw_monoisotopic'] is None: 
        mw_monoisotopic.append(   ''  )
    else:
        mw_monoisotopic.append( float(data['molecule_properties']['mw_monoisotopic']) )     
        
    
    if data['molecule_properties']['num_lipinski_ro5_violations'] is None: 
        num_lipinski_ro5_violations.append(   ''  )
    else:
        num_lipinski_ro5_violations.append( float(data['molecule_properties']['num_lipinski_ro5_violations']) )  
    
    if data['molecule_properties']['num_ro5_violations'] is None: 
        num_ro5_violations.append(   ''  )
    else:
        num_ro5_violations.append( float(data['molecule_properties']['num_ro5_violations']) )  
        
    if data['molecule_properties']['psa'] is None: 
        psa.append(   ''  )
    else:
        psa.append( float(data['molecule_properties']['psa']) )   
        
    if data['molecule_properties']['qed_weighted'] is None: 
        qed_weighted.append(   ''  )
    else:
        qed_weighted.append( float(data['molecule_properties']['qed_weighted']) )     
        
                                                                    
    if data['molecule_properties']['ro3_pass'] == 'N':
        ro3_pass.append(  0.0  )
        
    else:
        ro3_pass.append(  1.0  )
            
    if data['molecule_properties']['rtb'] is None: 
        rtb.append(   ''  )
    else:
        rtb.append( float(data['molecule_properties']['rtb']) )  
        
       
    
features = {'alogp': alogp, 'aromatic_rings': aromatic_rings, 'cx_logd': cx_logd, 'cx_logp': cx_logp, 'cx_most_apka': cx_most_apka,
            'cx_most_bpka': cx_most_bpka, 'full_mwt': full_mwt, 'hba': hba, 'hba_lipinski': hba_lipinski, 'hbd': hbd, 'hbd_lipinski': hbd_lipinski, 'heavy_atoms': heavy_atoms,
            'molecular_species': molecular_species, 'mw_freebase': mw_freebase, 'mw_monoisotopic': mw_monoisotopic, 'num_lipinski_ro5_violations': num_lipinski_ro5_violations,
            'num_ro5_violations': num_ro5_violations, 'psa':  psa, 'qed_weighted': qed_weighted, 'ro3_pass': ro3_pass, 'rtb': rtb  }    
    
df_mol = pd.DataFrame(features)
print(df_mol)

# Concatenating more data (selected columns) from IC50_data_list and adding to df_mol
df_mol = pd.concat([  df_mol, IC50_data_list['Molecular Weight'], IC50_data_list['Ligand Efficiency BEI'], IC50_data_list['Ligand Efficiency LE'], IC50_data_list['Ligand Efficiency LLE'], IC50_data_list['Ligand Efficiency SEI'], IC50_data_list['pChEMBL Value'],  IC50_data_list['Smiles']   ], axis = 1)

# # Data Engineering
# # deleting all the rows corresponding to ''/'nan'/ 'None' after noting down the row index 

# # Finding the row indices where any column element is missing

# row_index_for_missing_col = []

# for i in range(len(mol_list)):
#     for j in range(df_mol.shape[1]):
        
#         if df_mol.iat[i,j] == '' or df_mol.iat[i,j] == 'None' or df_mol.iat[i,j] == 'nan': ## df_mol.iat[i,j] identifies a cell of the dataframe
            
#             if i not in row_index_for_missing_col:
            
#                 row_index_for_missing_col.append(i)
            

        
# # deleting the column with most of the missing data from df_mol and using the rest of the datatable 
# # dropping the column 'cx_most_bpka' as many values are missing
df_mol = df_mol.drop(['cx_most_bpka'], axis = 1)         
    
row_index_for_missing_col = []

for i in range(len(mol_list)):
    for j in range(df_mol.shape[1]):
        
        if df_mol.iat[i,j] == '' or df_mol.iat[i,j] == 'None' or pd.isnull(df_mol.iat[i,j]): ## df_mol.iat[i,j] identifies a cell of the dataframe
            
            if i not in row_index_for_missing_col:
            
                row_index_for_missing_col.append(i)

# New df_mol after dropping rows: row_index_for_missing_col
new_df_mol = df_mol.drop(row_index_for_missing_col)

# Collecting the indices of the final dataframe, new_df_mol
index = new_df_mol.index       

# Saving the index corresponding molecular images

# importing the necessary packages
from rdkit import Chem
from rdkit.Chem.Draw import MolToImage, MolDrawOptions
opts = MolDrawOptions()

# ## When the feature extraction is not in use
# new_df_mol = IC50_data_list

## Else: When feature extraction is in use
new_df_mol = new_df_mol
# Extracting the SMILES corresponding to the index values
df_SMILES = new_df_mol['Smiles'].values
print(df_SMILES)

# loop over to get the SMILES and corresponding images saved as .png
n = len(df_SMILES)
for i in range(n):
    
    mol = Chem.MolFromSmiles(df_SMILES[i])
    #print(df_SMILES[i])
    
    # for atom in mol.GetAtoms():
    #     j = atom.GetIdx()
    #     opts.atomLabels[j] = '' # '' getting rid of atom labels

    #     img = MolToImage(mol, options=opts)
    #     # display(img) 
        
    #     #img.save('Ki_dataset_all\Ki_{}.png'.format(df_target[i]))
    #     # img.save('IC50_targets_edited\m_{}.png'.format(i))
    
    img =  MolToImage(mol)
    # img.save('IC50_targets_atom_labels\m_{}.png'.format(i))
    img.save('Ki_targets_atom_labels_wo_Milvexian_BMS\m_{}.png'.format(i))
    print("Processed {}-th image".format(i))

# ## Final dataframe (new_df_mol) after dropping the SMILES column        
new_df_mol.drop(['Smiles'], axis = 1, inplace = True)

# ## Final dataframe (new_df_mol) after dropping the ro3_pass column        
new_df_mol.drop(['ro3_pass'], axis = 1, inplace = True)

# ## Saving the final dataframe as a csv
new_df_mol.to_csv('Ki_dataset_all_filtered_features_atom_labels_wo_Milvexian_BMS.csv', index=False)
