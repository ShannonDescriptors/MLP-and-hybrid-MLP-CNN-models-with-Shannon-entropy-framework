# This program predicts Ki values in the form of pCHEMBL/MW values of molecules for the target protein coagulation factor F11a with fractional Shannon entropy and MW as descriptors


# import the necessary packages
from imutils import paths
import random
import os
import numpy as np
import pandas as pd


from KiNet_mlp import KiNet_mlp

from tensorflow.keras.optimizers import Adam, Adagrad
from tensorflow.keras.optimizers import SGD
from sklearn.model_selection import train_test_split
from tensorflow.keras.layers import Dense
from tensorflow.keras.models import Model
from tensorflow.keras.layers import concatenate
from sklearn.preprocessing import MinMaxScaler
from tensorflow.keras.callbacks import EarlyStopping, TensorBoard
import matplotlib.pyplot as plt


from rdkit import Chem
import re
import math

import scipy
from sklearn.metrics import mean_absolute_percentage_error,mean_absolute_error,mean_squared_error


# Getting the list of molecules from the .csv
df = pd.read_csv('Ki_dataset_all_filtered_features_atom_labels_wo_Milvexian_BMS_with_smiles.csv', encoding='cp1252') 
df_MW = df['full_mwt'].values ## The MW for ach molecule
df_target = df['pChEMBL Value'].values ## Defining df_target
df_smiles = df['smiles'].values

### Calculating MW normalized df_target
df_target_norm = df['pChEMBL Value'].values ## Defining normalized df_target
for i in range(0,len(df)):
    df_target_norm[i] = df_target[i]/df_MW[i] 
df['pChEMBL Value'] = df_target_norm
print("Shape of the pCHEMBL labels array", df_target.shape)
    
# getting the max/min values of the target column
maxPrice = df.iloc[:,-1].max() 
minPrice = df.iloc[:,-1].min() 
print(maxPrice,minPrice)



### Smiles shannon as a feature
## Generate a new column with title 'shannon_smiles'. Evaluate the Shannon entropy for each smile string and store into 'shannon_smiles' column

### Inserting the new column as the 2nd column in df
df.insert(1, "shannon_smiles", 0.0)


###------------------------------------------------------Shannon Entropy Generation: SMILES--------------------------------------------------------------------------------------------------

### Generate a new column with title 'shannon_smiles'. Evaluate the Shannon entropy for each smile string and store into 'shannon_smiles' column

# smiles regex definition
SMI_REGEX_PATTERN = r"""(\[[^\]]+]|Br?|Cl?|N|O|S|P|F|I|b|c|n|o|s|p|\(|\)|\.|=|#|-|\+|\\|\/|:|~|@|\?|>>?|\*|\$|\%[0-9]{2}|[0-9])"""
regex = re.compile(SMI_REGEX_PATTERN)

shannon_arr = []
for p in range(0,len(df_smiles)):
    
    molecule = df_smiles[p]
    tokens = regex.findall(molecule)
    
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
    
    # import math
    shannon = 0
    
    for k in range(0,len(num_token)):
        
        pi = num_token[k]/total_tokens
        
        # print(num_token[k])
        # print(math.log2(pi))
        
        shannon = shannon - pi * math.log2(pi)
        
    print("shannon entropy: ", shannon)
    # shannon_arr.append(shannon)
    
    shannon_arr.append(math.exp(-shannon))
        
# print(shannon)     
df['shannon_smiles']= shannon_arr        


###---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

###------------------------------------------------------Shannon Entropy Generation: SMARTS--------------------------------------------------------------------------------------------------

### Generate a new column with title 'shannon_smarts'. Evaluate the Shannon entropy for each smile string and store into 'shannon_smarts' column

### smarts
#### Inserting the new column as the 2nd column in df
df.insert(2, "shannon_smarts", 0.0)

# smarts regex definition
smarts_REGEX_PATTERN =  r"""(\[[^\]]+]|Br?|Cl?|N|O|S|P|F|I|b|c|n|o|s|p|\(|\)|\.|=|#|-|\+|\\|\/|:|~|@|\?|>>?|\*|\$|\%[0-9]{2}|[0-9])"""
regex = re.compile(smarts_REGEX_PATTERN)


shannon_arr = []
for p in range(0,len(df['shannon_smarts'])):
    
    smiles = df_smiles[p]
    mol = Chem.MolFromSmiles(smiles)
    # mol = Chem.AddHs(mol)
    m_smarts = Chem.MolToSmarts(mol)
    
    molecule = m_smarts
    tokens = regex.findall(molecule)
    # print(tokens)
    # tokens = kmer_tokenizer(molecule, ngram=10)
    
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
        
    print("shannon entropy on smarts: ", shannon)
    shannon_arr.append(shannon)  
    
df['shannon_smarts']= shannon_arr  
###---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------



###------------------------------------------------------Shannon Entropy Generation: InChiKey---------------------------------------------------------------------------------------------------------------------------------------------

### Generate a new column with title 'shannon_smiles'. Evaluate the Shannon entropy for each smile string and store into 'shannon_smiles' column

### inchikey
### InChiKey: reading all letters except special characters
df.insert(3, "shannon_inchikey", 0.0)  

# InChiKey regex definition
InChiKey_REGEX_PATTERN = r"""([A-Z])"""
regex = re.compile(InChiKey_REGEX_PATTERN)

def shannon_entropy_inch(ik):
    
    molecule = ik
    tokens = regex.findall(molecule)
    # print(tokens)
    # tokens = kmer_tokenizer(molecule, ngram=10)
    
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
    
    # import math
    shannon = 0
    
    for k in range(0,len(num_token)):
        
        pi = num_token[k]/total_tokens
        
        # print(num_token[k])
        # print(math.log2(pi))
        
        shannon = shannon - pi * math.log2(pi)
    
    # shannon = math.exp(-shannon)    
        
    return shannon  


shannon_arr = []
for p in range(0,len(df['shannon_inchikey'])):

    smiles = df_smiles[p]
    mol = Chem.MolFromSmiles(smiles)
    # mol = Chem.AddHs(mol)
    m_inchkey = Chem.MolToInchiKey(mol)
    ik = m_inchkey.split("-")
    # print("InchiKey splitted", ik)
    # print("\n")

    shannon_entropy = 0
    for i in range(len(ik)):
    
        if i<=1:
            shannon_entropy_by_parts = shannon_entropy_inch(ik[i])
            shannon_entropy = shannon_entropy + shannon_entropy_by_parts  
            # print(shannon_entropy_by_parts)
        else:
            freq = 1/25 ### Inchikey contains total 25 characters, apart from 2 hyphens
            shannon_entropy_by_parts = - freq * math.log2(freq)
            shannon_entropy = shannon_entropy + shannon_entropy_by_parts  
            # print(shannon_entropy_by_parts)
        
        
    print("shannon entropy on inchikey: ", shannon_entropy)
    shannon_arr.append(shannon_entropy) 
    
df['shannon_inchikey']= shannon_arr    
###---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------


# generating a dictionary of atom occurrence frequencies
def freq_atom_list(atom_list_input_mol):
    
    atom_list = ['H', 'C', 'N', 'O', 'S', 'P','F', 'Cl', 'Br', 'I'] 
    dict_freq = {}
    
    ### adding keys
    for i in range(len(atom_list)):
        dict_freq[atom_list[i]] = 0  ### The values are all set 0 initially
    # print(dict_freq)
    
    ### update the value by 1 when a key in encountered in the string
    for i in range(len(atom_list_input_mol)):
        dict_freq[ atom_list_input_mol[i] ] = dict_freq[ atom_list_input_mol[i] ] + 1
    
    ### The dictionary values as frequency array
    freq_atom_list =  list(dict_freq.values())/ (  sum(  np.asarray (list(dict_freq.values()))  )    )
    
    # print(list(dict_freq.values()))
    # print(freq_atom_list )
    
    ### Getting the final frequency dictionary
    ### adding values to keys
    for i in range(len(atom_list)):
        dict_freq[atom_list[i]] = freq_atom_list[i]  
        
    # print(dict_freq)
    freq_atom_list = dict_freq
    
        
    return freq_atom_list


# getting the maximum length of SMILES strings in the database
len_smiles = []
for j in range(0,len(df_smiles)):
    
    mol = Chem.MolFromSmiles(df_smiles[j])
    
    ### No H atom is considered
    # mol = Chem.AddHs(mol)
    
    k=0
    for atom in mol.GetAtoms():
        k = k +1 
    len_smiles.append(k)

max_len_smiles = max(len_smiles)


# Constructing the padded array of partial or fractional Shannon entropy per molecule
def ps_padding(ps, max_len_smiles):
    
    len_ps = len(ps)
    
    len_forward_padding = int((max_len_smiles - len_ps)/2)
    len_back_padding = max_len_smiles - len_forward_padding - len_ps
    
    ps_padded = list(np.zeros(len_forward_padding))  + list(ps) + list(np.zeros(len_back_padding))
    
    return ps_padded 


# getting the total Shannon and partial or fractional shannon data for all molecules
fp_combined = []
shannon_arr = []

for i in range(0,len(df_smiles)):  
    
  mol = Chem.MolFromSmiles(df_smiles[i])
  
  # estimating the partial or fractional shannon for an atom type => the current node
  total_shannon = df['shannon_smiles'][i]
  shannon_arr.append( total_shannon )
  
  # The atom list as per rdkit in string form
  atom_list_input_mol = []
  for atom_rdkit in mol.GetAtoms():
     atom_list_input_mol.append(str(atom_rdkit.GetSymbol()))     
        
     
  freq_list_input_mol = freq_atom_list(atom_list_input_mol)
  
  # estimating the fractional/partial Shannon entropy
  ps = []
  for atom_rdkit in mol.GetAtoms():
      atom_symbol = atom_rdkit.GetSymbol()
      atom_type = atom_symbol # atom symbol in atom type
      
      partial_shannon = freq_list_input_mol[atom_type] * total_shannon
      ps.append( partial_shannon )
      # ps.append( freq_list_input_mol[atom_type] )

  ps_arr = ps_padding(ps, max_len_smiles)     
  fp_combined.append(ps_arr)

# Smiles shannon as a feature
shannon_smiles = pd.DataFrame(shannon_arr, columns = ["shannon_smiles"] )

# partial shannon_entropy as feature
fp_mol = pd.DataFrame(fp_combined)

# constructing a new df table as df_new containng only MW, shannon values and target labels. Shannon entropy values could include shannon_smiles, shannon_smarts, shannon_inchikey and partial shannon
df_1 = pd.DataFrame( df['full_mwt'].values)
df_2 = pd.DataFrame(df['shannon_smiles'].values)
df_3 = pd.DataFrame(df['shannon_smarts'].values)
df_4 = pd.DataFrame(df['shannon_inchikey'].values)
df_5 = pd.DataFrame(df['pChEMBL Value'].values)

df_new = pd.concat([ df_1, df_2, df_3, df_4, fp_mol, df_5 ], axis = 1)

print("[INFO] constructing training/ testing split")
split = train_test_split( df_new, test_size = 0.15, random_state = 10)

# Distribute the input data columns in train & test splits  
(XtrainTotalData, XtestTotalData) = split  # split format always is in (data_train, data_test, label_train, label_test)
 
# Normalizing the test or labels column
XtrainLabels = (XtrainTotalData.iloc[:,-1])/ (maxPrice)
XtestLabels = (XtestTotalData.iloc[:,-1])/ (maxPrice)  

XtrainData = (XtrainTotalData.iloc[:,0:-1])
XtestData = (XtestTotalData.iloc[:,0:-1])

print("XtrainData shape",XtrainData.shape)
print("XtestData shape",XtestData.shape)


# perform min-max scaling each continuous feature column to the range [0 1]
cs = MinMaxScaler()
trainContinuous = cs.fit_transform(XtrainTotalData.iloc[:,0:XtrainTotalData.shape[1]-1])
testContinuous = cs.transform(XtestTotalData.iloc[:,0:XtestTotalData.shape[1]-1])

print("[INFO] processing input data after normalization....")
XtrainData, XtestData = trainContinuous,testContinuous
print("XtrainData shape",XtrainData.shape)
print("XtestData shape",XtestData.shape)


# create the MLP model
mlp = KiNet_mlp.create_mlp(XtrainData.shape[1], regress = False) # the input dimension to mlp would be shape[1] of the matrix i.e. column features
print("shape of mlp", mlp.output.shape)


# Processing output of MLP model
combinedInput = mlp.output
print("shape of combinedInput",combinedInput.shape)

# Defining the final FC (Dense) layers 
x = Dense(100, activation = "relu") (combinedInput)
x = Dense(10, activation = "relu") (combinedInput)
x = Dense(1, activation = "linear") (x)
print("shape of x", x.shape)


# Defining the FINAL MODEL as 'model'
model = Model(inputs = mlp.input, outputs = x)
print("shape of mlp input", mlp.input.shape)


# initialize the optimizer model and compile the model
print("[INFO] compiling model...")
opt = SGD(lr= 1.15e-5, decay = 1.15e-5/200)
model.compile(loss = "mean_absolute_percentage_error", optimizer = opt)

#train the network
print("[INFO] training network...")
trainY = XtrainLabels
testY = XtestLabels


epoch_number = 200 ## Larger epoch number would increase the accuracy, e.g. epoch number = 2000
BS = 5;

# Defining the early stop to monitor the validation loss to avoid overfitting.
early_stop = EarlyStopping(monitor='val_loss', min_delta=0, patience=1000, verbose=1, mode='auto')

# shuffle = False to reduce randomness and increase reproducibility
H = model.fit( x = XtrainData, y = trainY, validation_data = (XtestData, testY), batch_size = BS, epochs = epoch_number, verbose=1, shuffle=False, callbacks = [early_stop]) 

# evaluate the network
print("[INFO] evaluating network...")
preds = model.predict(XtestData,batch_size=BS)

# compute the difference between the predicted and actual values, then compute the % difference and abs % difference
diff = preds.flatten() - testY
PercentDiff = (diff/testY)*100
absPercentDiff = (np.abs(PercentDiff)).values

# compute the mean and standard deviation of absolute percentage difference
mean = np.mean(absPercentDiff)
std = np.std(absPercentDiff)
print("[INF0] mean {: .2f},  std {: .2f}".format(mean, std) )

## Plotting Predicted vs Actual pCHEMBL/MW values
N = len(testY)
colors = np.random.rand(N)
x = testY * maxPrice
y =  (preds.flatten()) * maxPrice
plt.scatter(x, y, c=colors)
plt.plot( [0,maxPrice],[0,maxPrice] )
plt.xlabel('Actual Ki: pCHEMBL Value/MW', fontsize=10)
plt.ylabel('Predicted Ki: pCHEMBL Value/MW', fontsize=10)
plt.savefig('MLP_only_Ki_pred_atom_labels_wo_Milvexian_BMS_shannon_partial_shannon.png')
plt.show()


# ### MAE as a function
# def mae(y_true, predictions):
#     y_true, predictions = np.array(y_true), np.array(predictions)
#     return np.mean(np.abs(y_true - predictions))

# ### The MAE estimated
# print("The mean absolute error estimated: {}".format( mae(x, y) )) 

###-------------------------------------------------------------------------------------------------------Evaluating standard statistics of model performance--------------------------------------------------------------------
### The MAPE
print("The mean absolute percentage error: {}".format( mean_absolute_percentage_error(x, y) ) )   

### The MAE
print("The mean absolute error: {}".format( mean_absolute_error(x, y) ))    
    
### The MSE
print("The mean squared error: {}".format( mean_squared_error(x, y) ))  

### The RMSE
print("The root mean squared error: {}".format( mean_squared_error(x, y, squared= False) ) )    

### General stats
slope, intercept, r_value, p_value, std_err = scipy.stats.linregress(x, y)
print("The R^2 value between actual and predicted target:", r_value**2)

# plot the training loss and validation loss
plt.style.use("ggplot")
plt.figure()
plt.plot(np.arange(0, epoch_number ), H.history["loss"], label="train_loss")
plt.plot(np.arange(0, epoch_number ), H.history["val_loss"], label="val_loss")
plt.title("Training loss and Validation loss")
plt.xlabel("Epoch #")
plt.ylabel("Loss")
plt.legend()
plt.savefig('loss@epochs_pChEMBL_MW_MLP_only_shannon_partial_shannon')
####------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------