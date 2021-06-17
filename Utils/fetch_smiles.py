# Fetches a file with SMILES string from Ligand Expo data
# then it creates a picke dictionary
import requests
import pandas as pd
import pickle

# Download Ligand expo data and convert it to a dataframe
file_name = "Components-smiles-stereo-oe.smi"
try:
    df = pd.read_csv(file_name, sep="\t",
                     header=None, 
                     names=["SMILES", "ID", "Name"])
except FileNotFoundError:
    url = f"http://ligand-expo.rcsb.org/dictionaries/{file_name}"
    print(url)
    r = requests.get(url, allow_redirects=True)
    open('Components-smiles-stereo-oe.smi', 'wb').write(r.content)
    df = pd.read_csv(file_name, sep="\t",
                     header=None, 
                     names=["SMILES", "ID", "Name"])
df.dropna(inplace=True)
df.drop("Name", axis=1, inplace=True)
df.set_index("ID", inplace=True)

# Tranfrom dataframe into a dictionary and store it as as pickle file 
smiles_dict = df.to_dict()

with open('smiles', 'ab') as f:
    pickle.dump(smiles_dict['SMILES'], f)

del smiles_dict

# Load the dictionary to see if it works
with open('smiles', 'rb') as f:
    smiles_dict = pickle.load(f)
    
print(smiles_dict['W11'])