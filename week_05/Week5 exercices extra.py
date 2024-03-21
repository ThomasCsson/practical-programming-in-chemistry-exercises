import rdkit
from rdkit import Chem
from rdkit.Chem import Draw
from rdkit.Chem import AllChem
import pandas as pd
from rdkit.Chem import Descriptors
from mordred import Calculator, AtomCount



drugs = pd.read_csv("/Users/thomaschristiansson/Documents/GitHub/practical-programming-in-chemistry-exercises/week_05/chembl_drugs.csv", sep= ";")
df = pd.DataFrame(drugs)
app1 = df[df["Phase"]==4]
 
df = pd.DataFrame(app1)
app2 =df[pd.notna(df["Smiles"])]

df = pd.DataFrame(app2)
app3 = df[df["Smiles"]!='.']




app2["My_Mol"] = app2["Smiles"].apply(Chem.MolFromSmiles)
app2["n_atoms"] = app2["My_Mol"].apply(lambda x: x.GetNumAtoms())
app2["n_atoms_from_smiles"] = app2["Smiles"].apply(lambda x: Chem.MolFromSmiles(x).GetNumAtoms())
app2['NumRotBonds']= app2["My_Mol"].apply(lambda x: Descriptors.NumRotatableBonds(x))
app2['NumHDonors']= app2["My_Mol"].apply(lambda x: Descriptors.NumHDonors(x))
app2['NumHAcceptors']= app2["My_Mol"].apply(lambda x: Descriptors.NumHAcceptors(x))
app2['NumHDonors']= app2["My_Mol"].apply(lambda x: Descriptors.NumHDonors(x))
app2['MolWt']= app2["My_Mol"].apply(lambda x: Descriptors.MolWt(x))
app2['MolLogP']= app2["My_Mol"].apply(lambda x: Descriptors.MolLogP(x))
app2['NumValenceElectrons']= app2["My_Mol"].apply(lambda x: Descriptors.NumValenceElectrons(x))
app2['NumF']= app2["Smiles"].count('F')


app2.head(2)
print(app2)