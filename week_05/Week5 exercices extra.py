import rdkit
from rdkit import Chem
from rdkit.Chem import Draw
from rdkit.Chem import AllChem
import pandas as pd

drugs = pd.read_csv("/Users/thomaschristiansson/Documents/GitHub/practical-programming-in-chemistry-exercises/week_05/chembl_drugs.csv", sep= ";")
df = pd.DataFrame(drugs)
app1 = df[df["Phase"]==4]

df = pd.DataFrame(app1)
app2 =df[pd.notna(df["Smiles"])]

df = pd.DataFrame(app2)
app3 = df[df["Smiles"]!='.']




app2["My_Mol"] = app2["Smiles"].apply(Chem.MolFromSmiles)
app2.head(2)
print(app2)