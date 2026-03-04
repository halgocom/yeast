smiles ="CCCCCC1=CC(=C2[C@@H]3C=C(CC[C@H]3C(OC2=C1)(C)C)C)O"
from rdkit import Chem
from rdkit.Chem import EnumerateStereoisomers
molecule = Chem.MolFromSmiles(smiles)
options = EnumerateStereoisomers.StereoEnumerationOptions(unique=True, tryEmbedding=True)
isomers = tuple(EnumerateStereoisomers.EnumerateStereoisomers(
      molecule,
      options=options)
)
for smiles in isomers:
  print(Chem.MolToSmiles(smiles, isomericSmiles=True))
