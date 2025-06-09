import streamlit as st
from rdkit import Chem
from rdkit.Chem import MACCSkeys, AllChem, Draw, DataStructs
import matplotlib.pyplot as plt

# Title
st.title("🌿 Scopolamine vs Cocaine: Structural Similarity Analysis")
st.markdown("**Tanimoto Similarity** based on MACCS and Morgan fingerprints.")

# SMILES
scopolamine_smiles = 'CN1C2CCC3C(C2C(=O)C4=C1C=CC(=C4)O)OC(C3)C5CC5'
cocaine_smiles = 'CN1C(=O)C2C(C1C(=O)OC)C3=CC=CC=C3C2'

# Create molecules
scopolamine = Chem.MolFromSmiles(scopolamine_smiles)
cocaine = Chem.MolFromSmiles(cocaine_smiles)

# Display molecules
col1, col2 = st.columns(2)
with col1:
    st.subheader("Scopolamine")
    st.image(Draw.MolToImage(scopolamine, size=(300, 300)))
with col2:
    st.subheader("Cocaine")
    st.image(Draw.MolToImage(cocaine, size=(300, 300)))

# Similarity calculation
def tanimoto_similarity(mol1, mol2, method='MACCS'):
    if method == 'MACCS':
        fp1 = MACCSkeys.GenMACCSKeys(mol1)
        fp2 = MACCSkeys.GenMACCSKeys(mol2)
    elif method == 'Morgan':
        fp1 = AllChem.GetMorganFingerprintAsBitVect(mol1, 2, nBits=2048)
        fp2 = AllChem.GetMorganFingerprintAsBitVect(mol2, 2, nBits=2048)
    return DataStructs.TanimotoSimilarity(fp1, fp2)

sim_maccs = tanimoto_similarity(scopolamine, cocaine, 'MACCS')
sim_morgan = tanimoto_similarity(scopolamine, cocaine, 'Morgan')

# Show results
st.header("📊 Similarity Scores")
st.write(f"**MACCS Keys Similarity**: {sim_maccs:.2f}")
st.write(f"**Morgan Fingerprints Similarity**: {sim_morgan:.2f}")

# Plot
st.header("📈 Similarity Visualization")
fig, ax = plt.subplots()
methods = ['MACCS', 'Morgan']
scores = [sim_maccs, sim_morgan]
ax.bar(methods, scores, color=['skyblue', 'lightgreen'])
ax.set_ylim(0, 1)
ax.set_ylabel("Tanimoto Similarity")
ax.set_title("Structural Similarity Between Scopolamine and Cocaine")
st.pyplot(fig)
