import streamlit as st
from rdkit import Chem
from rdkit.Chem import MACCSkeys, AllChem, Draw, DataStructs
import matplotlib.pyplot as plt

# Title
st.title("🌿 한방 약재 vs 마약성 물질: 구조적 유사성 분석")
st.markdown("""
이 앱은 **스코폴라민(Scopolamine)**과 **코카인(Cocaine)** 간의 구조적 유사성을 분석하고 시각화합니다.  
사용된 기법은 **MACCS Keys**와 **Morgan Fingerprints** 기반의 **Tanimoto Similarity**입니다.
""")

# SMILES
scopolamine_smiles = 'CN1C2CCC3C(C2C(=O)C4=C1C=CC(=C4)O)OC(C3)C5CC5'
cocaine_smiles = 'CN1C(=O)C2C(C1C(=O)OC)C3=CC=CC=C3C2'

# Molecules
scopolamine = Chem.MolFromSmiles(scopolamine_smiles)
cocaine = Chem.MolFromSmiles(cocaine_smiles)

# Draw Molecules
st.header("🧪 화합물 구조 시각화")
col1, col2 = st.columns(2)
with col1:
    st.subheader("Scopolamine")
    st.image(Draw.MolToImage(scopolamine, size=(300, 300)))
with col2:
    st.subheader("Cocaine")
    st.image(Draw.MolToImage(cocaine, size=(300, 300)))

# Fingerprint-based similarity
def calculate_similarity(mol1, mol2, method):
    if method == 'MACCS':
        fp1 = MACCSkeys.GenMACCSKeys(mol1)
        fp2 = MACCSkeys.GenMACCSKeys(mol2)
    elif method == 'Morgan':
        fp1 = AllChem.GetMorganFingerprintAsBitVect(mol1, radius=2, nBits=2048)
        fp2 = AllChem.GetMorganFingerprintAsBitVect(mol2, radius=2, nBits=2048)
    return DataStructs.TanimotoSimilarity(fp1, fp2)

sim_maccs = calculate_similarity(scopolamine, cocaine, 'MACCS')
sim_morgan = calculate_similarity(scopolamine, cocaine, 'Morgan')

# Display results
st.header("📊 Tanimoto 유사도 계산")
st.markdown(f"- **MACCS Keys 기반 유사도**: `{sim_maccs:.2f}`")
st.markdown(f"- **Morgan Fingerprints 기반 유사도**: `{sim_morgan:.2f}`")

# Visualization
st.header("📈 유사도 시각화")
fig, ax = plt.subplots()
methods = ['MACCS', 'Morgan']
scores = [sim_maccs, sim_morgan]
ax.bar(methods, scores)
ax.set_ylim(0, 1)
ax.set_ylabel("Tanimoto Similarity")
ax.set_title("Scopolamine vs Cocaine 구조 유사도")
st.pyplot(fig)

# Interpretation
def interpret_score(score):
    if score > 0.7:
        return "매우 높은 구조적 유사성"
    elif score > 0.5:
        return "부분적인 구조적 유사성"
    else:
        return "낮은 구조적 유사성"

st.header("📌 해석")
st.markdown(f"**MACCS 해석**: {interpret_score(sim_maccs)}")
st.markdown(f"**Morgan 해석**: {interpret_score(sim_morgan)}")

