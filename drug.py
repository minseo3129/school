import streamlit as st
from rdkit import Chem
from rdkit.Chem import Draw
from rdkit.Chem import MACCSkeys, DataStructs

# -------------------------------
# 화합물 SMILES 정의
# -------------------------------
scopolamine_smiles = "CN1C2CCC1CC(C2)OC(=O)C3=CC=CC=C3O"
cocaine_smiles = "CN1C2CCC1CC(C2)OC(=O)C3=CC=CC=C3C(=O)OC"

# RDKit 분자 객체 생성
mol_scopolamine = Chem.MolFromSmiles(scopolamine_smiles)
mol_cocaine = Chem.MolFromSmiles(cocaine_smiles)

# MACCS Keys Fingerprint 생성
fp_scopolamine = MACCSkeys.GenMACCSKeys(mol_scopolamine)
fp_cocaine = MACCSkeys.GenMACCSKeys(mol_cocaine)

# Tanimoto 유사도 계산
tanimoto_sim = DataStructs.FingerprintSimilarity(fp_scopolamine, fp_cocaine)

# -------------------------------
# Streamlit UI 구성
# -------------------------------
st.set_page_config(page_title="Scopolamine vs Cocaine 유사성 분석", layout="centered")
st.title("💊 Scopolamine vs Cocaine 구조 유사성 분석")

st.markdown("""
이 앱은 **스코폴라민(Scopolamine)** (한방 약재 낭탕근 유래)과 **코카인(Cocaine)** 간의  
구조적 유사성을 **MACCS Keys** 및 **Tanimoto Similarity**를 기반으로 분석합니다.
""")

# 구조 이미지 출력
col1, col2 = st.columns(2)
with col1:
    st.subheader("Scopolamine")
    st.image(Draw.MolToImage(mol_scopolamine, size=(300, 300)))
    st.code(scopolamine_smiles)

with col2:
    st.subheader("Cocaine")
    st.image(Draw.MolToImage(mol_cocaine, size=(300, 300)))
    st.code(cocaine_smiles)

# 유사성 출력
st.markdown("---")
st.subheader("🔗 Tanimoto 유사도 (MACCS Keys 기반)")
st.metric(label="Tanimoto Similarity", value=f"{tanimoto_sim:.3f}")

st.markdown("""
- **0.0**: 완전 불일치  
- **1.0**: 완전 일치  
- 일반적으로 **0.5 이상**이면 구조적 유사성이 있다고 간주됩니다.
""")
