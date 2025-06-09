import streamlit as st
from rdkit import Chem
from rdkit.Chem import MACCSkeys, AllChem, Draw, DataStructs
from rdkit.Chem.Draw import rdMolDraw2D
import matplotlib.pyplot as plt
import numpy as np

# App Title
st.title("한방 약재와 마약성 물질 간 구조적 유사성 분석")
st.markdown("""
본 애플리케이션은 **스코폴라민(Scopolamine)**과 **코카인(Cocaine)** 간의 구조적 유사성을 MACCS Keys 및 Morgan Fingerprints 기반으로 분석하고, **Tanimoto 유사도**를 계산합니다.
""")

# SMILES 입력
scopolamine_smiles = 'CN1C2CCC3C(C2C(=O)C4=C1C=CC(=C4)O)OC(C3)C5CC5'
cocaine_smiles = 'CN1C(=O)C2C(C1C(=O)OC)C3=CC=CC=C3C2'

# 분자 생성
scopolamine = Chem.MolFromSmiles(scopolamine_smiles)
cocaine = Chem.MolFromSmiles(cocaine_smiles)

st.header("1️⃣ 분자 구조")
col1, col2 = st.columns(2)

with col1:
    st.subheader("스코폴라민")
    st.image(Draw.MolToImage(scopolamine), use_column_width=True)

with col2:
    st.subheader("코카인")
    st.image(Draw.MolToImage(cocaine), use_column_width=True)

# 유사성 계산 함수
def calculate_similarity(mol1, mol2, method='MACCS'):
    if method == 'MACCS':
        fp1 = MACCSkeys.GenMACCSKeys(mol1)
        fp2 = MACCSkeys.GenMACCSKeys(mol2)
    elif method == 'Morgan':
        fp1 = AllChem.GetMorganFingerprintAsBitVect(mol1, 2, nBits=2048)
        fp2 = AllChem.GetMorganFingerprintAsBitVect(mol2, 2, nBits=2048)
    return DataStructs.TanimotoSimilarity(fp1, fp2)

# 유사성 계산
similarity_maccs = calculate_similarity(scopolamine, cocaine, method='MACCS')
similarity_morgan = calculate_similarity(scopolamine, cocaine, method='Morgan')

# 결과 출력
st.header("2️⃣ 구조적 유사성 (Tanimoto Similarity)")
st.markdown(f"- **MACCS Keys 유사도**: `{similarity_maccs:.2f}`")
st.markdown(f"- **Morgan Fingerprints 유사도**: `{similarity_morgan:.2f}`")

# 해석
st.header("3️⃣ 유사성 해석")
def interpret_similarity(score):
    if score > 0.7:
        return "매우 높은 구조적 유사성"
    elif score > 0.5:
        return "부분적인 구조적 유사성"
    else:
        return "구조적으로 다름"

st.write(f"**MACCS 해석:** {interpret_similarity(similarity_maccs)}")
st.write(f"**Morgan 해석:** {interpret_similarity(similarity_morgan)}")

# 부가 설명
st.markdown("""
---
### 🔬 추가 분석 제안
- Morphine과의 유사성 분석
- 약리 작용 비교: 콜린성/도파민성 영향
- 머신러닝 기반 예측 모델 구축
- PubChem 및 CNS 데이터 기반 약물 반응 예측

📚 참고 자료:
- RDKit: https://www.rdkit.org/
- PubChem: https://pubchem.ncbi.nlm.nih.gov/
""")
