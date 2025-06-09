import streamlit as st
from rdkit import Chem
from rdkit.Chem import Draw
from rdkit.Chem.AllChem import GetMorganFingerprintAsBitVect
from rdkit.DataStructs import TanimotoSimilarity
import pandas as pd

# 내부 비교용 SMILES 데이터 (예시)
data_smiles = [
    ("Aspirin", "CC(=O)OC1=CC=CC=C1C(=O)O"),
    ("Caffeine", "CN1C=NC2=C1C(=O)N(C(=O)N2C)C"),
    ("Acetaminophen", "CC(=O)NC1=CC=C(O)C=C1"),
    ("Ibuprofen", "CC(C)CC1=CC=C(C=C1)C(C)C(=O)O"),
    ("Naproxen", "CC1=CC=C(C=C1)C(O)=C(C)C(=O)O")
]

st.title("화합물 구조 유사성 비교기")
user_smiles = st.text_input("SMILES 화학식을 입력하세요:")

if user_smiles:
    try:
        user_mol = Chem.MolFromSmiles(user_smiles)
        user_fp = GetMorganFingerprintAsBitVect(user_mol, radius=2, nBits=2048)

        results = []
        for name, smiles in data_smiles:
            mol = Chem.MolFromSmiles(smiles)
            fp = GetMorganFingerprintAsBitVect(mol, radius=2, nBits=2048)
            sim = TanimotoSimilarity(user_fp, fp)
            results.append((name, smiles, sim, mol))

        results.sort(key=lambda x: x[2], reverse=True)

        st.subheader("유사한 화합물:")
        for name, smiles, sim, mol in results[:5]:
            st.write(f"**{name}** - Tanimoto 유사도: {sim:.3f}")
            st.image(Draw.MolToImage(mol, size=(200,200)))

    except Exception as e:
        st.error(f"입력 오류: {e}")
