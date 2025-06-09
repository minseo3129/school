import streamlit as st

# 화합물 기본 정보
scopolamine_smiles = "CN1C2CCC1CC(C2)OC(=O)C3=CC=CC=C3O"
cocaine_smiles = "CN1C2CCC1CC(C2)OC(=O)C3=CC=CC=C3C(=O)OC"
tanimoto_sim = 0.76  # 실제 유사도 값 예시

# 이미지 URL (미리 생성된 구조 이미지)
scopolamine_img = "https://raw.githubusercontent.com/your-repo/scopolamine.png"
cocaine_img = "https://raw.githubusercontent.com/your-repo/cocaine.png"

# 앱 타이틀
st.set_page_config(page_title="Scopolamine vs Cocaine 유사성", layout="centered")
st.title("💊 Scopolamine vs Cocaine 구조 유사성 분석")

st.markdown("""
이 앱은 **스코폴라민(Scopolamine)** (한방 약재 낭탕근 유래)과 **코카인(Cocaine)** 간의  
구조적 유사성을 **Tanimoto Similarity** 값을 통해 분석합니다.
""")

# 구조 이미지 및 SMILES 출력
col1, col2 = st.columns(2)
with col1:
    st.subheader("Scopolamine")
    st.image(scopolamine_img, caption="Scopolamine 구조", use_column_width=True)
    st.code(scopolamine_smiles)

with col2:
    st.subheader("Cocaine")
    st.image(cocaine_img, caption="Cocaine 구조", use_column_width=True)
    st.code(cocaine_smiles)

# 유사성 출력
st.markdown("---")
st.subheader("🔗 Tanimoto 유사도 (예시값)")
st.metric(label="Tanimoto Similarity", value=f"{tanimoto_sim:.2f}")

st.markdown("""
- **0.0**: 완전 불일치  
- **1.0**: 완전 일치  
- 일반적으로 **0.5 이상**이면 구조적 유사성이 있다고 간주됩니다.
""")

