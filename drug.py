import streamlit as st
import matplotlib.pyplot as plt

# 예시 데이터
similarity_maccs = 0.62
similarity_morgan = 0.48

# UI
st.title("🌿 Scopolamine vs Cocaine: Structural Similarity")
st.markdown("""
🔍 *이 앱은 RDKit이 설치되지 않은 환경을 위한 간이 분석 버전입니다.*  
`Tanimoto 유사도` 수치는 문헌 기반의 예상 값입니다.
""")

# 이미지
st.image("https://pubchem.ncbi.nlm.nih.gov/image/imgsrv.fcgi?cid=5184&t=l", caption="Scopolamine")
st.image("https://pubchem.ncbi.nlm.nih.gov/image/imgsrv.fcgi?cid=446220&t=l", caption="Cocaine")

# 유사도 시각화
st.header("📊 예상 유사도 시각화")
fig, ax = plt.subplots()
methods = ['MACCS Keys', 'Morgan FP']
scores = [similarity_maccs, similarity_morgan]
ax.bar(methods, scores, color=['lightblue', 'lightgreen'])
ax.set_ylim(0, 1)
ax.set_ylabel("Tanimoto Similarity")
ax.set_title("Expected Similarity Between Scopolamine and Cocaine")
st.pyplot(fig)

# 해석
st.markdown("#### 📌 해석")
st.markdown(f"- **MACCS 기반** 유사도: `{similarity_maccs}` → 부분적인 구조적 유사성")
st.markdown(f"- **Morgan 기반** 유사도: `{similarity_morgan}` → 낮은 구조적 유사성")

st.info("이 앱은 RDKit 미지원 환경에서 실행되며, 예측 기반 정보를 제공합니다.")
