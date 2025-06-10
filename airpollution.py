
# streamlit_app.py

import streamlit as st
import pandas as pd
import networkx as nx
import matplotlib.pyplot as plt

st.title("🌐 오염물질 상관 그래프 탐색기")
st.write("**여러 오염물질 간 상관관계를 시각화하고 탐색해보세요!**")

# 1. CSV 파일 업로드
uploaded_file = st.file_uploader("CSV 파일을 업로드하세요", type="csv")

if uploaded_file is not None:
    df = pd.read_csv(uploaded_file)

    # 2. 상관 분석용 변수 추출
    pollutants = ['CO', 'CO2', 'NO2', 'SO2', 'O3', 'PM2.5', 'PM10']
    pollutants = [p for p in pollutants if p in df.columns]
    df_pollutants = df[pollutants].dropna()

    # 3. 상관계수 계산
    corr_matrix = df_pollutants.corr()

    # 4. 그래프 생성
    G = nx.Graph()
    for col in corr_matrix.columns:
        for row in corr_matrix.index:
            if col != row:
                weight = corr_matrix.loc[row, col]
                if abs(weight) >= 0.6:
                    G.add_edge(row, col, weight=round(weight, 2))

    # 5. 그래프 시각화 출력
    st.subheader("📈 상관 그래프 (|r| ≥ 0.6)")
    fig, ax = plt.subplots(figsize=(8, 6))
    pos = nx.spring_layout(G, seed=42)
    nx.draw(G, pos, with_labels=True, node_color='lightblue', node_size=3000, font_size=12, font_weight='bold')
    edge_labels = {(u, v): d['weight'] for u, v, d in G.edges(data=True)}
    nx.draw_networkx_edge_labels(G, pos, edge_labels=edge_labels, font_color='red')
    st.pyplot(fig)

    # 6. 사용자 입력: 기체 이름 → 관련 상관 오염물질 출력
    st.subheader("🔎 특정 오염물질과 관련된 기체 찾기")
    gas = st.text_input("기체 이름을 입력하세요 (예: NO2, PM2.5, CO 등)")

    if gas in corr_matrix.columns:
        st.write(f"**'{gas}'와 상관관계 높은 오염물질 목록 (|r| ≥ 0.6):**")
        result = corr_matrix[gas].drop(gas)
        result = result[abs(result) >= 0.6].sort_values(ascending=False)
        st.dataframe(result)
    elif gas != "":
        st.warning("해당 기체는 데이터에 존재하지 않습니다.")