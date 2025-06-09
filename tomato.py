# streamlit_app.py

import streamlit as st
import numpy as np
from PIL import Image
import os

# ===== 트리 노드 클래스 정의 =====
class TreeNode:
    def __init__(self, feature_name=None, threshold=None, label=None):
        self.feature_name = feature_name
        self.threshold = threshold
        self.label = label
        self.left = None
        self.right = None

    def is_leaf(self):
        return self.label is not None

# ===== 트리 샘플 구조 정의 (간단한 임계값 분류) =====
def build_tree():
    root = TreeNode("Red", 110)
    root.left = TreeNode("Green", 85)
    root.right = TreeNode("Blue", 130)

    root.left.left = TreeNode(label="Tomato_Early_blight")
    root.left.right = TreeNode(label="Tomato_Leaf_Mold")
    root.right.left = TreeNode(label="Tomato_Bacterial_spot")
    root.right.right = TreeNode(label="Tomato_healthy")

    return root

# ===== 이미지 RGB 평균값 추출 함수 =====
def extract_rgb(image):
    img = image.resize((128, 128))
    arr = np.array(img)
    if arr.ndim == 3:
        r, g, b = arr[:, :, 0], arr[:, :, 1], arr[:, :, 2]
        return np.mean(r), np.mean(g), np.mean(b)
    else:
        return 0, 0, 0

# ===== 트리 순회 및 경로 추적 =====
def traverse_tree(node, features, path, order='pre'):
    if not node:
        return

    if order == 'pre':
        path.append(node)

    if not node.is_leaf():
        value = features[node.feature_name]
        if value < node.threshold:
            traverse_tree(node.left, features, path, order)
        else:
            traverse_tree(node.right, features, path, order)

    if order == 'in':
        path.append(node)

    if order == 'post':
        path.append(node)

# ===== Streamlit 앱 실행 =====
st.set_page_config(layout="wide")
st.title("🌿 토마토 질병 분류 시각화 - 이진 탐색 트리 기반")

uploaded = st.file_uploader("이미지를 업로드하세요 (jpg, png)", type=["jpg", "jpeg", "png"])
order = st.radio("트리 순회 방식 선택", ["pre", "in", "post"], horizontal=True)

if uploaded:
    img = Image.open(uploaded)
    st.image(img, caption="업로드된 이미지", use_column_width=True)

    r, g, b = extract_rgb(img)
    st.markdown(f"**RGB 평균값**: 🔴 `{r:.2f}`, 🟢 `{g:.2f}`, 🔵 `{b:.2f}`")

    features = {"Red": r, "Green": g, "Blue": b}

    # 트리 생성 및 순회
    tree = build_tree()
    path = []
    traverse_tree(tree, features, path, order=order)

    st.subheader("🧭 분류 경로 (트리 순회 결과)")
    for i, node in enumerate(path):
        if node.is_leaf():
            st.success(f"[{i+1}] ✅ 최종 예측: **{node.label}**")
        else:
            val = features[node.feature_name]
            st.info(f"[{i+1}] `{node.feature_name} < {node.threshold}` → 현재값: `{val:.2f}` → 이동 방향: {'왼쪽' if val < node.threshold else '오른쪽'}")
else:
    st.warning("좌측에 토마토 잎 이미지를 업로드해 주세요.")
