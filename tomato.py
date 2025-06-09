import streamlit as st
import numpy as np
from PIL import Image
import os
import matplotlib.pyplot as plt
from io import BytesIO

# ===== 트리 구조 정의 =====
class TreeNode:
    def __init__(self, feature_name=None, threshold=None, label=None):
        self.feature_name = feature_name
        self.threshold = threshold
        self.label = label
        self.left = None
        self.right = None

    def is_leaf(self):
        return self.label is not None

# ===== 간단한 트리 구성 예시 (색상 기반) =====
def build_sample_tree():
    root = TreeNode("Red Intensity", 100)
    root.left = TreeNode("Green Intensity", 80)
    root.right = TreeNode("Blue Intensity", 120)

    root.left.left = TreeNode(label="Tomato_Bacterial_spot")
    root.left.right = TreeNode(label="Tomato_Late_blight")
    root.right.left = TreeNode(label="Tomato_Leaf_Mold")
    root.right.right = TreeNode(label="Tomato_healthy")

    return root

# ===== 이미지에서 RGB 평균 추출 =====
def extract_rgb_features(image):
    image = image.resize((128, 128))
    arr = np.array(image)
    if arr.ndim == 3:
        red_mean = np.mean(arr[:, :, 0])
        green_mean = np.mean(arr[:, :, 1])
        blue_mean = np.mean(arr[:, :, 2])
    else:
        red_mean = green_mean = blue_mean = 0
    return red_mean, green_mean, blue_mean

# ===== 트리 순회 시각화 =====
def traverse_and_visualize(node, features, traversal=[], order='pre'):
    if node is None:
        return

    if order == 'pre':
        traversal.append((node, features))

    if not node.is_leaf():
        value = features.get(node.feature_name, 0)
        if value < node.threshold:
            traverse_and_visualize(node.left, features, traversal, order)
        else:
            traverse_and_visualize(node.right, features, traversal, order)

    if order == 'in':
        traversal.append((node, features))

    if order == 'post':
        traversal.append((node, features))

    return traversal

# ===== Streamlit 앱 구성 =====
st.set_page_config(page_title="Tomato Tree Classifier", layout="wide")
st.title("🌱 이진 탐색 트리 기반 토마토 질병 분류 시각화")

uploaded_file = st.file_uploader("🔍 질병이 의심되는 토마토 잎 이미지를 업로드하세요", type=["jpg", "png", "jpeg"])
order = st.selectbox("트리 순회 방식 선택", ["pre", "in", "post"], index=0)

if uploaded_file:
    image = Image.open(uploaded_file)
    st.image(image, caption="업로드된 이미지", use_column_width=True)

    red_mean, green_mean, blue_mean = extract_rgb_features(image)
    st.write(f"**Red 평균:** {red_mean:.2f}, **Green 평균:** {green_mean:.2f}, **Blue 평균:** {blue_mean:.2f}")

    features = {
        "Red Intensity": red_mean,
        "Green Intensity": green_mean,
        "Blue Intensity": blue_mean
    }

    # 트리 구축 및 순회
    tree = build_sample_tree()
    path = traverse_and_visualize(tree, features, [], order=order)

    # 순회 시각화
    st.subheader("🧭 순회 경로 시각화")
    for idx, (node, feats) in enumerate(path):
        if node.is_leaf():
            st.success(f"[{idx+1}] 예측 결과: **{node.label}**")
        else:
            val = feats[node.feature_name]
            direction = "왼쪽(작음)" if val < node.threshold else "오른쪽(큼)"
            st.info(f"[{idx+1}] `{node.feature_name} < {node.threshold}` → {val:.2f} → {direction}")

else:
    st.info("왼쪽에서 이미지를 업로드하면 분류 흐름이 시작됩니다.")
