# tree_disease_app.py
import streamlit as st
import graphviz

# -------------------
# 트리 구조 정의
# -------------------
class TreeNode:
    def __init__(self, severity, disease, description):
        self.severity = severity
        self.disease = disease
        self.description = description
        self.left = None
        self.right = None

class BST:
    def __init__(self):
        self.root = None

    def insert(self, severity, disease, description):
        def _insert(node, severity, disease, description):
            if node is None:
                return TreeNode(severity, disease, description)
            if severity < node.severity:
                node.left = _insert(node.left, severity, disease, description)
            else:
                node.right = _insert(node.right, severity, disease, description)
            return node
        self.root = _insert(self.root, severity, disease, description)

    def inorder(self):
        result = []
        def _inorder(node):
            if node:
                _inorder(node.left)
                result.append((node.disease, node.severity))
                _inorder(node.right)
        _inorder(self.root)
        return result

    def preorder(self):
        result = []
        def _preorder(node):
            if node:
                result.append((node.disease, node.severity))
                _preorder(node.left)
                _preorder(node.right)
        _preorder(self.root)
        return result

    def postorder(self):
        result = []
        def _postorder(node):
            if node:
                _postorder(node.left)
                _postorder(node.right)
                result.append((node.disease, node.severity))
        _postorder(self.root)
        return result

    def visualize(self):
        dot = graphviz.Digraph()
        def add_nodes(node):
            if node:
                dot.node(f"{node.severity}", f"{node.disease}\n(위험도 {node.severity})")
                if node.left:
                    dot.edge(f"{node.severity}", f"{node.left.severity}")
                    add_nodes(node.left)
                if node.right:
                    dot.edge(f"{node.severity}", f"{node.right.severity}")
                    add_nodes(node.right)
        add_nodes(self.root)
        return dot

# -------------------
# 샘플 데이터 (질병 + 위험도)
# -------------------
sample_diseases = [
    (2, "Tomato_healthy", "정상적인 토마토 잎입니다."),
    (5, "Tomato_Leaf_Mold", "잎에 곰팡이가 생깁니다. 습도 조절 필요."),
    (3, "Tomato_Early_blight", "잎 끝이 검게 마르며 퍼집니다."),
    (6, "Tomato_Septoria_leaf_spot", "작은 점들이 나타나며 병반이 커집니다."),
    (7, "Tomato_Spider_mites", "해충에 의해 점점 황변되고 말라갑니다."),
    (9, "Tomato_YellowLeaf_Curl_Virus", "심각한 바이러스성 질병입니다."),
]

# -------------------
# 웹앱 UI
# -------------------
st.title("🌿 토마토 질병 진단 트리 시각화 웹앱")
st.markdown("질병 진단 결과를 이진 탐색 트리로 구조화하고 시각적으로 탐색합니다.")

bst = BST()
for sev, name, desc in sample_diseases:
    bst.insert(sev, name, desc)

# 트리 시각화
st.subheader("🔎 트리 구조 시각화")
st.graphviz_chart(bst.visualize())

# 순회 결과
st.subheader("📋 트리 순회 결과")
col1, col2, col3 = st.columns(3)
with col1:
    st.markdown("**전위 순회** (Root → Left → Right)")
    st.write([d[0] for d in bst.preorder()])
with col2:
    st.markdown("**중위 순회** (Left → Root → Right)")
    st.write([d[0] for d in bst.inorder()])
with col3:
    st.markdown("**후위 순회** (Left → Right → Root)")
    st.write([d[0] for d in bst.postorder()])

# 질병 상세 정보
st.subheader("💡 질병 상세 보기")
selected = st.selectbox("질병을 선택하세요", [d[0] for d in bst.inorder()])
desc = next((d[2] for d in sample_diseases if d[1] == selected), "정보 없음")
st.info(f"**{selected} 설명:** {desc}")
