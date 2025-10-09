# (2) Starter notebook outline: simulation → descriptors → DR/cluster → plots
# Tip: keep functions in a .py module (src/descriptors.py), import here.

# --- 0. Imports & setup ---
import json, pathlib, numpy as np, pandas as pd
from pathlib import Path
from datetime import datetime

# math & stats
from scipy.fft import fft
from scipy.signal import find_peaks
from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA
from sklearn.cluster import KMeans
from sklearn.mixture import GaussianMixture
from sklearn.metrics import silhouette_score
import umap  # pip install umap-learn

# viz
import matplotlib.pyplot as plt

ROOT = Path("./vesicle_ml")  # adjust to your dataset root
PROFILE_DIR = ROOT / "profiles"
FEATURES_PATH = ROOT / "features" / "features.parquet"
META_DIR = ROOT / "metadata"
FIG_DIR = ROOT / "figures"
FIG_DIR.mkdir(parents=True, exist_ok=True)

# --- 1. Load data ---
df = pd.read_parquet(FEATURES_PATH)  # one row per solution
ids = df["solution_id"].tolist()

# optional: filter by solver status, constraints, etc.
df = df.query("constraints_ok == 1").reset_index(drop=True)

# --- 2. Build feature matrix (choose a view) ---
# Geometry-only baseline (robust for unsupervised discovery)
GEOM_COLS = [
    # spectral
    *[f"fd_re_{k}" for k in range(1,21)],
    *[f"fd_im_{k}" for k in range(1,21)],
    # differential geometry summaries
    "H_mean","H_std","H_skew","H_kurt","H_L1","H_L2","H_pos_area","H_neg_area",
    "K_mean","K_std","K_L1","K_L2","HK_corr",
    # landmarks/topology
    "rb","psib","H_jump_proxy","curvature_min","curvature_max",
    "curvature_min_loc","curvature_max_loc","num_lobes","has_invagination",
    "tubulation_score","asymmetry_score"
]
X = df[GEOM_COLS].to_numpy()

# Auxiliary features for annotation (not used in clustering)
AUX = df[["E_total","E_bend","E_line","P_osm","x1","H0_1","H0_2","kappa_ratio","sigma"]].copy()

# --- 3. Scale ---
scaler = StandardScaler()
Xz = scaler.fit_transform(X)

# --- 4. Dimensionality reduction ---
pca = PCA(n_components=10, random_state=0)
Xp = pca.fit_transform(Xz)

um = umap.UMAP(n_neighbors=30, min_dist=0.2, n_components=2, random_state=0)
Xu = um.fit_transform(Xz)

# --- 5. Clustering (k sweep + selection) ---
def kmeans_sweep(X, k_range=range(3,13), n_init=20, random_state=0):
    results = []
    for k in k_range:
        km = KMeans(n_clusters=k, n_init=n_init, random_state=random_state)
        lbl = km.fit_predict(X)
        sil = silhouette_score(X, lbl)
        results.append((k, sil, km.inertia_))
    return pd.DataFrame(results, columns=["k","silhouette","inertia"])

sweep = kmeans_sweep(Xp, k_range=range(3,13))
best_k = int(sweep.sort_values("silhouette", ascending=False).iloc[0]["k"])

km = KMeans(n_clusters=best_k, n_init=50, random_state=0)
labels = km.fit_predict(Xp)
df["cluster"] = labels

# Optional alternative: GMM or HDBSCAN
#gmm = GaussianMixture(n_components=best_k, covariance_type="full", random_state=0).fit(Xp)
#labels = gmm.predict(Xp)

# --- 6. Visualizations ---
# 6a) PCA scatter
plt.figure(figsize=(6,5))
plt.scatter(Xp[:,0], Xp[:,1], c=labels, s=12)
plt.xlabel("PCA1"); plt.ylabel("PCA2"); plt.title("PCA clusters")
plt.tight_layout(); plt.savefig(FIG_DIR/"pca_clusters.png", dpi=180)

# 6b) UMAP scatter
plt.figure(figsize=(6,5))
plt.scatter(Xu[:,0], Xu[:,1], c=labels, s=12)
plt.xlabel("UMAP1"); plt.ylabel("UMAP2"); plt.title("UMAP clusters")
plt.tight_layout(); plt.savefig(FIG_DIR/"umap_clusters.png", dpi=180)

# --- 7. Cluster inspection: medoids & representative profiles ---
def cluster_medoids(Xemb, labels):
    medoid_idx = []
    for c in np.unique(labels):
        idx = np.where(labels==c)[0]
        Xc = Xemb[idx]
        # pick point minimizing sum of distances within cluster
        D = np.linalg.norm(Xc[:,None,:]-Xc[None,:,:], axis=2)
        j = idx[np.argmin(D.sum(axis=1))]
        medoid_idx.append(j)
    return medoid_idx

medoid_idx = cluster_medoids(Xp, labels)
df_med = df.iloc[medoid_idx].copy()

# Plot a few medoid profiles
def load_profile(solution_id):
    arr = np.load(PROFILE_DIR / f"{solution_id}.npz")
    return {k: arr[k] for k in arr.files}

for _, row in df_med.iterrows():
    prof = load_profile(row.solution_id)
    r, z = prof["r"], prof["z"]
    plt.figure(figsize=(4,5))
    plt.plot(r, z, lw=2)
    plt.gca().set_aspect('equal', adjustable='box')
    plt.title(f"Cluster {row.cluster} • {row.solution_id}")
    plt.tight_layout()
    plt.savefig(FIG_DIR / f"medoid_{row.cluster}.png", dpi=180)

# --- 8. Parameter-plane maps (dominant cluster over H0 grid) ---
# Assumes many points over (H0_1, H0_2). Create a coarse heatmap of dominant cluster.
df_grid = df.copy()
# bin the plane
bins = 30
df_grid["H1_bin"] = pd.cut(df_grid["H0_1"], bins=bins)
df_grid["H2_bin"] = pd.cut(df_grid["H0_2"], bins=bins)
dom = (df_grid
       .groupby(["H1_bin","H2_bin","cluster"])
       .size()
       .reset_index(name="n")
       .sort_values(["H1_bin","H2_bin","n"], ascending=[True,True,False]))
dom = dom.drop_duplicates(subset=["H1_bin","H2_bin"], keep="first")

# simple image-like plot (index bins)
H1_labels = {k:i for i,k in enumerate(sorted(dom["H1_bin"].unique(), key=lambda iv: iv.left))}
H2_labels = {k:i for i,k in enumerate(sorted(dom["H2_bin"].unique(), key=lambda iv: iv.left))}
M = np.full((len(H1_labels), len(H2_labels)), fill_value=-1)
for _, r_ in dom.iterrows():
    i = H1_labels[r_["H1_bin"]]; j = H2_labels[r_["H2_bin"]]
    M[i,j] = int(r_["cluster"])

plt.figure(figsize=(6,5))
plt.imshow(M.T, origin="lower", aspect="auto")
plt.title("Dominant cluster over (H0_1, H0_2)")
plt.xlabel("H0_1 bins"); plt.ylabel("H0_2 bins")
plt.colorbar(label="cluster id")
plt.tight_layout(); plt.savefig(FIG_DIR/"dominant_cluster_plane.png", dpi=180)

# --- 9. Hysteresis path visualization (if branch_id or path order is known) ---
# Example: plot energy vs step for forward vs backward along same parameter path
def plot_energy_hysteresis(subset_df, title):
    # subset_df should include columns: step (int), E_total (float), branch_id (str)
    plt.figure(figsize=(6,4))
    for b, g in subset_df.groupby("branch_id"):
        g = g.sort_values("step")
        plt.plot(g["step"], g["E_total"], marker='o', label=b)
    plt.xlabel("Path step"); plt.ylabel("Total energy")
    plt.title(title); plt.legend()
    plt.tight_layout()

# --- 10. Optional: supervised quick labeler for new solutions ---
from sklearn.model_selection import cross_val_score
from sklearn.ensemble import RandomForestClassifier

rf = RandomForestClassifier(n_estimators=400, random_state=0, n_jobs=-1)
cv_scores = cross_val_score(rf, Xp, labels, cv=5, scoring="f1_macro")
rf.fit(Xp, labels)

# feature importances (on PCA space for demo; prefer training on raw GEOM_COLS)
imp = getattr(rf, "feature_importances_", None)

print("PCA variance explained:", pca.explained_variance_ratio_[:5].round(3))
print("Best k by silhouette:", best_k)
print("CV macro-F1 (RF):", np.mean(cv_scores).round(3))
