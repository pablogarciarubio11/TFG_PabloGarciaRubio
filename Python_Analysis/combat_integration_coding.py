import scanpy as sc
import anndata as ad
import gc
import time
import os
import matplotlib
matplotlib.use("Agg")

import numpy as np
from scipy import sparse

import celltypist
from celltypist import models

batch_key = "Sample"

input_path = "adatas_raw_coding_logcounts_before_combat.h5ad"
output_path = "adata_coding_combat_all_genes_umap_leiden.h5ad"

# Carpeta para guardar figuras
os.makedirs("Figures_coding", exist_ok=True)
sc.settings.figdir = "Figures_coding"

def safe_min_max(x):
    return x.min(), x.max()

def safe_copy(x):
    return x.copy()

def safe_astype(x, dtype):
    if sparse.issparse(x):
        return x.astype(dtype)
    else:
        return np.asarray(x).astype(dtype)


print("Leyendo objeto original...")
# adatas_raw = sc.read_h5ad(input_path)
adatas_raw = sc.read(input_path)

print(adatas_raw)
xmin, xmax = safe_min_max(adatas_raw.X)
print(xmin, xmax)
print(adatas_raw.layers.keys())
print(adatas_raw.obs.columns)

print("Matriz logcounts:")
print(adatas_raw.layers["logcounts"].min(), adatas_raw.layers["logcounts"].max())

print("Matriz counts:")
print(adatas_raw.layers["counts"].min(), adatas_raw.layers["counts"].max())

print("Matriz X:")
xmin, xmax = safe_min_max(adatas_raw.X)
print(xmin, xmax)
adatas_raw.X = adatas_raw.layers["logcounts"].copy()

print("Batch key:", batch_key)

# ✅ USAR ESTA CELDA
batch_key="Sample"

start = time.time()

# No uso scib.integration.combat porque hace una copia completa del objeto
# y puede provocar MemoryError. Uso scanpy directamente.
sc.pp.combat(adatas_raw, key=batch_key)

adatas_full_combat = adatas_raw

end = time.time()
print(f"⏱️ TIME: {end - start:.2f} seconds.")

# Guardar matriz corregida
adatas_full_combat.X = safe_astype(adatas_full_combat.X, "float32")
adatas_full_combat.layers["combat"] = safe_copy(adatas_full_combat.X)

# ⏱️ TIME: ~80 seconds.

print(adatas_full_combat)
print("combat min/max:", adatas_full_combat.layers["combat"].min(), adatas_full_combat.layers["combat"].max())
print("layers:", adatas_full_combat.layers.keys())

# Usar la matriz corregida por ComBat como matriz activa
adatas_full_combat.X = adatas_full_combat.layers["combat"].copy()

sc.tl.pca(
    adatas_full_combat,
    n_comps=50,
    svd_solver="arpack"
)

adatas_full_combat.obsm["X_pca_combat"] = adatas_full_combat.obsm["X_pca"].copy()
sc.pp.neighbors(
    adatas_full_combat,
    n_neighbors=15,
    n_pcs=30,
    use_rep="X_pca_combat"
)

sc.tl.umap(adatas_full_combat)

sc.tl.leiden(
    adatas_full_combat,
    resolution=0.5,
    key_added="leiden_combat_0_5"
)

adatas_full_combat.obs["leiden_final"] = adatas_full_combat.obs["leiden_combat_0_5"].copy()

print(adatas_full_combat.obs["leiden_final"].value_counts().sort_index())

print("Número de clusters:")
print(adatas_full_combat.obs["leiden_final"].nunique())

sc.pl.umap(
    adatas_full_combat,
    color=["leiden_final", "Sample", "Condition", "preannotation"],
    wspace=0.4,
    show=False,
    save="_combat_leiden_sample_condition_preannotation.png"
)

adatas_full_combat.X = adatas_full_combat.layers["logcounts"].copy()

# ✅ USAR ESTA CELDA
start = time.time()
# --- Prepara la matriz para CellTypist ---
adatas_full_combat.X = adatas_full_combat.layers["logcounts"].copy() # Usar 'logcounts' ya es log1p(normalize_total)
adatas_full_combat.raw = None   # evita conflictos con raw

# --- Cargar el modelo preentrenado ---
ref_model_Tabula_CellTypist = models.Model.load(
    model="Reference_Pablo_celltypist.pkl"
)
# --- Anotación con CellTypist ---
# predictions = [predict_cells(ad.copy()) for ad in adatas_Global] en el caso que quisieramos hacelo por separado.
# Emplear .copy para no usar la forma normalizada de la función 'def predict_cells(adata)' y eliminar el raw counts. 
# Eso lo haremos sobre una copia de anndataobjects almaceadas en adatas

predictions = celltypist.annotate(
    adatas_full_combat,
    model=ref_model_Tabula_CellTypist,
    majority_voting=False
)
# Nota: majority_voting=True suaviza las predicciones a nivel de cluster. Alternativa si quieres que cuando se ejecute se vea... majority_voting=False, verbose=True)

# --- Guardar las etiquetas en obs ---
adatas_full_combat.obs["annotation"] = predictions.predicted_labels
pmat = predictions.probability_matrix  # DataFrame (células × tipos)
adatas_full_combat.obs["annotation_score"] = pmat.max(axis=1).values # --- Probabilidad máxima por celda

# --- Restaurar counts en .X para integración posterior
# Para el objeto final es más seguro dejar X como logcounts.
# Los counts siguen guardados en layers["counts"].
adatas_full_combat.X = adatas_full_combat.layers["logcounts"].copy()
adatas_full_combat.raw = None

# --- Guardar
adatas_full_combat.write('adata_coding_combat_celltypist_final.h5ad', compression="gzip")

end = time.time()
print(f"Integrated annotation took {end - start:.2f} seconds.")

sc.pl.umap(
    adatas_full_combat,
    color=["leiden_final", "annotation"],
    wspace=0.4,
    show=False,
    save="_combat_leiden_annotation.png"
)

sc.pl.umap(
    adatas_full_combat,
    color=["annotation"],
    wspace=0.4,
    show=False,
    save="_celltypist_annotation.png"
)

sc.pl.umap(
    adatas_full_combat,
    color=["annotation_score"],
    wspace=0.4,
    show=False,
    save="_celltypist_annotation_score.png"
)

adatas_full_combat.X = adatas_full_combat.layers["logcounts"].copy()

adatas_full_combat.write(
    "adata_coding_combat_celltypist_final_final.h5ad",
    compression="gzip"
)

# Guardar también el objeto integrado antes/final con UMAP y Leiden
adatas_full_combat.write(
    output_path,
    compression="gzip"
)

print("✅ Terminado.")
print("Objeto final guardado en: adata_coding_combat_celltypist_final.h5ad")
print("Objeto final adicional guardado en: adata_coding_combat_celltypist_final_final.h5ad")
print("Objeto integrado guardado en:", output_path)
print("Figuras guardadas en la carpeta: Figures/")
