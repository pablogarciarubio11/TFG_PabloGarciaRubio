import scanpy as sc
import anndata as ad
import scib
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

input_path = "adatas_raw_noncoding_logcounts_before_baw.h5ad"
output_path = "adata_noncoding_baw_all_genes_umap_leiden.h5ad"

# Carpeta para guardar figuras
os.makedirs("Figures_noncoding", exist_ok=True)
sc.settings.figdir = "Figures_noncoding"

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
adatas_raw_noncoding_logcounts_before_baw2 = sc.read(input_path)

print(adatas_raw_noncoding_logcounts_before_baw2)
xmin, xmax = safe_min_max(adatas_raw_noncoding_logcounts_before_baw2.X)
print(xmin, xmax)
print(adatas_raw_noncoding_logcounts_before_baw2.layers.keys())
print(adatas_raw_noncoding_logcounts_before_baw2.obs.columns)

print("Matriz logcounts:")
print(
    adatas_raw_noncoding_logcounts_before_baw2.layers["logcounts"].min(),
    adatas_raw_noncoding_logcounts_before_baw2.layers["logcounts"].max()
)

print("Matriz counts:")
print(
    adatas_raw_noncoding_logcounts_before_baw2.layers["counts"].min(),
    adatas_raw_noncoding_logcounts_before_baw2.layers["counts"].max()
)

print("Matriz X:")
xmin, xmax = safe_min_max(adatas_raw_noncoding_logcounts_before_baw2.X)
print(xmin, xmax)

adatas_raw_noncoding_logcounts_before_baw2.X = safe_copy(
    adatas_raw_noncoding_logcounts_before_baw2.layers["logcounts"]
)

print("Batch key:", batch_key)

# ✅ USAR ESTA CELDA
batch_key = "Sample"
condition_key = "Condition" # (control / lesión-tipo)
label_key = "celltypist" # provisional label para el analisis de parámetros biologicos en integración

start = time.time()

# Escalado batch-aware con datos preanotados
adatas_batch_aware = adatas_raw_noncoding_logcounts_before_baw2.copy()
# scib.preprocessing.scale_batch(adatas_batch_aware, batch_key)
adatas_batch_aware = scib.preprocessing.scale_batch(adatas_batch_aware, batch_key)

end = time.time()
print(f"⏱️ Batch-aware scaling took {end - start:.2f} seconds.")

print("Tipo devuelto por scale_batch:", type(adatas_batch_aware))
print("Shape tras scale_batch:", adatas_batch_aware.shape)

# Dado que scale_batch puede eliminar los metadatos introducidos, los restauro en las siguientes líneas:
adatas_batch_aware = adatas_batch_aware[
    adatas_raw_noncoding_logcounts_before_baw2.obs_names
].copy()

adatas_batch_aware.obs = adatas_raw_noncoding_logcounts_before_baw2.obs.loc[
    adatas_batch_aware.obs_names
].copy()

print("Samples restaurados:")
print(adatas_batch_aware.obs[batch_key].value_counts())

# Dado que también puede generar algún problema en los datos debido al escalado, hay que limpiarlo también:
if sparse.issparse(adatas_batch_aware.X):
    print("NaN antes:", np.isnan(adatas_batch_aware.X.data).sum())
    print("Inf antes:", np.isinf(adatas_batch_aware.X.data).sum())

    adatas_batch_aware.X.data = np.nan_to_num(
        adatas_batch_aware.X.data,
        nan=0.0,
        posinf=0.0,
        neginf=0.0
    )

    adatas_batch_aware.X = adatas_batch_aware.X.astype(np.float32)

    print("NaN después:", np.isnan(adatas_batch_aware.X.data).sum())
    print("Inf después:", np.isinf(adatas_batch_aware.X.data).sum())

else:
    print("NaN antes:", np.isnan(adatas_batch_aware.X).sum())
    print("Inf antes:", np.isinf(adatas_batch_aware.X).sum())

    adatas_batch_aware.X = np.nan_to_num(
        adatas_batch_aware.X,
        nan=0.0,
        posinf=0.0,
        neginf=0.0
    ).astype(np.float32)

    print("NaN después:", np.isnan(adatas_batch_aware.X).sum())
    print("Inf después:", np.isinf(adatas_batch_aware.X).sum())

# Guardar matriz batch-aware
adatas_batch_aware.layers["batch_aware"] = safe_copy(adatas_batch_aware.X)

print(adatas_batch_aware)
print(
    "batch_aware min/max:",
    adatas_batch_aware.layers["batch_aware"].min(),
    adatas_batch_aware.layers["batch_aware"].max()
)
print("layers:", adatas_batch_aware.layers.keys())

# PCA sobre datos corregidos (en .X)
sc.tl.pca(
    adatas_batch_aware,
    use_highly_variable=True
) # Para asegurarte de usar únicamente HVGs si no se puede con un simple: sc.tl.pca(adatas_batch_aware)

# Guardar como X_batch_aware si lo necesitas explícitamente
adatas_batch_aware.obsm["X_BAW"] = adatas_batch_aware.obsm["X_pca"].copy()

# Usar esa representación para vecinos y UMAP
sc.pp.neighbors(adatas_batch_aware, use_rep="X_BAW")
sc.tl.umap(adatas_batch_aware)

# Clustering
start = time.time()
sc.tl.leiden(adatas_batch_aware, resolution=0.5, key_added="best_res")
end = time.time()
print(f"Leiden clustering took {end - start:.2f} seconds.")

adatas_batch_aware.obs["leiden_final"] = adatas_batch_aware.obs["best_res"].copy()

print(adatas_batch_aware.obs["leiden_final"].value_counts().sort_index())

print("Número de clusters:")
print(adatas_batch_aware.obs["leiden_final"].nunique())

sc.pl.umap(
    adatas_batch_aware,
    color=["leiden_final", "Sample", "Condition", "preannotation"],
    wspace=0.4,
    show=False,
    save="_baw_leiden_sample_condition_preannotation.png"
)

adatas_batch_aware.X = safe_copy(adatas_batch_aware.layers["logcounts"])

# ✅ USAR ESTA CELDA
start = time.time()
# --- Prepara la matriz para CellTypist ---
adatas_batch_aware.X = safe_copy(adatas_batch_aware.layers["logcounts"]) # Usar 'logcounts' ya es log1p(normalize_total)
adatas_batch_aware.raw = None   # evita conflictos con raw

# --- Cargar el modelo preentrenado ---
ref_model_Tabula_CellTypist = models.Model.load(
    model="Reference_Pablo_celltypist.pkl"
)
# --- Anotación con CellTypist ---
# predictions = [predict_cells(ad.copy()) for ad in adatas_Global] en el caso que quisieramos hacelo por separado.
# Emplear .copy para no usar la forma normalizada de la función 'def predict_cells(adata)' y eliminar el raw counts. 
# Eso lo haremos sobre una copia de anndataobjects almaceadas en adatas

predictions = celltypist.annotate(
    adatas_batch_aware,
    model=ref_model_Tabula_CellTypist,
    majority_voting=False
)
# Nota: majority_voting=True suaviza las predicciones a nivel de cluster. Alternativa si quieres que cuando se ejecute se vea... majority_voting=False, verbose=True)

# --- Guardar las etiquetas en obs ---
adatas_batch_aware.obs["annotation"] = predictions.predicted_labels
pmat = predictions.probability_matrix  # DataFrame (células × tipos)
adatas_batch_aware.obs["annotation_score"] = pmat.max(axis=1).values # --- Probabilidad máxima por celda

# --- Restaurar counts en .X para integración posterior
# Para el objeto final es más seguro dejar X como logcounts.
# Los counts siguen guardados en layers["counts"].
adatas_batch_aware.X = safe_copy(adatas_batch_aware.layers["logcounts"])
adatas_batch_aware.raw = None

# --- Guardar
adatas_batch_aware.write(
    "adata_noncoding_baw_celltypist_final.h5ad",
    compression="gzip"
)

end = time.time()
print(f"Integrated annotation took {end - start:.2f} seconds.")

sc.pl.umap(
    adatas_batch_aware,
    color=["leiden_final", "annotation"],
    wspace=0.4,
    show=False,
    save="_baw_leiden_annotation.png"
)

sc.pl.umap(
    adatas_batch_aware,
    color=["annotation"],
    wspace=0.4,
    show=False,
    save="_celltypist_annotation.png"
)

sc.pl.umap(
    adatas_batch_aware,
    color=["annotation_score"],
    wspace=0.4,
    show=False,
    save="_celltypist_annotation_score.png"
)

adatas_batch_aware.X = safe_copy(adatas_batch_aware.layers["logcounts"])

adatas_batch_aware.write(
    "adata_noncoding_baw_celltypist_final_final.h5ad",
    compression="gzip"
)

# Guardar también el objeto integrado antes/final con UMAP y Leiden
adatas_batch_aware.write(
    output_path,
    compression="gzip"
)

# Batch-aware scaled representation
# ESTO ESTÁ MAL, EJECUTAR TODO DE NUEVO
# AHORA YA ESTÁ BIEN ESTA PARTE
# adatas_batch_aware.write("d:/TFG/h5ad_files_coding/adatas_raw_wo_batch_aware.h5ad")
# adatas_batch_aware.write("adatas_raw_wo_batch_aware.h5ad")
os.makedirs("Noncoding", exist_ok=True)
adatas_batch_aware.write(
    "Noncoding/adatas_raw_batch_aware.h5ad",
    compression="gzip"
)

print("✅ Terminado.")
print("Objeto final guardado en: adata_noncoding_baw_celltypist_final.h5ad")
print("Objeto final adicional guardado en: adata_noncoding_baw_celltypist_final_final.h5ad")
print("Objeto integrado guardado en:", output_path)
print("Objeto batch-aware guardado en: Noncoding/adatas_raw_batch_aware.h5ad")
print("Figuras guardadas en la carpeta: Figures_noncoding/")