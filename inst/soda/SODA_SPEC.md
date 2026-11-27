# SODA — Single Omics Data Archive

**Status:** schema `1.0.0` (iCellR2 2.x)  
**Working form:** a **directory** (folder). Zip only for transfer.  
**Goal:** a portable, language-friendly on-disk layout for single-cell projects.

## Quick open

```r
library(iCellR2)
obj <- open.soda("MyProject.soda")
validate.soda("MyProject.soda")
```

## Required layout

```text
MyProject.soda/
  index.html             # static archive home page (auto-written)
  inventory.json         # per-file inventory (auto-written on refresh)
  soda.json              # manifest (required)
  soda.info              # key=value twin (optional but written)
  README.md              # human intro (auto-written)
  MANIFEST.txt           # core file list (auto-written)
  architecture.txt       # design notes (optional)
  SODA_SPEC.md           # this spec (copied when available)
  metadata/
    cell_data.tsv        # cells × metadata (required)
    sample_data.tsv
    project_info.json
    version.json
    params.json
    commands.json
    <mod>_feature_data.tsv
    <mod>_gene_stats.tsv
    <mod>_gene_model.tsv
  modalities/
    <mod>/               # e.g. rna, atac, adt
      assays/
        <assay>/         # e.g. counts, data, scale.data
          *.h5           # tenx_h5 or iCellR2_csc_v1 / dense
          assay.json     # optional assay sidecar
      reductions/        # preferred home for PCA/UMAP
      filters/
  reductions/            # legacy / shared reductions (optional)
  cache/                 # disposable derived files
  logs/
```

## `soda.json` contract

| Field | Meaning |
|-------|---------|
| `format` | Must be `"SODA"` |
| `soda_version` | Archive format version (e.g. `1.0.0`) |
| `schema_version` | Schema version (e.g. `1.0.0`) |
| `package_version` | Writer package version |
| `root_dir` | Prefer `"."` (portable). Absolute paths discouraged |
| `project_name` | Human name |
| `modalities` | Map of modality → assays / reductions |
| `metadata_paths` | **Relative** paths to tables |
| `portable` | `true` when paths are relative |
| `created` | Timestamp |

Assay `path` entries and reduction embedding paths must be **relative to the SODA root** so the folder can be moved or uploaded intact.

## Assay formats

| Backend | Description |
|---------|-------------|
| `tenx_h5` | Native 10x HDF5 (`/matrix` CSC-style) + `tenx_meta.json` |
| `hdf5_sparse` | `iCellR2_csc_v1` CSC HDF5 |
| `hdf5_dense` | Dense HDF5 (e.g. skinny `scale.data`) |
| `lazy_norm` | Size factors only; matrix lives in source assay |

Readers should honor `soda.json` + on-disk discovery under `modalities/<mod>/assays/`.



## Inventory

On `refresh.soda()` / manifest write, iCellR2 writes:

| File | Purpose |
|------|---------|
| `inventory.json` | Machine-readable rollup of every substantive file |
| `index.html` | Browser home page summarizing the archive |

HDF5 peeks record format/kind, group, dimensions, nnz (sparse), dtype, and dataset names — without loading full matrices.

## Stability

- **Do not break `schema_version: 1.0.0`** without a migration path.
- Additive fields are fine; renames/removals need a new schema version.
- Prefer TSV/JSON/Parquet for tables — avoid requiring RDS to open a SODA.

## Validation

```r
validate.soda("MyProject.soda")           # report
validate.soda("MyProject.soda", error = TRUE)
refresh.soda(obj)                         # rewrite soda.json + README + MANIFEST
refresh.soda("MyProject.soda")            # path form: relativize + README
```

## Distribution

1. Compute / analyze on the **folder**.
2. For sharing: zip the folder (`MyProject.soda.zip`) or future `pack.soda()`.
3. Recipients unzip → open `index.html` in a browser, or `open.soda("MyProject.soda")` in R.

## Non-R readers (future / community)

A minimal opener only needs to:

1. Parse `soda.json`
2. Resolve relative paths against the folder root
3. Read `metadata/cell_data.tsv` and assay HDF5 with h5py / anndata bridges

---

*SODA = Single Omics Data Archive. Spec owned with iCellR2; layout is intentionally tool-agnostic.*
