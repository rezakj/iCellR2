# iCellR2

Tools for large-scale single-cell and spatial omics analysis in R. Built around **SODA** archives (HDF5 assays/embeddings, TSV tables, JSON meta) with an R frontend and a C++ core for on-disk sparse storage, chunked HDF5 I/O, and multithreading.

Supports scRNA-seq, scATAC-seq, multiome, scVDJ-seq, CITE-seq/ADT, Visium, Xenium, and related modalities — import, QC, filtering, normalization, dimensionality reduction, clustering, markers, batch correction, annotation, pseudotime, and multimodal integration.

Python companion: **[iCellPy](https://github.com/rezakj/iCellPy)** (same SODA layout; mix R and Python on one archive).

## Quick start

```r
library(iCellR2)

obj <- load10x("sample1")          # or open an existing .soda
obj <- qc.stats(obj)
obj <- cell.filter(obj, max.mito = 0.1, min.genes = 200)
obj <- norm.data(obj)
obj <- gene.stats(obj)
obj <- make.gene.model(obj)
obj <- run.pca(obj)
obj <- run.umap(obj)
obj <- iclust(obj)
obj <- findMarkers(obj)
cluster.plot(obj)
```

## Bash CLI

```bash
# after install, or from the source tree:
./exec/iCellR2 --version
./exec/iCellR2 help
```

## License

GPL-3 — see [LICENSE](LICENSE).

## Links

- Issues: https://github.com/rezakj/iCellR2/issues
- Homepage: https://github.com/rezakj/iCellR2
