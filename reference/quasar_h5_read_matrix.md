# Read an anndata matrix group / dataset (X, raw/X, layers, obsm) as an R matrix

Sparse CSR/CSC groups become a Matrix sparse matrix; dense datasets a
base matrix. Orientation matches the R single-cell convention (features
x observations for X). A dataframe stored where a matrix is expected is
read and coerced to a numeric matrix.

## Usage

``` r
quasar_h5_read_matrix(filename, name)
```

## Arguments

- filename:

  path to the .h5ad file.

- name:

  group / dataset name, e.g. `"X"`, `"raw/X"`, `"obsm/X_pca"`.

## Value

a dense base matrix or a sparse `dgCMatrix`.
