# Read an anndata dataframe group (obs / var) into a data.frame

Handles the modern dataframe encoding (per-column categorical groups,
nullable/masked integer & boolean groups, string and numeric arrays) as
well as two legacy factor layouts (a shared `__categories` group and
`/uns/<col>_categories`). The index dataset becomes the rownames.

## Usage

``` r
quasar_h5_read_dataframe(filename, name, index_as_column = FALSE)
```

## Arguments

- filename:

  path to the .h5ad file.

- name:

  group name, e.g. `"obs"`, `"var"` or `"raw/var"`.

- index_as_column:

  keep the index as a column in addition to using it as rownames.
  Default `FALSE` (cleaner colData / meta.data).

## Value

a data.frame with rownames taken from the anndata index.
