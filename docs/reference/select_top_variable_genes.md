# Select top variable genes

Internal helper used by \[normalize_bulks()\] to retain the most
variable genes based on row variance.

## Usage

``` r
select_top_variable_genes(mat, top_n = 0, top_perc = 0)
```
