# MAFPC

Multi-Ancestry Fine-Mapping with functional Principal Components: joint fine-mapping across ancestries with functional annotations. *Preliminary R package.*

## Installation

```r
# Install devtools if needed
install.packages("devtools")

# Install from GitHub
devtools::install_github("gchen4422/MAFPC")


# Load the package
library(MAFPC)
```
## Core function

### `mafpc_core()`

```r
mafpc_core <- function(
  R_mat_list,              # named list of LD correlation matrices by ancestry (SNP x SNP)
  summary_stat_list,       # named list of data.frames with SNP, A1, A2, BETA, SE (harmonized, same SNP order)
  L,                       # integer: max number of causal signals
  residual_variance = NULL,
  prior_weights = NULL,    # optional numeric vector length m (per-SNP prior multipliers)
  ancestry_weight = NULL,  # optional named numeric vector (e.g., c(EUR = 0.6, AFR = 0.4))
  optim_method = "optim",
  estimate_residual_variance = FALSE,
  max_iter = 100,
  cor_method = "min.abs.corr",
  cor_threshold = 0.5,
  annot = NULL,            # optional SNP x feature matrix/data.frame aligned to SNP order
  annot_method = NULL,     # e.g., "mlk" for sparsePCA, "glm", or NULL
  est_annot_prior = "fixed"
)
```

#### Arguments

- **`R_mat_list`**  
  Named list of ancestry-specific **LD correlation matrices** (symmetric; identical SNP IDs and identical SNP order across ancestries).

- **`summary_stat_list`**  
  Named list of GWAS summary-stat data.frames aligned to LD. Required columns: `SNP`, `A1`, `A2`, `BETA`, `SE` (optionally `CHR`, `BP`). Alleles should be uppercase and harmonized across ancestries.

- **`L`**  
  Integer: maximum number of causal signals to model for the region.

- **`residual_variance`**  
  Numeric or `NULL`. If `NULL`, use default or estimate when `estimate_residual_variance = TRUE`.

- **`prior_weights`**  
  Optional numeric vector of length *m* (per-SNP prior multipliers), aligned to SNP order.

- **`ancestry_weight`**  
  Optional **named** numeric vector weighting ancestries in the joint objective (e.g., `c(EUR = 0.6, AFR = 0.4)`).

- **`optim_method`**  
  Optimizer choice (default `"optim"`).

- **`estimate_residual_variance`**  
  Logical. If `TRUE`, estimate residual variance.

- **`max_iter`**  
  Integer: maximum optimization iterations.

- **`cor_method`, `cor_threshold`**  
  Internal correlation/LD-consistency handling (defaults `"min.abs.corr"` and `0.5`).

- **`annot`**  
  Annotation matrix/data.frame (SNP × features), rows aligned to SNP order.

- **`annot_method`**  
  How to use annotations (e.g., `"mlk"` for sparse PCA, `"glm"`, or `NULL`).

- **`est_annot_prior`**  
  Annotation-prior mode (e.g., `"fixed"`).


