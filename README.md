# Nstar – R Function for Estimating N*

## Description

`Nstar` is an R function for estimating the species accumulation-based index **N\***, based on the `specaccum` function from the **vegan** package (Oksanen et al. 2012).  
It allows estimation using either all available samples or multiple random subsets of a community data matrix. 

> **Note**: The **vegan** package must be installed before using `Nstar`.

## Installation

Save the file `Nstar.R` in your working directory. Then in R:

```r
source("Nstar.R")
```

## Usage

```r
Nstar(datum, nset = 0, rep = 500, ExpAll = FALSE, perms = 1000, method = "exact")
```

### Arguments

- **datum**: A community data set — a matrix with N samples (rows) and S species (columns). Column names should represent species.
- **nset**: Integer. Number of samples to use. If `nset = 0`, all samples are used.
- **rep**: Integer. Number of random subsets to use when `nset > 0`.
- **ExpAll**: Logical. If `TRUE`, returns all N\* estimates; otherwise, returns only summary stats.
- **perms**: Integer. Number of permutations (only used when `method = "random"`).
- **method**: Species accumulation method. Default is `"exact"` (see `specaccum` documentation).

## Output

A matrix or list of N\*-related metrics, including:
- `N*`: Estimated number of samples for asymptotic richness
- `N*_Rnd`: Estimate under a random matrix
- `NL`, `NU`: Lower and upper bounds for `N*`
- `mean_alpha`: Mean species per sample
- `SN*`: Estimated number of species at `N*`
- Whittaker and Harrison diversity indices

## Examples

### Example 1: Use full dataset
```r
data(dune, package = "vegan")
Nstar(dune)
```

### Example 2: One subset of 4 samples
```r
Nstar(dune, 4, 1)
```

### Example 3: 100 random subsets of 4 samples
```r
Nstar(dune, 4, 100)
```

> Use `warnings()` to view any warnings issued during execution.

## Authors

Maria Lazarina, Vasiliki Sgardeli, Athanasios S. Kallimanis, Stefanos P. SgardelisStefanos Sgardelis  

## Publication

This function was published in:

Lazarina, M., Sgardeli, V., Kallimanis, A. S., & Sgardelis, S. P. (2013).  
**An effort‐based index of beta diversity**. *Methods in Ecology and Evolution*, 4(3), 217–225.  
https://doi.org/10.1111/2041-210X.12005

## References

Oksanen, J., Blanchet, G., Kindt, R., Minchin, P.R., Legendre, P., O’Hara, B., Simpson, G.L., Solymos, P., Stevens, M.H.H. & Wagner, H. (2012) vegan: Community Ecology Package. R package Version 2.0-3. Available at: http://cran.r-project.org/.
