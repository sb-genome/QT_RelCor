---
title: 'R Manual: Identifying Comorbid Endophenotype for a Binary Disease Using Clustered Data'
author: "Prasun Panja"
---

# Overview

This document is an R manual for the functions developed to identify comorbid associations between a quantitative trait (QT) and a binary disease trait. The QT can be identified as an endophenotype or quantitative precursor for the disease. Use the built-in table of contents to jump to sections.


# Required packages

```r
if (!requireNamespace("LaplacesDemon", quietly = TRUE)) {
  stop("Please install LaplacesDemon package before running (contains dmvn).")
}
```


# Introduction

This manual documents the functions to perform the followings (each section below contains: Description, Usage, Arguments, Details, Value, and Examples):

- A function to compute numerical differentiation
- A function to obtain the numerical Jacobian of a given function
- Functions to solve a system of equations using Newton-Raphson
- A function to perform a novel clustering based covariate adjustment technique
- Functions to create some basic matrices
- A function to identify comorbid endophenotype or quantitative precursor

---

# numeric_gradient (numeric differentiation of a function)

**Description**

Compute a numeric gradient for a function `f` at point `x`. This implementation handles functions that return either a scalar or a numeric vector. The finite-difference scheme used here is a forward difference per component (user-supplied small `h1`), so choose `h1` carefully.

**Usage**

```r
numeric_gradient <- function(f, x, h1 = 1e-8)
```

**Arguments**

- `f`: a function that accepts a numeric vector `x` and returns either a numeric scalar or a numeric vector.
- `x`: numeric vector at which the gradient is evaluated.
- `h1`: finite difference step used for each component (default `1e-8`).

**Details**

For each index `i` the function evaluates `f(x + h * e_i)` and computes the difference `(f(x + h)[i] - f(x)[i]) / h`. If `f` returns a scalar, the length of `fx` will be 1 and `grad` will be length 1. If `f` returns a vector, `grad` is computed component-wise.

**Value**

Numeric vector of partial derivatives (length equals length of `f(x)`).

**Example**

```r
f1 <- function(x) sum((x - 2)^2)
numeric_gradient(f1, c(1,1))
```


# jacobian (numeric Jacobian calculation of a function)

**Description**

Wrapper around `numeric_gradient` intended to produce a Jacobian matrix. Rows correspond to entries of `f(x)`, columns to entries of `x`.

**Usage**

```r
jacobian <- function(f, x, h = 1e-6)
```

**Arguments**

- `f`: function returning scalar or vector.
- `x`: numeric vector.
- `h`: finite-difference step forwarded to `numeric_gradient`.

**Details**

If `numeric_gradient` returns a simple numeric vector (scalar `f`), `jacobian` coerces it to a 1-row matrix. If `numeric_gradient` already returned matrix-like output, that is returned as-is.

**Value**

A numeric matrix (rows = length of `f(x)`, columns = length of `x`).

**Example**

```r
fvec <- function(x) c(x[1] + x[2], x[1] * x[2])
jacobian(fvec, c(2,3))
```


# slfn (multivariate Newton-Raphson system solver)

**Description**

Solve a system of (nonlinear) equations `f(x) = 0` where `f` returns a numeric vector of length `p`. Uses the Jacobian estimated from `jacobian` and performs Newton updates.

**Usage**

```r
slfn <- function(f, x0, tol = 1e-6, max_iter = 100)
```

**Arguments**

- `f`: function mapping `x` (numeric vector) to numeric vector `f(x)`.
- `x0`: initial numeric guess vector.
- `tol`: tolerance for convergence measured by Euclidean change in `x`.
- `max_iter`: maximum number of Newton iterations.

**Details**

At each iteration the code computes `fx <- f(x)` and the Jacobian `J <- jacobian(f, x)`. The Newton step solves `J * delta = fx` and updates `x <- x - delta`. The routine uses `tryCatch` around `solve()` to capture singular Jacobians and stops with an informative error if `J` is singular.

**Value**

Returns the converged numeric vector `x_new` when the Euclidean change falls below `tol`. Throws an error if iteration limit is exceeded or Jacobian is singular.

**Example**

```r
# Solve x^2 - 2 = 0 (2D toy system: each component identical)
fn <- function(x) c(x[1]^2 - 2, x[2]^2 - 3)
slfn(fn, c(1.4, 1.7))
```


# sl1 (1-D Newton root finder with bracketing)

**Description**

Find a root for a scalar function `m(x)` in 1D. The function first searches for a sign change over a grid between `lower` and `upper`. If none is found, it falls back to the midpoint. Then it runs Newton iterations using the numeric derivative.

**Usage**

```r
sl1 <- function(m, lower = -1, upper = 1, by = 0.01, max_newton_iter = 200)
```

**Arguments**

- `m`: scalar function of one variable.
- `lower`, `upper`: interval limits used to search for a sign change.
- `by`: grid spacing used while searching for sign change.
- `max_newton_iter`: maximum Newton steps.

**Details**

Uses `numeric_gradient` to estimate the derivative at the current guess. If derivative is too small or NA, iteration halts and the current guess is returned.

**Value**

Numeric scalar — the estimated root (or fallback midpoint / last Newton iterate).

**Example**

```r
m <- function(x) x^3 - 2*x - 5
sl1(m, -5, 5, by = 0.1)
```


# incc (intra-cluster correlation coefficient)

**Description**

This is an adaptation of classical Intra-class Correlation Coefficient. Computes a pairwise intra-cluster/ intra-family correlations. Accepts a matrix (rows = families, columns = individuals per family) or a numeric vector which is coerced to a matrix.

**Usage**

```r
incc <- function(x)
```

**Arguments**

- `x`: numeric matrix (or object coercible to matrix).

**Details**

Let `global_mean` be the overall mean across the matrix. The numerator accumulates pairwise products over distinct individual pairs within each row (family). The denominator uses a scaled total-sum-of-squares. The implementation is robust to missing values via `mean(..., na.rm = TRUE)` when computing `global_mean` but note that pairwise loops do not explicitly skip NA entries — passing NA-containing matrices can create NA results.

**Value**

Numeric scalar summary `r`. If denominators collapse to zero or family sizes are 1, returns 0.

**Example**

```r
mat <- matrix(c(1,2,3, 2,3,4), nrow = 2, byrow = TRUE)
incc(mat)
```


# clust_reg (non-parametric covariate adjustment by clustering)

**Description**

Performs a clustering-based non-parametric covariate adjustment of the response. The response for each observation is centered by the mean of its cluster; clustering is performed on the predictor columns via hierarchical clustering through `hclust`.

**Usage**

```r
clust_reg <- function(df, response = 'y', predictor = c('x1','x2'), clust_met = 'average', nclust = NULL)
```

**Arguments**

- `df`: data.frame containing response and predictor columns.
- `response`: column name (or numeric index) for the response variable.
- `predictor`: character vector of predictor column names (or numeric indices). Multiple predictors supported.
- `clust_met`: clustering linkage method passed to `hclust` (default `'average'`).
- `nclust`: desired number of clusters. If `NULL`, a heuristic `max(2, floor(nrow(df) * 0.1))` is used.

**Details**

- The function extracts predictor columns and applies `dist()` followed by `hclust()` and `cutree()` to obtain cluster assignments.
- Cluster means for the response are computed and subtracted from each response to produce non-parametric residuals.

**Value**

A list with components:

- `residuals`: numeric vector of response minus cluster mean.
- `clusters`: integer vector of cluster assignments per row.
- `predicted`: numeric cluster means aligned to rows.

**Example**

```r
set.seed(123)
x1 <- x2 <- x3 <- y <- e <- vector()
a <- 0.3
b <- 0.5
c <- 0.6
for (i in 1:50) {
  x1[i] <- rnorm(1, 2, 2)
  x2[i] <- rchisq(1, 1)
  x3[i] <- rpois(1, 1.3)
  e[i] <- rnorm(1, 0, 1)
  y[i] <- a * x1[i] + b * x2[i] + c * x3[i] + e[i]
}
df_dj <- data.frame(y, x1, x2, x3)
out <- clust_reg(df_dj, response = 'y', predictor = c('x1', 'x2', 'x3), nclust = 7)
head(out$residuals)
```


# idm and jm

**Description**

Convenience constructors:

- `idm(n)`: identity matrix `n x n`.
- `jm(n)`: all-ones matrix `n x n`.

**Usage**

```r
idm <- function(n) diag(1, n, n)
jm  <- function(n) matrix(1, nrow = n, ncol = n)
```

**Example**

```r
idm(3)
jm(3)
```


# cmrbd (hypothesis testing for location parameter based on clustered data)

**Description**

`cmrbd` implements a complex family-based test comparing a quantitative variable across affected/unaffected statuses using a model that accounts for within-family correlation. It estimates parameters (means, variance, intra-family correlation `r`) under a null model (no case/control mean difference) and an alternative model (separate means for cases/controls), computes a likelihood ratio test and a permutation-based non-parametric statistic.

**Usage**

```r
cmrbd <- function(df, value_col = "value", affected_label = "Affected", n_perm = 1000)
```

**Arguments**

- `df`: data.frame — **must** contain `Family_Number` and `Category` columns. `Category` encodes affected status.
- `value_col`: column name or index of the quantitative variable to test.
- `affected_label`: the label in `Binary Trait` corresponding to affected individuals (default: `'Affected'`). For the code above, numeric `1` is accepted in examples.
- `n_perm`: number of permutations used for non-parametric inference (default 1000).

**Details**

- The function groups observations by family, constructs per-family vectors of `x` (quantitative values) and binary `y` (affected indicator), and removes families that are non-informative.
- The function tests in both parametric and non-parametric approach
- It iteratively estimates `r` (intra-family correlation) and mean/variance under H0 and H1 using internal root-finding (`sl1` and `slfn`) and closed-form updates (based on normal likelihood).
- Log-likelihoods under both models are computed using `LaplacesDemon::dmvn()` with covariance matrices of the form `var * ((1 - r) I + r J)`.
- The LRT statistic (`-2*(lhd0 - lhd1)`) and its chi-square p-value are returned.
- The non-parametric approach divided into two parts
- First part do a Wilcoxo rank-sum test to compare between independent `Affected` individuals with Family based `Unaffected`
- The second part do a Wilcoxon signed-rank test between the mean of cases and controls of a family
- The two test statistics are summed and the p-value is obtained by permutations

**Value**

A `data.frame` with columns:

- `lhd0`: log-likelihood under H0
- `lhd1`: log-likelihood under H1
- `lrt`: likelihood ratio test statistic
- `p_val`: chi-square p-value for `lrt`
- `D`: qualitative decision from parametric p-value (`Significant` / `Insignificant`)
- `T`: observed permutation statistic (sum of two Wilcoxon-type pieces)
- `p_np`: permutation p-value
- `Dnp`: qualitative decision from permutation (`Significant` / `Insignificant`)

**Important implementation notes and caveats**

- The code uses nested loops and `apply` patterns; for very large family sets and/or large `n_perm` this can be slow.
- Missing data handling: the function converts values to `numeric` but does not comprehensively handle `NA` entries everywhere. You may need to pre-clean `df` (remove or impute `NA` values) before calling `cmrbd`.
- Stability of root-finding for `r` is guarded by limiting `r` search to (-0.999, 0.999) and by fallbacks to `0` when the numeric solver fails.

**Examples**

```r
# toy example (construct a small sib-ship-like data.frame)
df <- data.frame(
  Family_Number = c(1,1,2,2,2,3,3),
  Category = c(1,0,1,0,0,1,0),
  Trait_val = c(2.3, 1.8, 3.2, 2.9, 2.5, 4.1, 3.9)
)
res <- cmrbd(df, value_col = "Trait_val", affected_label = 1, n_perm = 200)
print(res)
```


# Final notes

- This manual documents the logic and usage for the functions included. This work is part of PhD Thesis of Prasun Panja and Supervised by Dr. Samsiddhi Bhattacharjee.



