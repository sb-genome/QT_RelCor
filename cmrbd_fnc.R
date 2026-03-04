# Required packages
if (!requireNamespace("LaplacesDemon", quietly = TRUE)) {
  stop("Please install LaplacesDemon package before running (contains dmvn).")
}

# ---------------------------
# Numeric gradient (central difference)
# Works for scalar or vector-valued f.
# h default chosen small but not too small.
numeric_gradient <- function(f, x, h1 = 1e-8) {
  fx <- f(x)
  grad <- numeric(length(fx))
  for (i in seq_along(fx)) {
    h <- numeric(length(x))
    h[i] <- h1
    grad[i] <- (f(x + h)[i] - fx[i]) / h1
  }
  return(grad)
}

# ---------------------------
# Jacobian: wrapper that returns matrix (rows = length(f(x)), cols = length(x))
jacobian <- function(f, x, h = 1e-6) {
  g <- numeric_gradient(f, x, h = h)
  # If g is vector => scalar f => return as column matrix
  if (is.vector(g) && is.null(dim(g))) {
    # numeric_gradient returned gradient (for scalar f)
    return(matrix(g, nrow = 1, ncol = length(g), byrow = TRUE))
  } else {
    return(g) # already matrix (rows = length(fx), cols = length(x))
  }
}

# ---------------------------
# Newton-Raphson system solver (multivariate)
# f: function(x) -> numeric vector (length p)
# x0: initial guess vector (length p)
slfn <- function(f, x0, tol = 1e-6, max_iter = 100) {
  x <- as.numeric(x0)
  for (iter in seq_len(max_iter)) {
    fx <- f(x)
    J <- jacobian(f, x)
    # Ensure J is square (p x p)
    if (nrow(J) != ncol(J)) {
      if (nrow(J) == length(x)) {
        # typical: f returns length p = length(x) so J is p x n; ensure it's square
        # proceed only if square
        # else try to return NA
      }
    }
    # try solve; use tryCatch to capture singularity
    s <- tryCatch({
      delta <- solve(J, fx)
      x_new <- x - as.numeric(delta)
      list(ok = TRUE, x_new = x_new)
    }, error = function(e) {
      list(ok = FALSE, msg = conditionMessage(e))
    })
    if (!s$ok) stop("Jacobian is singular or ill-conditioned: ", s$msg)
    if (sqrt(sum((s$x_new - x)^2)) < tol) return(s$x_new)
    x <- s$x_new
  }
  stop("slfn: did not converge within max_iter")
}

# ---------------------------
# 1-D Newton root finder (robusted from sl1)
# m: scalar function m(x)
# range: numeric vector c(lower, upper) where we search for sign change
sl1 <- function(m, lower = -1, upper = 1, by = 0.01, max_newton_iter = 200) {
  # Create bracket grid
  guess <- seq(lower, upper, by = by)
  # find interval where sign changes
  br_found <- FALSE
  for (o in seq_len(length(guess) - 1)) {
    vl1 <- m(guess[o])
    vl2 <- m(guess[o + 1])
    if (!is.na(vl1) && !is.na(vl2) && vl1 * vl2 <= 0) {
      a <- guess[o]; b <- guess[o + 1]; br_found <- TRUE; break
    }
  }
  if (!br_found) {
    # fallback: try midpoint of full interval
    a <- (lower + upper) / 2
  }
  # initial guess
  x0 <- (a + b) / 2
  # Newton iterations (1D)
  for (iter in seq_len(max_newton_iter)) {
    g <- numeric_gradient(function(xx) m(xx), x0)
    if (is.na(g) || abs(g) < .Machine$double.eps) break
    x1 <- x0 - m(x0) / g
    if (is.nan(x1)) break
    if (abs(x1 - x0) < 1e-7) return(x1)
    x0 <- x1
  }
  return(x0)
}

# ---------------------------
# incc: corrected to compute pairwise in-family covariance-like measure
# Accepts matrix or numeric vector; if matrix, rows = families, columns = individuals? 
# I keep your original formula intent but make it robust.
incc <- function(x) {
  x <- as.matrix(x)
  nrow_x <- nrow(x); ncol_x <- ncol(x)
  global_mean <- mean(x, na.rm = TRUE)
  sc <- 0
  s <- 0
  for (i in seq_len(nrow_x)) {
    for (j in seq_len(ncol_x)) {
      for (k in seq_len(ncol_x)) {
        if (j != k) sc <- sc + (x[i, j] - global_mean) * (x[i, k] - global_mean)
      }
    }
  }
  for (i in seq_len(nrow_x)) {
    for (j in seq_len(ncol_x)) {
      s <- s + (x[i, j] - global_mean)^2
    }
  }
  len <- length(x)
  s2 <- ifelse(len == 0, 0, s / len)
  if (s2 != 0 && ncol_x > 1) {
    r <- (1 / (len * (ncol_x - 1) * s2)) * sc
  } else {
    r <- 0
  }
  return(r)
}

########### Non-parametric Covariate Adjustment of Data (response centered by cluster) ###########

clust_reg <- function(df, response = 'y', predictor = c('x1','x2'), clust_met = 'average', nclust = NULL) {
  
  # Resolve response column
  if (is.numeric(response)) {
    response_col_idx <- response
  } else {
    response_col_idx <- which(names(df) == response)
  }
  
  # Resolve predictor columns (can be multiple)
  if (is.numeric(predictor)) {
    predictor_col_idx <- predictor
  } else {
    predictor_col_idx <- match(predictor, names(df))
  }
  
  # Safety checks
  if (length(response_col_idx) != 1) stop("Response column ambiguous or not found")
  if (any(is.na(predictor_col_idx))) stop("Some predictor columns not found")
  if (response_col_idx %in% predictor_col_idx) stop("Response and predictors must be different")
  if (length(predictor_col_idx) >= ncol(df)) stop("Trivial number of predictors")
  
  # Set number of clusters if not provided
  if (is.null(nclust)) {
    nc <- max(2, floor(nrow(df) * 0.1))  # heuristic: 10% of samples
  } else if (nclust < 2 || nclust == nrow(df)) {
    stop("Trivial number of clusters for covariate adjustment")
  } else {
    nc <- nclust
  }
  
  # Clustering on predictor values (multivariate)
  df2 <- df[, predictor_col_idx, drop = FALSE]
  h1 <- hclust(dist(df2), method = clust_met)
  ct <- cutree(h1, k = nc)
  
  # Cluster means of the response
  cluster_means <- tapply(df[, response_col_idx], ct, mean)
  
  # Map each response to its cluster mean
  resp_means <- cluster_means[ct]
  
  # Nonparametric residual (response adjusted by its cluster mean)
  np_res <- df[, response_col_idx] - resp_means
  
  return(list(
    residuals = np_res,
    clusters = ct,
    predicted = resp_means
  ))
}

# ---------------------------
# Simple matrix constructors
idm <- function(n) diag(1, n, n)
jm  <- function(n) matrix(1, nrow = n, ncol = n)

# ---------------------------
# cmrbd function
# df: data.frame with columns Family_Number, Category, and value_col (name or index)
# value_col: name or index of quantitative variable
# affected_label: value in Category indicating affected (default "Affected")
cmrbd <- function(df, value_col = "value", affected_label = "Affected", n_perm = 1000) {
  # defensive checks
  if (!("Family_Number" %in% names(df))) stop("df must contain 'Family_Number' column")
  if (!("Category" %in% names(df))) stop("df must contain 'Category' column")
  # allow either name or index for value_col
  if (is.character(value_col) && !(value_col %in% names(df))) {
    stop("value_col name not found in df")
  }
  if (is.numeric(value_col)) value_col_idx <- value_col else value_col_idx <- which(names(df) == value_col)
  if (length(value_col_idx) != 1) stop("value_col ambiguous")
  
  # build family-wise lists
  famid <- unique(df$Family_Number)
  nfam <- length(famid)
  y <- vector("list", nfam)  # category (0/1) per family
  x <- vector("list", nfam)  # numeric values per family
  
  for (j in seq_along(famid)) {
    sel <- which(df$Family_Number == famid[j])
    cats <- as.character(df$Category[sel])
    vals <- as.numeric(df[[value_col_idx]][sel])
    y[[j]] <- ifelse(cats == affected_label, 1, 0)
    x[[j]] <- vals
  }
  
  # drop families with no informative variability in y (all 0 or all 1) but only if family size >1
  keep_idx <- sapply(seq_len(nfam), function(k) {
    mk <- length(x[[k]])
    if (mk <= 1) return(TRUE)
    ok <- !(all(y[[k]] == 1) || all(y[[k]] == 0))
    return(ok)
  })
  if (!all(keep_idx)) {
    x <- x[keep_idx]; y <- y[keep_idx]; famid <- famid[keep_idx]; nfam <- length(famid)
  }
  
  # m[k] lengths
  m <- sapply(x, length)
  if (any(is.na(m))) m[is.na(m)] <- 0
  
  # Prepare fam_num: how many families of each size (up to max(m))
  max_m <- ifelse(length(m) == 0, 0, max(m))
  fam_num <- integer(max_m + 1)
  for (i in seq_len(max_m)) fam_num[i] <- sum(m == i)
  
  ## -------------- Estimation under H0 ----------------
  # initial r
  r <- 0; z <- 1; v2 <- 0
  mu0 <- mean(unlist(x), na.rm = TRUE); var0 <- var(unlist(x), na.rm = TRUE)
  if (is.na(var0) || var0 == 0) var0 <- 1e-8
  
  while (z > 1e-4 && v2 <= 30) {
    cvec <- numeric(nfam); hvec <- numeric(nfam)
    for (k in seq_len(nfam)) {
      mk <- m[k]
      cvec[k] <- (r * (2 - mk) - 1) / ((r - 1) * (mk * r + 1 - r))
      hvec[k] <- ifelse(mk > 1, r / ((r - 1) * (mk * r + 1 - r)), 0)
    }
    # mu0 update
    ns0 <- 0; ds0 <- 0
    for (k in seq_len(nfam)) {
      for (i in seq_len(m[k])) {
        ns0 <- ns0 + cvec[k] * x[[k]][i] + hvec[k] * (m[k] - 1) * x[[k]][i]
      }
      ds0 <- ds0 + m[k] * cvec[k] + hvec[k] * m[k] * (m[k] - 1)
    }
    mu0_new <- ifelse(ds0 != 0, ns0 / ds0, mu0)
    # var0 update
    nsv0 <- 0
    for (k in seq_len(nfam)) {
      for (i in seq_len(m[k])) {
        nsv0 <- nsv0 + cvec[k] * ((x[[k]][i] - mu0_new)^2)
        if (i < m[k]) {
          for (j in seq((i + 1), m[k])) {
            nsv0 <- nsv0 + 2 * hvec[k] * (x[[k]][i] - mu0_new) * (x[[k]][j] - mu0_new)
          }
        }
      }
    }
    var0_new <- ifelse(sum(m) != 0, nsv0 / sum(m), var0)
    # rh0 function
    rh0 <- function(rval) {
      s0 <- 0
      for (k in seq_len(nfam)) {
        mk <- m[k]
        if (mk < 2) next
        term1 <- ((mk) * rval * (mk - 1)) / (2 * (1 - rval) * (1 + rval * mk - rval))
        s0 <- s0 + term1
        for (i in seq_len(mk)) {
          xi <- x[[k]][i]
          s0 <- s0 - (1 / (2 * var0_new)) * ((rval * (mk - 1) * (2 - 2 * rval + mk * rval)) / (((mk - 1) * (rval^2) + rval * (2 - mk) - 1)^2)) * ((xi - mu0_new)^2)
          if (i < mk) {
            for (j in seq((i + 1), mk)) {
              xj <- x[[k]][j]
              if (mk > 1) {
                s0 <- s0 - (1 / var0_new) * ((-1 - (mk - 1) * rval^2) / (((mk - 1) * (rval^2) + rval * (2 - mk) - 1)^2)) * ((xi - mu0_new) * (xj - mu0_new))
              }
            }
          }
        }
      }
      return(s0)
    }
    # find r0
    r0 <- tryCatch(sl1(rh0, lower = -0.999, upper = 0.999, by = 0.01), error = function(e) 0)
    v2 <- v2 + 1
    z <- abs(r0 - r)
    if (z > 1e-4) r <- r0
    mu0 <- mu0_new; var0 <- ifelse(var0_new > 0, var0_new, var0)
    if (v2 > 30) { r <- 0; break }
  }
  r0 <- r
  
  # Likelihood under H0
  lhd0 <- 0
  for (k in seq_len(nfam)) {
    mk <- m[k]
    if (mk == 0) next
    m0_vec <- rep(mu0, mk)
    cov_mat <- var0 * ((1 - r0) * idm(mk) + r0 * jm(mk))
    lhd0 <- lhd0 + LaplacesDemon::dmvn(x[[k]], m0_vec, cov_mat, log = TRUE)
  }
  
  ## ----------------- Estimation under H1 ----------------
  r <- 0; z <- 1; v2 <- 0
  mu1 <- c(mu0, mu0)
  while (z > 1e-4 && v2 <= 30) {
    cvec <- numeric(nfam); hvec <- numeric(nfam)
    for (k in seq_len(nfam)) {
      mk <- m[k]
      cvec[k] <- (r * (2 - mk) - 1) / ((r - 1) * (mk * r + 1 - r))
      hvec[k] <- ifelse(mk > 1, r / ((r - 1) * (mk * r + 1 - r)), 0)
    }
    # define function for mu vector root (2 parameters)
    muh1 <- function(mu_vec) {
      mu_vec <- as.numeric(mu_vec)
      ns1 <- 0; ns2 <- 0
      for (k in seq_len(nfam)) {
        for (i in seq_len(m[k])) {
          xi <- x[[k]][i]
          yi <- y[[k]][i]
          resid <- xi - mu_vec[1] * yi + mu_vec[2] * yi - mu_vec[2]
          ns1 <- ns1 + (cvec[k]) * resid * (-yi)
          ns2 <- ns2 + (cvec[k]) * resid * (yi - 1)
          if (i < m[k]) {
            for (j in seq((i + 1), m[k])) {
              xj <- x[[k]][j]; yj <- y[[k]][j]
              resid_j <- xj - mu_vec[1] * yj + mu_vec[2] * yj - mu_vec[2]
              ns1 <- ns1 + hvec[k] * (resid * (-yj) + resid_j * (-yi))
              ns2 <- ns2 + hvec[k] * (resid * (yj - 1) + resid_j * (yi - 1))
            }
          }
        }
      }
      return(c(ns1, ns2))
    }
    # initial guess: two group means
    try_mu_init <- c(
      ifelse(sum(unlist(y)) > 0, mean(unlist(x)[unlist(y) == 1], na.rm = TRUE), mu0),
      ifelse(sum(1 - unlist(y)) > 0, mean(unlist(x)[unlist(y) == 0], na.rm = TRUE), mu0)
    )
    # solve for mu1 using Newton (slfn)
    mu1 <- tryCatch(slfn(muh1, try_mu_init, tol = 1e-6, max_iter = 100), error = function(e) try_mu_init)
    # variance under H1
    nsv1 <- 0
    for (k in seq_len(nfam)) {
      for (i in seq_len(m[k])) {
        resid_i <- x[[k]][i] - y[[k]][i] * mu1[1] + y[[k]][i] * mu1[2] - mu1[2]
        nsv1 <- nsv1 + cvec[k] * (resid_i^2)
        if (i < m[k]) {
          for (j in seq((i + 1), m[k])) {
            resid_j <- x[[k]][j] - y[[k]][j] * mu1[1] + y[[k]][j] * mu1[2] - mu1[2]
            nsv1 <- nsv1 + 2 * hvec[k] * (resid_i * resid_j)
          }
        }
      }
    }
    var1 <- ifelse(sum(m) != 0, nsv1 / sum(m), var0)
    # rh1 function (same structure, using mu1)
    rh1 <- function(rval) {
      s1 <- 0
      for (k in seq_len(nfam)) {
        mk <- m[k]
        if (mk < 2) next
        s1 <- s1 + ((mk) * rval * (mk - 1)) / (2 * (1 - rval) * (1 + rval * mk - rval))
        for (i in seq_len(mk)) {
          xi <- x[[k]][i]; yi <- y[[k]][i]
          resid <- xi - yi * mu1[1] + yi * mu1[2] - mu1[2]
          s1 <- s1 - (1 / (2 * var0)) * ((rval * (mk - 1) * (2 - 2 * rval + mk * rval)) / (((mk - 1) * (rval^2) + rval * (2 - mk) - 1)^2)) * (resid^2)
          if (i < mk) {
            for (j in seq((i + 1), mk)) {
              xj <- x[[k]][j]; yj <- y[[k]][j]
              resid_j <- xj - yj * mu1[1] + yj * mu1[2] - mu1[2]
              if (mk > 1) {
                s1 <- s1 - (1 / var0) * ((-1 - (mk - 1) * rval^2) / (((mk - 1) * (rval^2) + rval * (2 - mk) - 1)^2)) * (resid * resid_j)
              }
            }
          }
        }
      }
      return(s1)
    }
    r1 <- tryCatch(sl1(rh1, lower = -0.999, upper = 0.999, by = 0.01), error = function(e) 0)
    v2 <- v2 + 1
    z <- abs(r1 - r)
    if (z > 1e-4) r <- r1
    if (v2 > 30) { r <- 0; break }
  }
  r1 <- r
  
  # Likelihood under H1
  lhd1 <- 0
  for (k in seq_len(nfam)) {
    mk <- m[k]
    if (mk == 0) next
    m1_vec <- ifelse(y[[k]] == 1, mu1[1], mu1[2])
    cov_mat <- var1 * ((1 - r1) * idm(mk) + r1 * jm(mk))
    lhd1 <- lhd1 + LaplacesDemon::dmvn(x[[k]], m1_vec, cov_mat, log = TRUE)
  }
  
  lrt <- -2 * (lhd0 - lhd1)
  p_val <- pchisq(lrt, df = 1, lower.tail = FALSE)
  D <- ifelse(p_val <= 0.05, "Significant", "Insignificant")
  
  # ---------------- Non-parametric permutation test --------------
  st <- numeric(n_perm)
  # precompute counts used inside permutation
  # sample-by-family
  for (t in seq_len(n_perm)) {
    x_cpy <- lapply(x, function(v) if (length(v) > 1) sample(v) else v)
    # compute r_in_cpy
    rs0_cpy <- 0; n2_cpy <- 0
    for (i in seq(2, max(m))) {
      dm_cpy <- do.call(rbind, lapply(seq_len(nfam)[m == i], function(k) x_cpy[[k]]))
      n1_cpy <- ifelse(is.null(dim(dm_cpy)), 0, nrow(dm_cpy))
      if (n1_cpy > 0) rs0_cpy <- rs0_cpy + n1_cpy * incc(dm_cpy)
      n2_cpy <- n2_cpy + sum(m > 1)
    }
    r_in_cpy <- ifelse(n2_cpy != 0, rs0_cpy / n2_cpy, 0)
    cs1_cpy <- unlist(lapply(seq_len(nfam)[m == 1], function(k) x_cpy[[k]]))
    cs_cpy <- numeric(0); cn_cpy <- numeric(0)
    for (k in seq_len(nfam)) {
      if (m[k] > 1) {
        cs_k <- (t(y[[k]]) %*% (x_cpy[[k]])) / sum(y[[k]])
        cs_k <- cs_k / sqrt((1 / sum(y[[k]])) + r_in_cpy * (1 - 1 / sum(y[[k]])))
        cn_k <- (t(1 - y[[k]]) %*% (x_cpy[[k]])) / sum(1 - y[[k]])
        cn_k <- cn_k / sqrt((1 / sum(1 - y[[k]])) + r_in_cpy * (1 - 1 / sum(1 - y[[k]])))
        if (!is.na(cs_k)) cs_cpy <- c(cs_cpy, cs_k)
        if (!is.na(cn_k)) cn_cpy <- c(cn_cpy, cn_k)
      }
    }
    T1_cpy <- if (length(cs1_cpy) > 0 & length(cn_cpy) > 0) wilcox.test(cn_cpy, cs1_cpy, paired = FALSE)$statistic else 0
    T2_cpy <- if (length(cs_cpy) > 0 & length(cn_cpy) > 0) wilcox.test(cn_cpy, cs_cpy, paired = TRUE)$statistic else 0
    st[t] <- (as.numeric(T1_cpy) + as.numeric(T2_cpy))
  }
  
  # observed T
  # compute r_in for observed
  rs0 <- 0; n2 <- 0
  for (i in seq(2, max(m))) {
    dm <- do.call(rbind, lapply(seq_len(nfam)[m == i], function(k) x[[k]]))
    n1 <- ifelse(is.null(dim(dm)), 0, nrow(dm))
    if (n1 > 0) rs0 <- rs0 + n1 * incc(dm)
    n2 <- n2 + sum(m > 1)
  }
  r_in <- ifelse(n2 != 0, rs0 / n2, 0)
  cs1 <- unlist(lapply(seq_len(nfam)[m == 1], function(k) x[[k]]))
  cs <- numeric(0); cn <- numeric(0)
  for (k in seq_len(nfam)) {
    if (m[k] > 1) {
      cs_k <- (t(y[[k]]) %*% (x[[k]])) / sum(y[[k]])
      cs_k <- cs_k / sqrt((1 / sum(y[[k]])) + r_in * (1 - 1 / sum(y[[k]])))
      cn_k <- (t(1 - y[[k]]) %*% (x[[k]])) / sum(1 - y[[k]])
      cn_k <- cn_k / sqrt((1 / sum(1 - y[[k]])) + r_in * (1 - 1 / sum(1 - y[[k]])))
      if (!is.na(cs_k)) cs <- c(cs, cs_k)
      if (!is.na(cn_k)) cn <- c(cn, cn_k)
    }
  }
  T1 <- if (length(cs1) > 0 & length(cn) > 0) wilcox.test(cn, cs1, paired = FALSE)$statistic else 0
  T2 <- if (length(cs) > 0 & length(cn) > 0) wilcox.test(cn, cs, paired = TRUE)$statistic else 0
  T_obs <- as.numeric(T1) + as.numeric(T2)
  np <- sum(st > T_obs)
  p_np <- np / n_perm
  Dnp <- ifelse(p_np > 0.025 & p_np < 0.975, "Insignificant", "Significant")
  p_val_np <- 2 * min(c(p_np, 1-p_np))
  
  comp_sim <- data.frame(
    lhd0 = lhd0, lhd1 = lhd1, lrt = lrt, p_val = p_val, D = D,
    T = T_obs, p_np = p_val_np, Dnp = Dnp
  )
  return(comp_sim)
}

# Example usage:
# df <- read.csv("toydata.csv", header = TRUE)
result <- cmrbd(toydata, value_col = "Trait_val", affected_label = 1)
# print(result)
