test_that("multivariate sd two groups", {
  X <- data.frame(
    x1 = c(0.2, 0.6, 0.9, 0.8),
    x2 = c(0.6, 0.2, 0.8, 0.9)
  )
  y <- c(0.0, 1.0, 0.5, 0.6)
  w <- rep(1, length(y))

  groups <- stats::setNames(c(1, 1), c("x1", "x2"))
  orders <- c("sd" = 1)

  pars <- list(
    verbose = FALSE,
    eps_abs = 1e-6,
    eps_rel = 1e-6,
    max_iter = 10000L
  )

  fit <- idr(
    y = y,
    X = X,
    weights = w,
    groups = groups,
    orders = orders,
    stoch = "sd",
    pars = pars,
    progress = FALSE
  )

  preds <- predict(fit, X)

  # 4 predictions (one per row), 4 thresholds
  K <- 4L
  expect_equal(nrow(preds$cdf), 4L)
  expect_equal(ncol(preds$cdf), K)

  # Rows 1,2 collapse under SD (same sorted pair), as do rows 3,4
  expect_equal(preds$cdf[1, ], preds$cdf[2, ])
  expect_equal(preds$cdf[3, ], preds$cdf[4, ])

  # Extract unique CDFs (one per collapsed group)
  unique_cdfs <- preds$cdf[c(1, 3), ]
  N <- 2L

  for (i in seq_len(N)) {
    expect_true(all(diff(unique_cdfs[i, ]) >= 0))
    expect_true(all(unique_cdfs[i, ] >= 0 & unique_cdfs[i, ] <= 1))
  }
  expect_equal(unique_cdfs[, K], rep(1, N), tolerance = 1e-12)

  expect_equal(
    unique_cdfs,
    matrix(c(0.5, 0, 0.5, 0.5, 0.75, 0.75, 1.0, 1.0), nrow = 2),
    tolerance = 1e-6
  )
})

test_that("multivariate sd case larger", {
  # 6 observations that collapse to 3 unique rows under SD preparation
  X <- data.frame(
    x1 = c(0.2, 0.3, 0.4, 0.6, 0.8, 0.9),
    x2 = c(0.3, 0.2, 0.6, 0.4, 0.9, 0.8)
  )
  y <- c(0.1, 0.2, 0.2, 0.6, 0.7, 0.9)
  w <- rep(1, length(y))

  # One SD-ordered group
  groups <- stats::setNames(c(1, 1), c("x1", "x2"))
  orders <- c("sd" = 1)

  pars <- list(
    verbose = FALSE, eps_abs = 1e-8, eps_rel = 1e-8, max_iter = 10000L
  )

  fit <- idr(
    y = y,
    X = X,
    weights = w,
    groups = groups,
    orders = orders,
    stoch = "sd",
    pars = pars,
    progress = FALSE
  )

  preds <- predict(fit, X)

  # Thresholds are the sorted unique y values
  expect_equal(sort(unique(y)), preds$points)

  K <- length(preds$points)
  expect_equal(K, 5L)

  # 6 predictions; pairs (1,2), (3,4), (5,6) collapse under SD
  expect_equal(nrow(preds$cdf), 6L)
  expect_equal(ncol(preds$cdf), K)
  expect_equal(preds$cdf[1, ], preds$cdf[2, ])
  expect_equal(preds$cdf[3, ], preds$cdf[4, ])
  expect_equal(preds$cdf[5, ], preds$cdf[6, ])

  # Extract unique CDFs and sort by row mean (descending) for comparison
  unique_cdfs <- preds$cdf[c(1, 3, 5), ]
  K_inner <- K - 1L
  expect_equal(K_inner, 4L)

  row_order <- order(rowMeans(unique_cdfs[, 1:K_inner, drop = FALSE]),
    decreasing = TRUE
  )
  cdf_sorted <- unique_cdfs[row_order, 1:K_inner, drop = FALSE]

  expected <- rbind(
    c(0.5, 1.0, 1.0, 1.0), # A
    c(0.0, 0.5, 1.0, 1.0), # B
    c(0.0, 0.0, 0.0, 0.5) # C
  )

  expect_equal(cdf_sorted, expected, tolerance = 1e-7)
  expect_true(all(unique_cdfs >= 0 & unique_cdfs <= 1))
})

test_that("multivariate multi-ordering case", {
  # 6 observations that collapse to 3 unique rows under the mixed group orders
  # (SD, ICX, COMP)
  X <- data.frame(
    # G1 (SD): two exchangeable columns
    a1 = c(0.2, 0.1, 0.6, 0.4, 0.9, 0.8),
    a2 = c(0.1, 0.2, 0.4, 0.6, 0.8, 0.9),
    # G2 (ICX): two exchangeable columns
    b1 = c(0.3, 0.2, 0.5, 0.3, 0.8, 0.7),
    b2 = c(0.2, 0.3, 0.3, 0.5, 0.7, 0.8),
    # G3 (COMP): one component-wise column
    c  = c(0.2, 0.2, 0.5, 0.5, 0.9, 0.9)
  )

  y <- c(0.1, 0.2, 0.2, 0.7, 0.6, 0.9)
  w <- rep(1, length(y))

  # Groups and orders: G1=1 (sd), G2=2 (icx), G3=3 (comp)
  groups <- stats::setNames(c(1, 1, 2, 2, 3), colnames(X))
  orders <- c("sd" = 1, "icx" = 2, "comp" = 3)

  pars <- list(
    verbose = FALSE,
    eps_abs = 1e-8,
    eps_rel = 1e-8,
    max_iter = 10000L
  )

  fit <- idr(
    y = y,
    X = X,
    weights = w,
    groups = groups,
    orders = orders,
    stoch = "sd",
    pars = pars,
    progress = FALSE
  )

  preds <- predict(fit, X)

  # Thresholds are sorted unique y
  expect_equal(preds$points, sort(unique(y)))

  K <- length(preds$points)

  # 6 predictions; pairs (1,2), (3,4), (5,6) collapse
  expect_equal(nrow(preds$cdf), 6L)
  expect_equal(ncol(preds$cdf), K)
  expect_equal(preds$cdf[1, ], preds$cdf[2, ])
  expect_equal(preds$cdf[3, ], preds$cdf[4, ])
  expect_equal(preds$cdf[5, ], preds$cdf[6, ])

  # Extract unique CDFs
  unique_cdfs <- preds$cdf[c(1, 3, 5), ]

  K_inner <- K - 1L
  expect_equal(K_inner, 4L)

  # Align rows by decreasing mean over the first K_inner columns
  row_order <- order(rowMeans(unique_cdfs[, 1:K_inner, drop = FALSE]),
    decreasing = TRUE
  )
  cdf_sorted <- unique_cdfs[row_order, 1:K_inner, drop = FALSE]

  expected <- rbind(
    c(0.5, 1.0, 1.0, 1.0), # R1
    c(0.0, 0.5, 0.5, 1.0), # R2
    c(0.0, 0.0, 0.5, 0.5) # R3
  )

  expect_equal(cdf_sorted, expected, tolerance = 1e-7)

  # Additional spot checks at specific thresholds:
  cdf_aligned <- unique_cdfs[row_order, ]
  # z = 0.2 -> [1.0, 0.5, 0.0] in row order (R1, R2, R3)
  z_idx <- match(0.2, preds$points)
  expect_gte(z_idx, 1L)
  expect_lte(z_idx, K_inner)
  vals_at_02 <- cdf_aligned[, z_idx]
  expect_equal(vals_at_02, c(1.0, 0.5, 0.0), tolerance = 1e-7)

  # z = 0.7 -> [1.0, 1.0, 0.5]
  z_idx <- match(0.7, preds$points)
  expect_gte(z_idx, 1L)
  expect_lte(z_idx, K_inner)
  vals_at_07 <- cdf_aligned[, z_idx]
  expect_equal(vals_at_07, c(1.0, 1.0, 0.5), tolerance = 1e-7)

  # Boundedness and monotonicity in thresholds
  expect_true(all(preds$cdf >= 0.0 & preds$cdf <= 1.0))
  for (i in seq_len(nrow(preds$cdf))) {
    expect_true(all(diff(preds$cdf[i, ]) >= 0.0))
  }
})

test_that("multivariate SD edge cases are valid", {
  # Columns: SD0, SD1, COMP0, ICX0, ICX1, COMP1
  sd_pairs <- list(
    c(0.2, 0.6), c(0.6, 0.2),
    c(0.9, 0.8), c(0.8, 0.9),
    c(0.5, 0.5), c(0.4, 0.1)
  )
  comp_pairs <- list(c(0.1, 0.9), c(0.7, 0.3), c(0.4, 0.5))
  icx_pairs <- list(c(0.8, 0.3), c(0.3, 0.8), c(0.9, 0.1), c(0.1, 0.9))

  X <- matrix(nrow = 0, ncol = 6)
  for (t in seq_along(sd_pairs)) {
    sp <- sd_pairs[[t]]
    cp <- comp_pairs[[(t - 1) %% length(comp_pairs) + 1]]
    ip <- icx_pairs[[(t - 1) %% length(icx_pairs) + 1]]
    X <- rbind(X, c(sp[1], sp[2], cp[1], ip[1], ip[2], cp[2]))
    if ((t - 1) %% 2 == 0) {
      # Add a "paired" line that will collapse under SD/ICX preparation
      X <- rbind(X, c(sp[2], sp[1], cp[1], ip[2], ip[1], cp[2]))
    }
  }
  X <- as.data.frame(X)
  names(X) <- c("SD0", "SD1", "C0", "I0", "I1", "C1")

  n <- nrow(X)
  levels <- c(0.0, 0.3, 0.7, 1.0)
  y <- levels[((seq_len(n) - 1) %% length(levels)) + 1L]
  w <- rep(1.0, n)
  if (n >= 3) w[3] <- 2.0
  if (n >= 5) w[5] <- 0.5
  w[n] <- 1e-4

  # Groups: 1-based labels for R
  # group 1: SD0, SD1 -> "sd"
  # group 2: C0,  C1  -> "comp"
  # group 3: I0,  I1  -> "icx"
  groups <- setNames(c(1, 1, 2, 3, 3, 2), names(X))
  orders <- c("sd" = 1, "comp" = 2, "icx" = 3)

  # Fit via R implementation (multivariate + stoch == "sd")
  pars <- list(
    verbose = FALSE,
    eps_abs = 1e-6,
    eps_rel = 1e-6,
    max_iter = 20000L
  )
  fit <- idr(
    y = y,
    X = X,
    weights = w,
    groups = groups,
    orders = orders,
    stoch = "sd",
    pars = pars,
    progress = FALSE
  )

  preds <- predict(fit, X)

  # Basic validity checks on predicted CDFs
  expect_true(is.matrix(preds$cdf))
  expect_equal(ncol(preds$cdf), length(unique(y)))
  # last column exactly 1
  expect_true(all(abs(preds$cdf[, ncol(preds$cdf)] - 1) < 1e-12))
  # in [0,1]
  expect_true(all(preds$cdf >= 0 & preds$cdf <= 1))
  # nondecreasing in thresholds
  diffs <- t(apply(preds$cdf, 1, diff))
  expect_true(all(diffs >= 0))
})
