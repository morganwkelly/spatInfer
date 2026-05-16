library(spatInfer)

data(opportunity)

# Small reproducible subset to keep tests fast
set.seed(42)
opp <- opportunity |> dplyr::slice_sample(n = 150)

# Shared parameters derived from a quick optimal_basis run
SPLINES  <- 4
PC_NUM   <- 3
CLUSTERS <- 4

test_that("basis_regression runs without weights (BCH)", {
  result <- basis_regression(
    mobility ~ single_mothers + social_cap,
    opp,
    splines  = SPLINES,
    pc_num   = PC_NUM,
    clusters = CLUSTERS,
    cov      = "BCH"
  )
  expect_s3_class(result, "fixest")
  # Treatment variable should be present
  expect_true("single_mothers" %in% names(coef(result)))
})

test_that("basis_regression runs without weights (HC)", {
  result <- basis_regression(
    mobility ~ single_mothers + social_cap,
    opp,
    splines  = SPLINES,
    pc_num   = PC_NUM,
    clusters = CLUSTERS,
    cov      = "HC"
  )
  expect_s3_class(result, "fixest")
  expect_true("single_mothers" %in% names(coef(result)))
})

test_that("basis_regression runs with weights (BCH)", {
  set.seed(7)
  opp_w <- opp
  opp_w$wts <- runif(nrow(opp_w), 0.5, 2.0)

  result <- basis_regression(
    mobility ~ single_mothers + social_cap,
    opp_w,
    splines  = SPLINES,
    pc_num   = PC_NUM,
    clusters = CLUSTERS,
    weights  = TRUE,
    cov      = "BCH"
  )
  expect_s3_class(result, "fixest")
  expect_true("single_mothers" %in% names(coef(result)))
})

test_that("basis_regression runs with weights (HC)", {
  set.seed(7)
  opp_w <- opp
  opp_w$wts <- runif(nrow(opp_w), 0.5, 2.0)

  result <- basis_regression(
    mobility ~ single_mothers + social_cap,
    opp_w,
    splines  = SPLINES,
    pc_num   = PC_NUM,
    clusters = CLUSTERS,
    weights  = TRUE,
    cov      = "HC"
  )
  expect_s3_class(result, "fixest")
  expect_true("single_mothers" %in% names(coef(result)))
})

test_that("weighted and unweighted estimates differ", {
  set.seed(7)
  opp_w <- opp
  opp_w$wts <- runif(nrow(opp_w), 0.5, 2.0)

  unweighted <- basis_regression(
    mobility ~ single_mothers + social_cap,
    opp_w,
    splines  = SPLINES,
    pc_num   = PC_NUM,
    clusters = CLUSTERS,
    weights  = FALSE,
    cov      = "HC"
  )
  weighted <- basis_regression(
    mobility ~ single_mothers + social_cap,
    opp_w,
    splines  = SPLINES,
    pc_num   = PC_NUM,
    clusters = CLUSTERS,
    weights  = TRUE,
    cov      = "HC"
  )
  expect_false(
    isTRUE(all.equal(coef(unweighted)[["single_mothers"]],
                     coef(weighted)[["single_mothers"]]))
  )
})

test_that("basis_regression stops with missing X/Y", {
  opp_bad <- opp
  opp_bad$X[1] <- NA
  expect_error(
    basis_regression(
      mobility ~ single_mothers,
      opp_bad,
      splines  = SPLINES,
      pc_num   = PC_NUM,
      clusters = CLUSTERS
    ),
    "missing values in longitude"
  )
})
