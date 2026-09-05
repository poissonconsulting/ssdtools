test_that("vectorised quantile dispatch matches the elementwise path", {
  p <- c(NA, NaN, 0, 1e-6, 0.05, 0.5, 0.95, 1, 2, -1)
  for (dist in ssd_dists_all()) {
    qfun <- get(paste0("ssd_q", dist))
    pars <- get(paste0("ssd_e", dist))()
    vectorised <- do.call(qfun, c(list(p), pars))
    elementwise <- vapply(p, function(pi) do.call(qfun, c(list(pi), pars)), 1)
    expect_identical(vectorised, elementwise, info = dist)
  }
})

test_that("vectorised cumulative dispatch matches the elementwise path", {
  q <- c(NA, NaN, -Inf, -1, 0, 1e-6, 1, 10, Inf)
  for (dist in ssd_dists_all()) {
    pfun <- get(paste0("ssd_p", dist))
    pars <- get(paste0("ssd_e", dist))()
    vectorised <- do.call(pfun, c(list(q), pars))
    elementwise <- vapply(q, function(qi) do.call(pfun, c(list(qi), pars)), 1)
    expect_identical(vectorised, elementwise, info = dist)
  }
})

test_that("invalid scalar parameters give NaN for every element", {
  expect_identical(ssd_qlnorm(c(0.1, 0.5, NA), sdlog = -1), c(NaN, NaN, NA))
  expect_identical(ssd_plnorm(c(1, 2, NaN), sdlog = 0), c(NaN, NaN, NaN))
  expect_identical(ssd_qlnorm_lnorm(c(0.1, 0.5), pmix = 1), c(NaN, NaN))
})

test_that("parameter vectors still recycle elementwise", {
  expect_identical(
    ssd_qlnorm(c(0.25, 0.75), meanlog = c(0, 1), sdlog = c(1, 2)),
    c(ssd_qlnorm(0.25, 0, 1), ssd_qlnorm(0.75, 1, 2))
  )
  expect_identical(
    ssd_qlnorm(c(0.25, 0.75), meanlog = c(0, NA), sdlog = 1),
    c(ssd_qlnorm(0.25, 0, 1), NA_real_)
  )
})
