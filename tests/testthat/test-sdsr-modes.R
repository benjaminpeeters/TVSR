# ============================================================================
# Tests for the four-spec SDSR introduced in TVSR 0.1.2.
#
# Fixture: a deterministic n=15, T=80 DGP with slope-up rho, slope-up trend
# and constant variance. Small enough to keep test wall-time < 30s.
#
# The golden RDS snapshot for the default (tv, tv) spec is committed at
# tests/testthat/sdsr_default_golden.rds. Regenerate via:
#     Rscript -e 'testthat::source_test_helpers("tests/testthat");
#                 saveRDS(make_fixture_fit(),
#                         "tests/testthat/sdsr_default_golden.rds")'
# ============================================================================

make_fixture <- function() {
  set.seed(101)
  # kLargestW requires the input matrix to carry dimnames (it assigns by
  # name inside apply(...)).
  N   <- 15L
  rn  <- paste0("r", seq_len(N))
  cn  <- paste0("c", seq_len(N))
  EXP <- matrix(rep(stats::rexp(N, rate = 1), N), N, N)
  W   <- matrix(stats::runif(N * N), N, N, dimnames = list(rn, cn))
  W   <- W * EXP * t(EXP)
  W   <- TVSR::kLargestW(W, 7L)
  W   <- TVSR::normalizationMatrix(W, type = "eigenvalues")

  Tt  <- 80L
  rho <- seq(0.2, 0.7, length.out = Tt)
  trd <- seq(1,   5,   length.out = Tt)
  vr  <- rep(2, Tt)

  set.seed(202)
  In    <- diag(N); one_n <- rep(1, N)
  Y     <- matrix(NA_real_, N, Tt)
  for (t in seq_len(Tt)) {
    eps    <- stats::rnorm(N, 0, sqrt(vr[t]))
    Y[, t] <- trd[t] * one_n + solve(In - rho[t] * W, eps)
  }
  list(Y = Y, W = W, Wlist = replicate(Tt, W, simplify = FALSE),
       T = Tt, rho_true = rho, trd_true = trd, var_true = vr)
}

fit_default_sdsr <- function(fx) {
  set.seed(303)
  TVSR::SDSR(fx$Y, fx$Wlist, verbose = FALSE, mc.cores = 1L)
}

test_that("default SDSR output matches the committed golden snapshot", {
  skip_on_cran()
  fx  <- make_fixture()
  got <- fit_default_sdsr(fx)

  golden_path <- test_path("sdsr_default_golden.rds")
  skip_if_not(file.exists(golden_path),
              "golden snapshot not present; regenerate via README")
  want <- readRDS(golden_path)

  expect_equal(got$RHO,      want$RHO,      tolerance = 0)
  expect_equal(got$TRD,      want$TRD,      tolerance = 0)
  expect_equal(got$VAR,      want$VAR,      tolerance = 0)
  expect_equal(got$omegaRho, want$omegaRho, tolerance = 0)
  expect_equal(got$aRho,     want$aRho,     tolerance = 0)
  expect_equal(got$bRho,     want$bRho,     tolerance = 0)
  expect_equal(got$f1Rho,    want$f1Rho,    tolerance = 0)
  expect_equal(got$lik,      want$lik,      tolerance = 0)
})

test_that("all four (trd, var) spec combinations return well-formed output", {
  skip_on_cran()
  fx <- make_fixture()
  for (spec in list(c("tv", "tv"), c("const", "tv"),
                    c("tv", "const"), c("const", "const"))) {
    set.seed(303)
    fit <- TVSR::SDSR(fx$Y, fx$Wlist, verbose = FALSE,
                      model = list(trd = spec[1], var = spec[2]),
                      mc.cores = 1L)
    expect_length(fit$RHO, fx$T)
    expect_length(fit$TRD, fx$T)
    expect_length(fit$VAR, fx$T)
    expect_true(all(is.finite(fit$RHO)))
    expect_true(all(is.finite(fit$TRD)))
    expect_true(all(is.finite(fit$VAR) & fit$VAR > 0))
    expect_true(is.finite(fit$lik))
  }
})

test_that("const specs return flat TRD / VAR paths", {
  skip_on_cran()
  fx <- make_fixture()

  set.seed(303)
  fit_cc <- TVSR::SDSR(fx$Y, fx$Wlist, verbose = FALSE, mc.cores = 1L,
                       model = list(trd = "const", var = "const"))
  expect_equal(diff(range(fit_cc$TRD)), 0, tolerance = 1e-12)
  expect_equal(diff(range(fit_cc$VAR)), 0, tolerance = 1e-12)

  set.seed(303)
  fit_ct <- TVSR::SDSR(fx$Y, fx$Wlist, verbose = FALSE, mc.cores = 1L,
                       model = list(trd = "const", var = "tv"))
  expect_equal(diff(range(fit_ct$TRD)), 0, tolerance = 1e-12)
  expect_gt(diff(range(fit_ct$VAR)), 0)

  set.seed(303)
  fit_tc <- TVSR::SDSR(fx$Y, fx$Wlist, verbose = FALSE, mc.cores = 1L,
                       model = list(trd = "tv", var = "const"))
  expect_gt(diff(range(fit_tc$TRD)), 0)
  expect_equal(diff(range(fit_tc$VAR)), 0, tolerance = 1e-12)
})

test_that("loglikTVRhoCond_mixed(tv, tv) is bit-identical to loglikTVRhoCond", {
  fx <- make_fixture()
  args <- list(omegaRho = 0.1, aRho = 0.05, bRho = 0.8, f1Rho = atanh(0.4))
  orig <- do.call(TVSR:::loglikTVRhoCond,
                  c(list(Y = fx$Y, w = fx$W), args))
  mixd <- do.call(TVSR:::loglikTVRhoCond_mixed,
                  c(list(Y = fx$Y, w = fx$W, trd_mode = "tv", var_mode = "tv"),
                    args))
  expect_identical(mixd, orig)
})

test_that("loglikTVRhoCond_mixed_tvW(tv, tv) is bit-identical to loglikTVRhoCondtvW", {
  fx <- make_fixture()
  args <- list(omegaRho = 0.1, aRho = 0.05, bRho = 0.8, f1Rho = atanh(0.4))
  orig <- do.call(TVSR:::loglikTVRhoCondtvW,
                  c(list(Y = fx$Y, W = fx$Wlist), args))
  mixd <- do.call(TVSR:::loglikTVRhoCond_mixed_tvW,
                  c(list(Y = fx$Y, W = fx$Wlist, trd_mode = "tv", var_mode = "tv"),
                    args))
  expect_identical(mixd, orig)
})

test_that("SDSR rejects model$trd = \"gas\" (reserved for future release)", {
  fx <- make_fixture()
  expect_error(
    TVSR::SDSR(fx$Y, fx$Wlist, verbose = FALSE, mc.cores = 1L,
               model = list(trd = "gas", var = "tv")),
    "reserved for a future release"
  )
})
