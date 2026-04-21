# TVSR 0.1.1

## Backwards-compatible additions

- `SDSR()`, `LKSR()`, `SRstatic()`, and `crossValidLKSR()` now accept an
  `mc.cores` argument controlling the number of forked workers used in their
  internal grid searches. Default: `max(1L, parallel::detectCores() - 2L)`,
  which leaves 2 physical cores free for the rest of the workstation.
- Previously these functions hardcoded `parallel::detectCores()` internally.
  The new behaviour is a strict superset: callers that do not pass
  `mc.cores` still get a sensible default, and numerical outputs are
  unchanged (grid evaluations are deterministic, so worker count affects
  only wall-clock time). Verified by the numerical-equivalence test in
  `p8_spatial_regression_bias/tvsr_numerical_equivalence_tests/`.
- `LKSR()` and `crossValidLKSR()` now pass `mc.cores = 1L` to their internal
  `SRstatic()` calls to prevent nested parallelism when the outer loop is
  already forked.
