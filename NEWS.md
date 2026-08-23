# bilatr (development version)

## Development notes (internal, not part of the public API)

* Added an experimental Stan model variant, `inst/stan/bilatr_phi_logn.stan`,
  in which the Dirichlet-multinomial concentration `phi` is modeled as a
  per-dyad-period quantity (`log(phi[d,t]) = log_phi0[d] + beta_logn *
  centered_log_n[d,t]`) rather than a per-dyad constant. Shares the stable
  model's data interface exactly (no changes to `assemble_stan_data()` were
  needed) and nests the stable model's dispersion structure as
  `beta_logn -> 0`.
* Added an internal Stan model registry (`R/model-registry.R`,
  `.bilatr_stan_models`) mapping short names to `inst/stan/` files, plus a
  shared, cached compilation path (`.compile_stan_model()`) used by both
  the exported fitters and the new internal dev entry points below.
  `fit_dyad_ts()` and `fit_panel()` are unaffected: their signatures and
  behavior are unchanged, and they always resolve the registry's
  `"stable"` entry internally.
* Added internal (non-exported) development entry points
  `fit_dyad_ts_dev()`/`fit_panel_dev()` in `R/fit-dev.R`, identical to
  `fit_dyad_ts()`/`fit_panel()` but with a `stan_model` argument for
  selecting a registered model during active development. Not part of the
  public API.
