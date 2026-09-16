
# bilatr 0.5.2

## Changes

* `diagnose_and_extract_bilatr()` now returns `theta_filtered`/
  `theta_filtered_sd` (empty tibbles unless `compute_theta_filtered = 1`
  was set when the fit was assembled) alongside `theta`/`alpha`/
  `mu_intercept`. These columns were already being read structurally as
  part of the existing Tier 3 sweep (same bracket-index-count rule as
  `theta`/`log_lik`; see 0.5.0's NEWS entry), but were silently dropped
  by the extraction step, which was hardcoded to only pull out
  `"theta["`-prefixed rows -- no second CSV read is needed for this fix.
  Correctly translates the position-within-`filter_dyads` index these two
  variables carry back to the true `dyad_id` via `stan_data$filter_dyads`
  before joining `dyad_ids`, resolving the documented limitation from
  0.5.0's own NEWS entry (previously only valid when every dyad was
  filtered; now correct for an explicit subset too).

# bilatr 0.5.1

## New features

* Added `EventRootCode4`, a finer alternative to `EventRootCode3`:
  `041` ("Discuss by telephone") is split out of the "Consult: meet,
  discuss, or visit" grouping into its own class, and the "Threaten or
  exhibit force posture" and "Assault, fight, or mass violence" groupings
  are split back into their individual root codes (Threaten, Exhibit
  force posture; Assault, Fight, Use unconventional mass violence).
  `016`/`018`/`019`'s relocations into Reject/diplomatic-cooperation
  carry over unchanged from `EventRootCode3`. Classes are numbered
  1-23. `cameo_lookup`/`recode_cameo()` gain `EventRootCode4`,
  `EventRootCode4Name`, and `EventRootCode4RootCodes` alongside the new
  `assign_eventrootcode4()`, `eventrootcode4_name()`, and
  `eventrootcode4_rootcodes()` helpers.

# bilatr 0.5.0

## New features

* **`diagnose_category_merges()`**: for a fitted `alpha`, reports what
  merging any two (or more) of the `A` event categories would cost in
  Fisher information about `theta`. Two categories are redundant for the
  latent scale exactly when their discriminations are equal; the
  information lost by merging `j`/`k` is exactly
  `pi_j*pi_k*(alpha_j-alpha_k)^2/(pi_j+pi_k)` (Ward's linkage on `alpha`
  weighted by category share `pi`), computed per posterior draw and
  summarised (never at the posterior mean, since the loss is nonlinear
  and `alpha`'s sum-to-zero/RMS-1 constraints induce strong negative
  correlation between categories). Returns a `categories` table (with
  contraction/z-score secondary columns, computed against `alpha`'s
  actual per-index prior moments -- see `alpha_prior_moments()` below,
  not a flat `N(0,1)`), a `pairwise` table of every merge's cost, and a
  `ladder`: the full greedy agglomerative merge path from `A` down to 2,
  so the whole cost curve can be read at once rather than pricing one
  grouping at a time. `merge_cost()` prices an explicit grouping already
  in mind. Documents plainly that `alpha` is a discrimination, not a
  severity ranking, and that cheapness is a price, not a recommendation
  -- the function never suggests a grouping.
* **`alpha_prior_moments()`**: `alpha`'s prior is not exchangeable across
  action classes even though the underlying construction is -- folding by
  the sign of the reference class's own raw coordinate
  (`orientation_sign()`) breaks that symmetry, giving the reference class
  a materially different prior mean/sd than the rest (e.g. at `A = 10`,
  mean 0.820/sd 0.572 vs. mean -0.091/sd 0.996). Returns the per-index
  prior mean/sd by Monte Carlo simulation mirroring the Stan construction
  exactly, memoised per call signature.
* **Forward-filtered `theta`**: `theta_filtered`/`theta_filtered_sd` in
  `generated quantities` (both `stable` and `ou`), a West-Harrison/
  Fisher-scoring linear-Bayes filter -- conditional on each draw's
  hyperparameters, `p(theta_t | y_1:t)` rather than the smoothed `theta`
  every fit already returns (an approximation, since the hyperparameters
  were themselves fit on all `T` periods). Hand-derived score (no
  autodiff in `generated quantities`): with the likelihood unweighted
  since 0.4.6, `conc_0 = phi` is exactly constant in `theta`, so the DM
  score reduces to a clean sum with no special cases. Gated behind a new
  `compute_theta_filtered` data flag (default off: `D*T` iterations with
  `A` `digamma()` calls each, single-threaded, on the order of a minute
  or two per chain for a full production dyad set), with an optional
  `filter_dyads` argument on `assemble_stan_data()` (dyad names, matched
  against the `dyad_ids` attribute) to filter a cheap subset instead of
  every dyad.
* **`prior_only`** data flag (both `stable` and `ou`, plumbed through
  `assemble_stan_data()`): skips the `reduce_sum` likelihood call
  entirely, fitting the prior alone -- useful for prior-predictive checks
  and validates `alpha_prior_moments()` against an actual fit (tested).

## Changes

* `.bilatr_tier1_names` drops `mu_intercept_raw`/`mu_theta0`: neither
  exists in any currently registered model (dead names that could only
  ever silently match nothing). `theta_filtered`/`theta_filtered_sd` need
  no entry anywhere in `.classify_bilatr_tier()`: Tier 3 is matched
  structurally by bracket-index count, so they land there automatically,
  the same way `log_lik` already does. One documented limitation: when
  `filter_dyads` narrows the dyad set, these two variables' Tier 3
  `dyad_id` is the position within the filtered subset, not the true
  `D`-space `dyad_id` -- there is no translation layer back, out of scope
  for this release.
* `tests/testthat/helper_fixtures.R`'s `make_csv_diagnostics_fixture()`
  no longer builds the retired `dyad_weight`/`period_weight`/
  `action_weight` fields into its default `data_list` unconditionally
  (leftover from 0.4.6) -- added only when fitting one of the legacy
  `stable_soft_anchor`/`ou_soft_anchor` programs, which still declare
  them.

# bilatr 0.4.6

## Breaking changes

* **Retires the likelihood-weighting scheme** (`dyad_weight`,
  `period_weight`, `action_weight`) from the `stable`/`ou` Stan programs
  and from `assemble_stan_data()`. These were never used in production
  (`weighted = FALSE` throughout every runscript), and two of the three
  were pure likelihood tempering with no generative interpretation;
  `action_weight` was different in kind -- it multiplied inside the
  Dirichlet-multinomial concentration, making the total concentration
  depend on `theta` and complicating a forward-filter feature planned on
  top of this likelihood, which hand-differentiates the likelihood in
  `generated quantities` (no autodiff there).
* `partial_log_lik()`/`dyad_period_log_lik()` no longer take
  `dyad_weight`/`period_weight`/`action_weight`; the concentration is now
  simply `phi[d] * p`. `compute_default_weights()` and
  `parse_weighted_arg()` are removed. `assemble_stan_data()` no longer
  returns `dyad_weight`/`period_weight`/`action_weight`, and its
  `weighted` argument is now defunct: it must be `FALSE` or `"none"`
  (both accepted as no-ops, matching this project's own SLURM runscripts,
  which already map their `weighted = "none"` CLI argument to `FALSE`),
  and errors with a message pointing here for any other value.
* Verified this is a pure deletion, not a behavior change: at unit
  weights (the only weights ever used), fitting the pre-0.4.6 `stable`
  program and the new one on identical data/inits/seed produces
  **bit-identical** `lp__`, draws, and sampler diagnostics
  (`identical()`, not merely equal to tolerance).
* Older `stan_data.rds` files (still carrying the three weight fields)
  remain usable with the new programs -- CmdStan ignores data fields a
  program doesn't declare. The legacy `stable_soft_anchor`/
  `ou_soft_anchor` programs (and the two pre-0.4.0 retirees) are
  unchanged and still declare/apply all three, so pre-0.4.6 output made
  under them stays readable; they are reachable only via the dev-only
  `fit_dyad_ts_dev()`/`fit_panel_dev()`, never via `fit_dyad_ts()`/
  `fit_panel()`.

# bilatr 0.4.5

## Changes

* `EventRootCode3`'s class 7 ("Engage in diplomatic cooperation") now
  also picks up `019` ("Express accord"), alongside `018`. Previously
  `019` stayed grouped under root 01 (class 1, "Make a public
  statement") along with the rest of root 01's codes; it now joins `018`
  in relocating out of that root, since both behave more like
  diplomatic cooperation than a neutral statement. `assign_eventrootcode3()`
  and `eventrootcode3_rootcodes()` (class 7's root-codes string is now
  `"05, 018, 019"`) are updated accordingly; `EventRootCode2`/
  `BilatrClass` are unaffected.

# bilatr 0.4.4

## New features

* Added `EventRootCode3`, a finer alternative to `EventRootCode2`:
  `016`/`018` are relocated into the Reject/diplomatic-cooperation
  classes (rather than staying grouped under root 01), and the
  `"10"`/`"13"` groupings are split back into Investigate + Demand and
  Threaten-or-exhibit-force-posture + Protest. Classes are numbered
  1-19. `cameo_lookup`/`recode_cameo()` gain `EventRootCode3`,
  `EventRootCode3Name`, and `EventRootCode3RootCodes` (the underlying
  CAMEO root/event codes each class draws from) alongside the new
  `assign_eventrootcode3()`, `eventrootcode3_name()`, and
  `eventrootcode3_rootcodes()` helpers.

# bilatr 0.4.3

## Breaking changes

* **Reverts 0.4.2's `alpha[1] > 0` hard constraint** and replaces it with
  an orientation fold. A production cluster run showed the 0.4.2
  constraint doesn't work: a hard constraint on `alpha_raw_1` removes one
  of the two reflection-symmetric modes from the parameter space, but not
  the likelihood barrier between them, so a chain that would have started
  (or drifted, during warmup) into the excluded orientation cannot cross
  to the feasible one -- it slides to the constraint boundary and parks
  there instead (observed: `alpha[1] = 0.00046`, 5-95% interval
  `[0.000025, 0.00127]`, wandering on the log scale, `lp__` materially
  worse than a healthy chain; confirmed to be the mirrored solution
  projected onto `alpha[1] = 0` and re-fitted around that boundary, not a
  repairable mirror image). No hard constraint, on any element, can fix
  this -- a constraint identifies by removing ambiguity, it cannot move a
  chain across a barrier.
* The fix changes what is *reported*, not what is *reachable*: `alpha_raw`
  is a free `sum_to_zero_vector[A]` again (as before 0.4.2), and both
  `alpha` and `theta` are multiplied by `sign(alpha_raw[1])` in
  `transformed parameters`. Since the likelihood depends on `alpha`/`theta`
  only through their elementwise product, this leaves the target
  completely unchanged (no Jacobian adjustment -- it is not a change of
  variables, just a canonical relabelling of the output) while making
  `alpha[1] >= 0` always hold in every REPORTED fit, regardless of which
  raw-space basin a chain occupies. `anchor_scale`/the soft anchor stay
  removed, as in 0.4.2.
* The RAW parameters the fold consumes (`alpha_raw`, `z_theta0`/
  `mu_dyad_raw`, `theta_raw`) remain genuinely sign-ambiguous themselves
  -- if two chains land in opposite raw-space basins, those variables'
  own cross-chain Rhat is meaningless even though everything reported is
  fine. `diagnose_convergence()`/`diagnose_and_extract_bilatr()` now
  exclude them from the tiered diagnostics tables entirely (new
  `.bilatr_sign_ambiguous_raw_names()`, derived from the same source of
  truth as `.bilatr_flip_variables()` so the two lists cannot drift
  apart) rather than let a meaningless Rhat surface as a false Tier 1
  alarm.
* Verified with a two-chain test giving deliberately opposite `alpha_raw`
  inits: both chains report `alpha[1] > 0` and agree on `alpha`/`theta`
  within MCMC error (small Rhat), while `alpha_raw`/`theta_raw` show large
  Rhat -- confirming the fold, not coincidence, is doing the work. Also
  re-confirmed (same-data comparison, `stable_soft_anchor` vs. `stable`)
  that this is a re-parameterization, not a model change: `alpha`/`theta`
  agree after orienting the legacy fit, `mu_intercept`/`phi`/`sigma_*`
  agree directly.
* 0.4.2's own "Breaking changes"/"Changes" entries below otherwise stand:
  the orientation-retirement (`bilatr_orient()` deprecation), benchmark
  Pss correction, and worker/chunk trade-off reporting are unaffected by
  this reversion.

# bilatr 0.4.2

## Breaking changes

* `stable`/`ou` now identify `alpha[1]`'s sign BY CONSTRUCTION instead of
  with a soft anchor. Previously, `alpha_raw` was a free
  `sum_to_zero_vector[A]` with a soft penalty
  (`target += log_inv_logit(alpha[1] * inv(anchor_scale))`) nudging
  `alpha[1]` positive; that penalty made the target correctly specified
  but could not move a chain across the alpha/theta reflection
  symmetry's likelihood barrier (thousands of nats) once warmup had
  landed it in a basin, so independently-initialized chains could
  disagree on sign with no diagnostic recourse short of post-hoc
  relabeling. `alpha_raw` is now built from `real<lower=0> alpha_raw_1`
  (the reference/neutral class) plus free middle elements
  (`vector[A - 2] alpha_raw_mid`), making `alpha[1] > 0` structural: no
  reflection symmetry remains, cross-chain Rhat on `alpha`/`theta` is
  meaningful again, and no post-hoc relabeling is needed. `anchor_scale`
  is no longer declared by (or consumed by) `stable`/`ou`, though
  `assemble_stan_data()` keeps supplying it unconditionally (CmdStan
  ignores data a program doesn't declare).
* The pre-0.4.2 `stable`/`ou` programs (free `sum_to_zero_vector`
  `alpha_raw`, soft anchor) are retired to `inst/stan/legacy/` and
  registered as `stable_soft_anchor`/`ou_soft_anchor` (`status =
  "legacy"`), kept only so CmdStan output produced before 0.4.2 stays
  readable via `bilatr_orient()`. Fits made under the new `stable`/`ou`
  are **not parameter-comparable** to fits made under the retired
  programs -- this is a re-parameterization, not merely a relabeling
  (verified to agree on `alpha`/`theta`/`mu_intercept` after orienting
  the legacy fit, on a same-data comparison; see `dev/
  claude_code_prompt_0.4.2_identification.md`). The pre-0.4.0
  `"alphanorm"`/`"alphanorm_ou"` aliases now resolve to these legacy
  entries (previously `"stable"`/`"ou"`), since every fit ever made
  under those alias names necessarily predates this change.

## Changes

* `bilatr_orient()` is deprecated: retained only for reading/re-deriving
  CmdStan output from the retired `stable_soft_anchor`/`ou_soft_anchor`
  programs, slated for removal once those runs are gone.
  `.bilatr_flip_variables()` returns `character(0)` for `stable`/`ou`
  accordingly, and every call site (`extract_theta()`, `extract_alpha()`,
  `extract_mu_intercept()`, `.warn_if_wrong_basin()`) now skips
  `bilatr_orient()`/the `alpha[1]` read it needs entirely for those two
  models, rather than calling it as a no-op.
* Corrected the CSV-path memory model's benchmark methodology
  (`dev/bench_memory.R`): it previously summed RSS across the
  `parallel::mclapply()` fork tree, which double-counts copy-on-write
  shared pages a SLURM cgroup only charges once -- measured to overstate
  the cgroup-enforced peak by 1.8-2.1x at `n_cores` 2/4. The poller now
  sums Pss (`/proc/<pid>/smaps_rollup`) and additionally reads the
  cgroup's own peak-usage counter on Linux; on other platforms it falls
  back to the old sum-RSS figure, now explicitly labeled an upper bound.
  `.BILATR_CHUNK_CORES_STEP_FACTOR`/`.BILATR_CHUNK_CORES_PER_CORE_FACTOR`
  are not yet re-fit from this corrected methodology (needs a Linux run;
  see their updated docs) -- treat them as a documented-conservative
  placeholder for now. `.BILATR_CHUNK_BASELINE_MB` (previously a
  hardcoded 200) is now measured at runtime from the calling process's
  own memory footprint, floored at 200 for platforms where that fails.
* Added `.bilatr_worker_tradeoff()`: reports estimated wall time and
  core-seconds across candidate `n_workers` for a CSV-path chunked
  sweep, since chunk count rises with `n_workers` (chunk size shrinks to
  fit the same memory budget) and each chunk is a full re-parse of the
  largest chain file -- past some point more workers make both wall
  time and core-seconds worse, not just core-seconds. `diagnose_convergence()`/
  `extract_theta()`/`diagnose_and_extract_bilatr()` gain an optional
  `read_seconds` argument (from the caller's own logs) that, when
  supplied, checks a small neighborhood around the chosen `n_workers`
  and names a level that would give both lower wall time and lower
  core-seconds, if one exists.

# bilatr 0.4.1

## Bug fixes

* Removed the last intermediate file from the CSV-file-path read path.
  0.3.10 stopped copying through `tempdir()` but still wrote one
  comment-stripped `*.csv.nocomments` copy per chain file (by default
  alongside the source file, or under `scratch_dir`); every
  CSV-file-path read (`diagnose_convergence()`,
  `diagnose_and_extract_bilatr()`, `extract_theta()`, `extract_alpha()`,
  `extract_mu_intercept()`) now reads `data.table::fread(file = ...,
  skip = ..., nrows = ..., select = ...)` directly against the original
  CmdStan CSV, with no copy of any kind made at any point. `scratch_dir`
  is accordingly deprecated (accepted, with a warning, and ignored) on
  all five functions. Any leftover `*.csv.nocomments` files from a
  0.3.10 or 0.4.0 run are no longer read or written and can be deleted.
* Fixed a sign-orientation bug (B5) in the chunked Tier 3 read: a chunk
  containing both `theta`/`theta_raw` (which the `stable`/`ou`
  reflection symmetry flips) and `log_lik[d,t]` (which it must not,
  since `log_lik` isn't a location parameter) previously negated every
  column in the chunk uniformly whenever any of its columns needed
  flipping. Only the columns [.bilatr_flip_variables()] actually lists
  are negated now, so a `compute_log_lik = 1` fit's `log_lik` is never
  corrupted by orientation regardless of how it happens to fall into
  chunks alongside `theta`/`theta_raw`.
* Fixed a bug (B1) where an unrecognized `stan_model` name silently
  disabled sign orientation instead of raising an error, in
  [bilatr_orient()], `.warn_if_wrong_basin()` (`fit_dyad_ts()`/
  `fit_panel()`'s post-sampling basin check), and every function above
  that accepts `stan_model`. Unknown names now error immediately.
  Pre-0.4.0 model names (`"alphanorm"`, `"alphanorm_ou"`) are still
  accepted everywhere `stan_model` is read, resolved to their current
  names (`"stable"`, `"ou"`), with a message that now fires once per
  session per alias name rather than on every call. This includes the
  fitting path: `fit_dyad_ts_dev()`/`fit_panel_dev(stan_model =
  "alphanorm")` previously died inside the model-specific init
  generator (which has no entry for the old name) despite resolving the
  Stan file and the post-sampling basin check correctly; `stan_model` is
  now canonicalised once, ensuring all three agree.
* A chain file with fewer post-warmup rows than expected (e.g. a SLURM
  job that died mid-sampling) is now detected up front, by name, with a
  clear error (B7) -- previously this surfaced as `posterior`'s generic
  "must have the same length" error only once chunks were later
  assembled into a draws array.
* `diagnose_convergence()`'s in-memory branch, given a `CmdStanMCMC`-
  like fit object (as opposed to a `posterior::draws` object already in
  hand), now reads variable names via `$metadata()$variables` and calls
  `$draws(variables = keep_vars)` for only the requested `tiers` (B8),
  instead of reading every variable via `$draws()` and subsetting
  afterward -- `tiers = 1` no longer touches Tier 3 in memory at all,
  matching what the CSV-file-path branch already did.

## Changes

* `parallel`/`n_workers` (on `diagnose_convergence()`,
  `diagnose_and_extract_bilatr()`, and `extract_theta()`'s CSV-path
  mode) now drive `posterior::summarise_draws()`'s `.cores` argument
  over each chunk's Rhat/rank-normalised-ESS computation, not the read:
  reading a chunk is disk-bound and single-threaded regardless (see
  `.fast_read_post_warmup_draws()`), so chunks are always read strictly
  sequentially. `n_workers` now defaults to
  `parallelly::availableCores()` rather than
  `parallel::detectCores()`-based logic, so it respects a SLURM
  allocation's `SLURM_CPUS_PER_TASK` instead of reporting the whole
  node. **This is a genuine memory/wall-time trade-off, not a free
  choice**: forking `n_workers` processes for the summary step measures
  as costing memory roughly proportional to `n_workers` (not a fixed
  amount regardless of it), and that cost comes out of the same
  `max_memory_mb` budget the read uses -- so a larger `n_workers`
  indirectly means MORE, not fewer, passes over each chain file at a
  fixed `max_memory_mb`. Check the resolved chunk count in the
  pre-flight `message()` before committing a long run to a large
  `n_workers`.
* The `max_memory_mb`/`chunk_size` memory model was re-derived against
  measured peak RSS (`dev/bench_memory.R`, not shipped with the
  package), not just re-derived from reading the new code, and now
  accounts for two things the initial 0.4.1 model missed: a fixed
  ~200 MB floor for R and its loaded packages (`baseline`, alongside
  the existing file-sized term -- neither shrinks with `chunk_size`),
  and the per-worker forking cost above. The pre-flight `message()`
  reports the resolved chunk count/size and the estimated peak broken
  into its baseline, file-sized, and per-chunk terms separately.
* Added a minimum version bound, `posterior (>= 1.0.0)`, for the
  `.cores` argument this release relies on. Removed `future` from
  `Imports` (unused since the `parallel`/`n_workers` change above,
  which stopped this package's own code from calling it); `furrr`
  remains, still used by `ingest_icews()`.

# bilatr 0.4.0

## Breaking changes

* Promoted the experimental `alphanorm` model to `stable`: it is now the
  model `fit_dyad_ts()`/`fit_panel()`/`compile_bilatr_model()` fit by
  default, replacing the previous `stable` (the consolidated
  Dirichlet-multinomial model with `alpha[1]` hard-fixed to 1 and a
  `mu_theta0`-anchored `theta0` location), which is retired to
  `inst/stan/legacy/bilatr_dirmult_irt_pre_0.4.0.stan` (gitignored, kept
  for local reference only, no longer registered/fittable). Renamed the
  experimental `alphanorm_ou` model to `ou`, still experimental,
  replacing the previous `ou` (an OU/AR(1) variant of the old stable
  model), similarly retired to
  `inst/stan/legacy/bilatr_ou_pre_0.4.0.stan`. The promoted `stable`/`ou`
  identify `alpha`/`mu_intercept` via `sum_to_zero_vector`s with a soft
  sign anchor on `alpha[1]` (new `anchor_scale`/`compute_log_lik`
  arguments to [assemble_stan_data()], already present since 0.3.x for
  the experimental variants) rather than hard-fixing `alpha[1] = 1`; this
  leaves an exact alpha/theta reflection symmetry that
  [bilatr_orient()] -- already applied by default in
  [extract_theta()]/[extract_alpha()]/[extract_mu_intercept()] -- now
  corrects for the default model too, not just the experimental one.
  `theta0`'s population mean is pinned at exactly 0 (no separate
  `mu_theta0`). CSVs from a fit made under the pre-0.4.0 `stable`/`ou`
  are unaffected and remain fully readable by
  [diagnose_convergence()]/[extract_theta()]/etc. (nothing about reading
  raw CmdStan CSVs depends on the registry), but code assuming
  `mu_theta0`, `mu_intercept_raw`, `alpha[1] == 1`, or
  `mu_intercept[1] == 0` will need updating for any model *fit* under
  this release's `stable`/`ou`, since `stan_model = "stable"`/`"ou"` now
  names a different Stan program than before.

# bilatr 0.3.10

## Bug fixes

* Fixed a memory-model bug in every CSV-file-path read
  (`diagnose_convergence()`, `diagnose_and_extract_bilatr()`,
  `extract_theta()`, `extract_alpha()`, `extract_mu_intercept()`) that
  could OOM-kill a job however small `max_memory_mb`/`chunk_size` was
  set, or however many cores `parallel = TRUE` was given: internally,
  every one of these previously called
  `cmdstanr::read_cmdstan_csv(csv_files, variables = ...)`, which reads
  each chain's CSV via `data.table::fread(cmd = "grep -v '^#' <file>")`
  -- and per `?data.table::fread`, a `cmd=`/piped input is always
  written to a full temporary copy in `tempdir()` before being read "as
  normal," regardless of how few `variables` are requested. That
  temp-file write happened on *every* call -- once per Tier 3 chunk, not
  once per run -- so a sweep with hundreds of chunks repeated a
  near-complete copy of each (often many-GB) chain file hundreds of
  times; and if `tempdir()`/`$TMPDIR` resolves to a RAM-backed `tmpfs`
  (a common per-node HPC/SLURM default), each copy was a direct hit
  against the job's memory allocation, independent of anything
  `max_memory_mb`/`chunk_size`/`n_workers` controlled.
* All of these now read via two new internal helpers,
  `.prepare_fast_csv_read()`/`.fast_read_post_warmup_draws()`: each
  chain file has its comment lines stripped exactly once per call (via
  a portable, streaming `readr::read_lines_chunked()` pass -- no
  external `grep` dependency), written next to the source file by
  default rather than through `tempdir()`, and every subsequent
  Tier 1/2/3 read/chunk reads `variables` straight from that one cleaned
  file via `data.table::fread(file = ...)` (a real path, so `fread` can
  select columns without buffering the whole file), verified to match
  `cmdstanr::read_cmdstan_csv()`'s output exactly. This restores the
  memory bound `max_memory_mb`/`chunk_size` was always meant to provide,
  makes `parallel = TRUE` safe to use at its documented memory cost
  again, and removes the redundant per-chunk file-copying regardless. A
  new `scratch_dir` argument (default `NULL`, meaning alongside each
  source file) lets you redirect the one-time cleaned copies elsewhere
  if needed.

# bilatr 0.3.9

## New features

* Added `diagnose_and_extract_bilatr()`, the single-shared-read follow-up
  0.3.8 anticipated: given raw CmdStan CSV file paths, it reads Tier 1/2
  once and Tier 3 once (chunked, same as `diagnose_convergence()`),
  reorients whichever columns `bilatr_orient()` would flip as part of
  that one read, and derives `diagnose_convergence()`,
  `extract_theta()`, `extract_alpha()`, and `extract_mu_intercept()`'s
  outputs from those two summaries -- rather than calling the four
  public functions independently, which re-reads `theta`'s columns
  twice (once for Tier 3 diagnostics, once for `extract_theta()`) and
  pays a full separate file-scan for each of `extract_alpha()`/
  `extract_mu_intercept()`, on top of `diagnose_convergence()`'s own
  Tier 1/2 read of the same columns. Meaningful for CSV-path,
  I/O-bound reads specifically (see 0.3.8); an in-memory `fit` has no
  repeated-read cost to fuse away, so it isn't supported here.

# bilatr 0.3.8

## Changes

* A benchmark against a production-scale CmdStan CSV (~1.9M Tier 3
  columns) found `read_cmdstan_csv()` strongly I/O-bound, not
  parsing-bound: per-chunk wall-time barely depends on how many
  variables are requested (a 1000x range in chunk size changed per-call
  time by under 10%), so total sweep time scales with chunk COUNT, not
  memory saved. `diagnose_convergence()`'s chunking is unchanged, but
  its default-`max_memory_mb` message now explains this tradeoff
  explicitly, so a smaller budget isn't mistaken for a "safer" choice.
* `extract_theta()`'s `fit` argument now also accepts a character vector
  of raw CmdStan CSV file paths, reading/summarising `theta` in the same
  memory-bounded, optionally-parallel chunks as
  `diagnose_convergence()`'s CSV path (new `max_memory_mb`/`chunk_size`/
  `parallel`/`n_workers` arguments, sharing its underlying helper so a
  worker script can eventually feed both from one pass over the CSVs
  instead of two). Sign orientation (`bilatr_orient()`) is applied to
  raw draws before summarising, per chunk, matching the in-memory path
  exactly rather than adjusting already-computed quantiles post hoc.
  Only supports the default `probs`. `extract_alpha()`/
  `extract_mu_intercept()` also gained CSV-path support (no chunking --
  Tier 1/2 is small regardless of dyad-set size), for cross-chain
  roll-up scripts run after several independently-submitted SLURM jobs
  have completed, with no in-memory fit available.

# bilatr 0.3.7

## New features

* Added preliminary human-readable labels for `EventRootCode2`
  (`eventrootcode2_name()`, matching the existing `bilatr_class_name()`/
  `bilatr_class2_name()` convention), attached as a new
  `EventRootCode2Name` column in `cameo_lookup` and by `recode_cameo()`.
  The two root-04 splits (`"044"`/`"046"`) reuse `bilatr_class_name()`'s
  text for the same category; every other label is new, since
  `EventRootCode2` is finer-grained than `BilatrClass` elsewhere.

# bilatr 0.3.6

## Changes

* `diagnose_convergence()`'s `fit` argument now also accepts a character
  vector of raw CmdStan CSV file paths (one per chain), for diagnosing a
  completed run without loading it into memory first via
  `cmdstanr::as_cmdstan_fit()`. In this mode, Tier 1/2 variables are read
  in one small call as before, but Tier 3 (`theta`/`theta_raw`, typically
  the overwhelming majority of monitored quantities for production-sized
  panels) is read and summarised in memory-bounded chunks via
  `cmdstanr::read_cmdstan_csv()`'s `variables` argument, discarding each
  chunk's draws before reading the next -- the full draws array is never
  materialized at once. New arguments `max_memory_mb` (default `8192`;
  drives an automatically-derived chunk size, so callers don't have to
  guess a variable count directly), `chunk_size` (explicit override),
  `parallel` (process chunks concurrently via `furrr::future_map_dfr()`;
  trades the sequential path's memory bound for wall-clock speed, and
  falls back from `future::multicore` to `future::multisession` with a
  warning on Windows), and `n_workers`. Existing calls (an in-memory fit
  or draws object, the previous signature) are unaffected -- the new
  arguments are unused on that path.

# bilatr 0.3.5

## Changes

* `alphanorm`/`alphanorm_ou`'s `alpha[1]` sign orientation is now
  deterministic across runs: their soft sign anchor alone could not
  guarantee it (it corrects relative posterior mass between the two
  mirror modes, but the modes are separated by a likelihood barrier a
  single chain essentially never crosses, so whichever mode a chain's
  init landed in was the one it reported). `bilatr_init_fn()` now biases
  initial values toward the canonical (`alpha[1] > 0`) mode and gives the
  latent states real initial spread, `fit_bilatr()` warns after sampling
  if a fit still came back in the wrong mode, and the new
  `bilatr_orient()` (`R/orient.R`) deterministically relabels a fit's
  draws to the canonical orientation regardless of which mode the
  sampler found -- now wired into `extract_theta()`/`extract_alpha()`/
  `extract_mu_intercept()` via their new `stan_model` argument, so no
  consumer of those functions sees an unoriented fit. Fits produced
  before this change may be sign-flipped relative to fits produced after
  it; `bilatr_orient()` relabels them to match.

# bilatr 0.3.4

## Changes

* `alphanorm`/`alphanorm_ou` gained a soft sign anchor on `alpha[1]`
  (`target += log_inv_logit(alpha[1] * inv(anchor_scale))`, new
  `anchor_scale` data field on `assemble_stan_data()`, default `0.1`),
  breaking an exact reflection symmetry those two models' identification
  (`sum_to_zero_vector` alpha with no fixed element) otherwise leaves in
  place: negating `alpha` together with `theta` and its upstream raw
  parameters leaves the likelihood and every prior unchanged, so the
  posterior has two mirror modes of equal mass and chains could land in
  either, making Rhat on `alpha`/`theta` uninterpretable. The anchor
  orients the positive mode as canonical, so higher `theta` means better
  (less hostile) relations, matching `stable`/`ou`. `ou` is unaffected
  (its hard `alpha[1] = 1` already selects a mode) and does not get this
  data field.

# bilatr 0.3.3

## New features

* Added three EXPERIMENTAL Stan model variants (`R/model_registry.R`,
  `status = "experimental"`), fittable only via
  `fit_dyad_ts_dev()`/`fit_panel_dev(stan_model = ...)`; the stable model
  is unchanged and remains the only one the exported fitters use.
  * `alphanorm` (`inst/stan/bilatr_alphanorm.stan`) targets the residual
    affine ridge in stable's identification: with `alpha[1] = 1` and
    `mu_intercept[1] = 0`, location and scale are only softly pinned, and
    an affine reparameterization of `theta`/`alpha`/`mu_intercept` leaves
    the likelihood exactly invariant -- inflating posterior intervals and
    likely contributing to treedepth pathology. `alphanorm` pins location
    hard (`mu_theta0` removed) and pins scale on the alpha side instead
    of theta's (`sum_to_zero_vector` alpha, normalized to RMS 1).
  * `ou` (`inst/stan/bilatr_ou.stan`) targets the stable model's random
    walk having no restoring force: cross-dyad SD at t=1 is small
    relative to a panel's worth of accumulated drift, so any dyad can
    traverse the entire cross-dyad range within a decade and initial
    ordering carries little information about later ordering -- this is
    why dyads that should be clearly separated in level can invert during
    a temporary detente. `ou` replaces the random walk with an OU/AR(1)
    process with dyad-specific equilibria and a global persistence `rho`,
    giving cross-dyad ordering a permanent component.
  * `alphanorm_ou` (`inst/stan/bilatr_alphanorm_ou.stan`) combines both:
    `alphanorm`'s identification with `ou`'s dynamics.
  * All three are new and UNCALIBRATED: several priors carry over from
    stable with materially changed units/meaning (documented per-model in
    each `.stan` file's header comment), and none has been checked
    against prior predictive simulation yet.

## Changes

* `assemble_stan_data()` gained `rho_prior_a`/`rho_prior_b` (default `8`,
  `2`) and `compute_log_lik` (default `0`) arguments, threaded into the
  returned data list. Consumed only by the three new experimental
  variants (`rho_prior_a`/`rho_prior_b` by `ou`/`alphanorm_ou`;
  `compute_log_lik` gates a per-dyad-period `log_lik` in `generated
  quantities` for all three, off by default since it is `D x T x` draws);
  ignored by `stable`, so existing callers are unaffected.
* `bilatr_init_fn()` now branches on `stan_model` to build
  parameterization-appropriate inits for each registered model, instead
  of always returning stable's init list regardless of the argument.

## Internal

* The `reduce_sum` likelihood (`partial_log_lik`) shared by all four
  registered models is no longer duplicated inline in each `.stan` file.
  Its single canonical source is
  `inst/stan/include/partial_log_lik.stanfunctions`;
  `data-raw/sync_stan_functions.R` splices it into a marker-delimited
  block in each model file (`R/stan_includes.R`), checked for drift by
  `tests/testthat/test_stan_includes.R`. This is NOT a Stan `#include`:
  cmdstanr has a bug (upstream issue
  [stan-dev/cmdstanr#820](https://github.com/stan-dev/cmdstanr/issues/820),
  fixed but not yet released) where `#include` resolution breaks at
  `$sample()`-time whenever the include path contains a space, which
  bites this project's own `devtools::load_all()` working tree.

# bilatr 0.3.2

## New features

* Added `BilatrClass2`, a 9-level (0-8) coarsening of `BilatrClass` that
  merges the two adjacent pairs of hostile levels: `BilatrClass` 7
  ("Investigate, demand, reject, or reduce relations") + 8 ("Disapprove")
  -> `BilatrClass2` 7 ("Disapprove, demand, reject, or reduce relations"),
  and `BilatrClass` 9 ("Threaten or coerce") + 10 ("Assault, fight, or
  mass violence") -> `BilatrClass2` 8 ("Threaten, coerce, or use force").
  `recode_cameo()` now attaches `BilatrClass2` / `BilatrClass2Name`
  alongside the other recode columns, `cameo_lookup` gained the two
  columns, and new internal helpers `assign_bilatr_class2()` /
  `bilatr_class2_name()` (`R/cameo_recode.R`) encode the mapping.

## Changes

* `recode_cameo()` now checks whether any of the recode columns it would
  add are already present in `data`; those are left untouched (no more
  silent `.x`/`.y` suffixing from the join) and a warning names them.
* Retired the experimental `phi_logn` model. `stable` is now the only
  entry in the internal model registry; `.resolve_stan_model("phi_logn")`
  (and `fit_*_dev(stan_model = "phi_logn")`) now error like any other
  unknown name. Its final source is kept at
  `inst/stan/legacy/bilatr_phi_logn.stan` (gitignored, not shipped), with
  the pre-0.3.0 centered version alongside it as
  `bilatr_phi_logn_pre_0.3.0.stan`. `bilatr_init_fn()` no longer branches
  on model name (no `log_phi0_raw` / `beta_logn` init), and `beta_logn`
  was dropped from `diagnose_convergence()`'s Tier 1 name list. The
  exported fitters were already `stable`-only and are unaffected.

# bilatr 0.3.1

## Internal

* Added internal helpers for the `reduce_sum` threading benchmark sweep
  (`R/benchmark.R`, all `@keywords internal`, not exported):
  `build_reduce_sum_grid()` (sizes a cores × grainsize × replicate grid
  from a scenario's dyad count), `insert_reduce_sum_profile()` (wraps the
  model's `reduce_sum` likelihood call in a Stan `profile()` block),
  `summarise_benchmark()`, `pareto_front()`, and `plot_benchmark()`.
  Covered by `tests/testthat/test_benchmark.R`. The SLURM sweep that
  uses them — re-tuning `threads_per_chain` × grainsize for the `stable`
  model on the politically-relevant vs full dyad sets — lives in
  `runscripts/benchmark/` (not shipped).
* `ggplot2` added to `Suggests` (used by `plot_benchmark()`).

# bilatr 0.3.0

## New features

* `recode_cameo()` now also attaches `EventRootCode2`, `BilatrClass`, and
  `BilatrClassName` columns (alongside the existing
  `QuadClass`/`PentaClass`/`PentaClass_modified`). `EventRootCode2` is a
  coarser regrouping of the CAMEO root codes (root 04 "Consult" split
  three ways; roots 09/14/15/18/20 folded into related roots), and
  `BilatrClass` is an 11-level (0-10) action-class scheme intended as the
  model's default `grouping_var`. New internal helpers
  `assign_eventrootcode2()`, `assign_bilatr_class()`, and
  `bilatr_class_name()` in `R/cameo_recode.R` encode the mapping (derived
  from `original_code/cameo_df.csv`), and the `cameo_lookup` package data
  object gained the three columns.

## Breaking changes

* **Removed the `reference_hostile` argument** from
  `assemble_stan_data()` and `grouped_events_to_dyad_period()`, and the
  `reference_hostile` parameter from the internal `order_event_classes()`.
  The Stan models no longer anchor the known-hostile action class:
  `alpha[A]` is no longer fixed to `-alpha_hostile` (with
  `alpha_hostile > 0`) but is now freely estimated like every other
  `alpha[2:A]`. `alpha[1] = 1` remains the sole discrimination-scale/sign
  reference (with `mu_intercept[1] = 0` for location), which is still
  sufficient for identification. Existing calls that passed
  `reference_hostile = ...` must drop that argument; the resulting class
  ordering changes only in that the previously-last "hostile" class now
  sorts alphabetically with the rest.
* **Promoted the non-centered-`process_noise` models to the primary
  models and retired the previous ones.** `inst/stan/bilatr_dirmult_irt.stan`
  (registry `"stable"`) and `inst/stan/bilatr_phi_logn.stan` (registry
  `"phi_logn"`) now sample the per-dyad `process_noise` lognormal
  hierarchy in non-centered form (`log_process_noise_raw ~ std_normal()`,
  with `process_noise` built in transformed parameters), which broke the
  `process_noise` funnel that was driving non-convergence. The
  transitional `"stable_ncproc"` / `"phi_logn_ncproc"` registry entries
  (added in 0.2.2) are removed. The pre-0.3.0 centered sources are kept
  under `inst/stan/legacy/` (gitignored) for reference only and are not
  registered or shipped.

## Internal

* `bilatr_init_fn()` now initializes `log_process_noise_raw` (not
  `process_noise`) and an `alpha_raw` of length `A - 1` (not `A - 2`),
  and no longer sets `alpha_hostile`, to match the new parameterization.
* Dropped `"alpha_hostile"` from `diagnose_convergence()`'s Tier 1
  parameter-name list (the parameter no longer exists).

# bilatr 0.2.2

* Added two new experimental models

# bilatr 0.2.0

## New features

* Added `diagnose_convergence()` and its `print.bilatr_diagnostics()`
  method (`R/diagnose_convergence.R`), a post-hoc MCMC diagnostic triage
  helper for the panel model. It runs `posterior::summarise_draws()` and
  splits the result into three tiers -- global/shared parameters, per-dyad
  hierarchical parameters, and per-dyad-period latent states -- so that
  Rhat/ESS problems in shared parameters are always surfaced, while the
  wide, poorly-identified posteriors expected for sparse dyads' own
  `theta` trajectories are aggregated rather than flooding the report.
  Per-dyad tiers are screened against a supplied `n_dt` (per-dyad
  observation count) table to flag dyads whose diagnostics are worse than
  their sparsity alone would predict. A `tiers` argument (default `1:3`)
  restricts computation to a subset of tiers -- e.g. `tiers = 1` checks
  only the global/shared parameters, skips computing Rhat/ESS for the
  (typically far more numerous) Tier 2/3 quantities entirely via
  `posterior::subset_draws()`, and needs no `n_dt` argument. See
  `vignette("diagnostics")`.

## Internal

* Removed the internal `relevant_actors()` helper; its value (`c("GOV",
  "MIL", "SPY")`) is now the default for a new `relevant_actors`
  argument on [extract_all_relevant_gdelt()], matching how
  [ingest_icews()] already exposes its own actor filter via
  `relevant_sectors`. Backward compatible (existing calls with no
  argument are unaffected); callers can now override which GDELT
  actor-type codes count as "relevant" without editing package code.
* Renamed all kebab-case source, test, and vignette files to snake_case
  (e.g. `R/model-registry.R` -> `R/model_registry.R`,
  `vignettes/panel-model.Rmd` -> `vignettes/panel_model.Rmd`) for
  consistent file naming across the package. No functional changes;
  `vignette()` cross-references were updated to match the new vignette
  file names.
* Declared `stats` and `utils` in `Imports` (both already used via `::`
  elsewhere in the package, e.g. `stats::quantile()` in
  `extract_theta()`, but not previously declared).

## Development notes (internal, not part of the public API)

* Added an experimental Stan model variant, `inst/stan/bilatr_phi_logn.stan`,
  in which the Dirichlet-multinomial concentration `phi` is modeled as a
  per-dyad-period quantity (`log(phi[d,t]) = log_phi0[d] + beta_logn *
  centered_log_n[d,t]`) rather than a per-dyad constant. Shares the stable
  model's data interface exactly (no changes to `assemble_stan_data()` were
  needed) and nests the stable model's dispersion structure as
  `beta_logn -> 0`.
* Added an internal Stan model registry (`R/model_registry.R`,
  `.bilatr_stan_models`) mapping short names to `inst/stan/` files, plus a
  shared, cached compilation path (`.compile_stan_model()`) used by both
  the exported fitters and the new internal dev entry points below.
  `fit_dyad_ts()` and `fit_panel()` are unaffected: their signatures and
  behavior are unchanged, and they always resolve the registry's
  `"stable"` entry internally.
* Added internal (non-exported) development entry points
  `fit_dyad_ts_dev()`/`fit_panel_dev()` in `R/fit_dev.R`, identical to
  `fit_dyad_ts()`/`fit_panel()` but with a `stan_model` argument for
  selecting a registered model during active development. Not part of the
  public API.
