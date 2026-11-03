# Equivalence checks for 0.12.0

These checks establish that bilatr 0.12.0 builds Stan input `identical()`
to bilatr 0.10.2/0.11.0 for the published 16cat models. That is why the
existing production fits were not rerun.

The folder is tracked in git but excluded from the package build, since
`data-raw/` is in `.Rbuildignore`.

| file | what |
|---|---|
| `equivalence_0.12.0a.R` / `.md` | `prepare_stan_data(aggregate_events(...))` against 0.10.2's `assemble_stan_data()`, on the full production event data, both 16cat specs |
| `equivalence_0.12.0c.R` / `.md` | the same comparison, run through the ported replication scripts' `build_stan_data()` |

Both scripts run from the package root. They need a 0.10.2 worktree at
`../bilatr-0.10.2` (`git worktree add ../bilatr-0.10.2 4a4ffd5`) and the
GDELT event file used by the replication scripts, which is not in the
repository. The `.md` reports are the recorded results; the scripts
regenerate them.

The in-tree regression test that replaces these comparisons in the package
itself is `tests/testthat/test_prepare_stan_data_equivalence.R`. It compares
against frozen reference outputs made by `data-raw/make_equivalence_reference.R`.
