# sgev

R functions and simulation studies for second-generation e-values in sequential
Gaussian regression. The functions remove cohort-specific nuisance effects,
accumulate the remaining information, and assess practical significance against
a prespecified margin. This repository contains no manuscript, real data,
or precomputed simulation results.

Requires R 4.1 or later. All computations use base R and the recommended
`stats` package. No external R packages or data downloads are needed.

## Install and use

From the repository root:

```sh
R CMD INSTALL .
```

```r
library(sgev)

# One accumulated scalar score, with known error standard deviation.
sgev_scalar(score = 600, information = 5000, margin = 0.05)

# Two predictors: predictor 1 is the target; predictor 2 is a nuisance.
state <- sgev_start(p = 2, target = 1, sigma = 1)
gram <- 300 * matrix(c(1, 0.8, 0.8, 1), 2)
score <- drop(gram %*% c(0.12, 0.4))
state <- sgev_update(state, information = gram, score = score)
sgev_contrast(sgev_fit(state), contrast = 1, margin = 0.05)
```

Pass each new cohort to `sgev_update` once. Its information input is that
cohort's Gram matrix, not a covariance matrix or the cumulative Gram matrix.
In contrast, `sgev_scalar` takes an already accumulated reduced score and
information. The example above uses a deterministic mean score to illustrate
the interface; it is not a simulated data set.

`sgev_box` tests a practical-null box. `sgev_convex` evaluates a feasible dual
bound for a supplied support function. `sgev_information` bounds the Gaussian
information distance to a box, and `sgev_partition` compares local and pooled
information. See `help("sgev")` and `help("sgev_scalar")` for inputs, outputs,
assumptions, and treatment of singular information.

The method assumes Gaussian score increments and known error variance.
Designs and nuisance effects may depend on past observations, subject to the
stated conditional score model. The target effect, scientific margin, and
mixture scale are fixed during monitoring. Strict interval exclusion means
non-negligible; strict containment means negligible; other cases are unresolved.
Rows representing different targets require their own multiplicity allocation.

## Simulations

Run a short execution check or all manuscript settings from the repository root:

```sh
Rscript simulations/run.R smoke
Rscript simulations/run.R full
```

One or more studies can be selected, for example:

```sh
Rscript simulations/run.R full paired drift
```

The runner also accepts its absolute path from another working directory.
Results and figures are created under `output/smoke/` or `output/full/`.
Smoke mode reduces the replication and calibration counts; it is an execution
check, not a replacement for the full studies. Full mode uses these settings:

| Study | Purpose and comparators | Replications |
| :--- | :--- | ---: |
| `monitoring` | Local mixture, scheduled joint SGPV, spending joint SGPV, invariant t-mixture, and calibrated joint region; includes cohort-size comparison | 3,000 |
| `drift` | Local and pooled common-model mixtures and fixed-horizon oracle under common effects and nuisance drift | 5,000 |
| `paired` | First meaningful, negligible, or unresolved outcome; profile and interval decisions checked on each path | 5,000 |
| `stopping` | Local mixture and a supporting likelihood-ratio oracle for three null geometries | 3,000 |
| `near_margin` | Fixed scale, scale mixture, and stitched boundary near the practical margin | 5,000 |
| `adaptive_check` | Full-score and reduced-score likelihood identity at bounded stopping times | 30,000 |

The monitoring calibration uses 200,000 independent paths. The scheduled and
spending SGPV rules use custom simultaneous confidence regions; they are not
implementations of SeqSGPV or PRISM. The invariant t-mixture permits unknown
variance; the other monitoring comparators use known variance. The stopping
oracle uses the true alternative and serves as a benchmark. Paired profiles
and their common confidence sequence give the same decisions, not a power
improvement over one another. Detailed settings and comparator formulas are
in the scripts.

The studies share innovations across methods and record Monte Carlo uncertainty
and horizon censoring. Random-number settings are specified in the scripts.
The runner writes session information with the results. Repeating a study in
the same mode replaces its generated results. Generated files are excluded
from version control and from the submission archive.

## Package checks

```sh
R CMD build .
R CMD check --no-manual sgev_0.2.0.tar.gz
```

Package tests cover nuisance removal, singular information, convex bounds,
scalar formulas, paired decisions, and invalid inputs. Simulation scripts are
kept in the GitHub repository and submission ZIP, not in the installable R
source archive. Run simulations from the repository, rather than the installed
package directory.

## License

Copyright 2026 Siyi Chen. The existing author-controlled license is retained.
Permission to use, modify, or redistribute the package must be obtained from
the copyright holder. See `LICENSE`.
