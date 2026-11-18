# backbone-theory

Code and numerical results for

> C. E. Pedreira. *Backbone theory: inference from unpaired measurement blocks.*

The paper studies an experimental design in which a population is measured in
several blocks, no individual is measured in more than one block, and a common
set of attributes — the *backbone* — is measured in every block. This
repository contains the twelve numerical experiments reported in Section 9 and
Appendix C, in the exact form in which they produced the numbers in the paper.

## Requirements

* MATLAB. The experiments were developed and run on **R2025b**; the constructs
  used require R2016b or later (implicit expansion, `string`, `contains`).
* **Statistics and Machine Learning Toolbox**, for `knnsearch` and `quantile`.

No other toolbox, package or external data is needed. Every experiment
generates its own population from a seeded random stream, so all results are
reproducible to the last digit on any machine.

## Layout

```
exp01/Cod      code for one experiment, self-contained
exp01/Result   its output
...
exp10_12/      the later experiments and the post-processing scripts
run_all.m      driver: runs everything in dependency order
```

One folder per experiment, each self-contained. **The shared helpers are not
identical across experiments**: `bbParams.m` and `bbSample.m` evolved while the
work was done, and each experiment is kept with the version it ran with.
Merging the folders would silently change the first experiments. To run one
experiment, add its `Cod` folder to the path — one at a time:

```matlab
addpath(fullfile(pwd,'exp01','Cod'))
runE1          % full run, as in the paper
runE1(true)    % quick smoke test, about one minute
rmpath(fullfile(pwd,'exp01','Cod'))
```

Each driver writes to the sibling `Result` folder, creating it if absent; see
`bbOutDir.m`. Experiments E1 to E8 accept a `quick` flag that shrinks the
sample sizes to a smoke test. The later scripts accept the number of
replications as their first argument.

`run_all.m` does the whole sequence, handling the path and the order:

```matlab
run_all        % full run; hours, see the runtimes below
run_all(true)  % quick pass over everything
```

Two things to know about the quick pass. It **overwrites** the stored `Result`
files with smoke-test output, since each driver writes to its own `Result`
folder; under version control a checkout brings the stored files back. And it
**skips** `E10b_greedyOnly` and `E11_withinDesign`: a quick run of `runE6`
produces a universe of ten attributes instead of fourteen, and those two are
tied to the full universe — `E11` by the realized backbone it examines, `E10b`
by needing the pairing term of every subset of the complete enumeration. Both
run normally in the full pass, and both stop with an explicit message if given
an `E6.mat` from a quick run.

## What each script produces, and where it appears in the paper

| Script | Paper | Figure / Table | Note |
|---|---|---|---|
| `exp01/Cod/runE1` | 9.3, C.2 | Figure 1, Table 3 | nearest-donor rate, Theorem T1 |
| `exp02/Cod/runE2` | 9.4, C.3 | Figure 7, Table 4 | risk decomposition, Theorem T2 |
| `exp03/Cod/runE3` | 9.5, C.4 | Figure 2, Table 5 | prediction/reconstruction dichotomy |
| `exp04/Cod/runE4` | 9.6, C.5 | Figure 3, Table 6 | non-identifiability made operational |
| `exp05/Cod/runE5` | 9.7, C.6 | Figure 8, Table 7 | limit coupling with several blocks |
| `exp06/Cod/runE6` | 9.8, C.7 | Figure 4, Table 8 | exhaustive enumeration of the criterion |
| `exp07/Cod/runE7` | 9.10, C.8 | Figure 5, Table 9 | completed base, Theorem T5 |
| `exp08/Cod/runE8` | 9.11, C.9 | Figure 6, Table 10 | design under cost and capacity |
| `exp10_12/Cod/baselinesE8` | 9.11, C.9 | Table 11 | cost of omitting one component |
| `exp09/Cod/E9_pilotSelection` | 9.9, C.10 | Table 12 | selection from a finite pilot |
| `exp10_12/Cod/E10_fullyEmpirical` | 9.9, C.10 | Table 13 | every component estimated |
| `exp10_12/Cod/E10b_greedyOnly` | 9.9, C.10 | Table 13 | search cost isolated from geometry cost |
| `exp10_12/Cod/E11_withinDesign` | 9.9, C.10 | Table 14 | selection inside the realized design |
| `exp10_12/Cod/E12_nonGaussian` | 9.11, C.11 | Table 15 | two non-Gaussian populations |
| `exp10_12/Cod/postE12_dispersao` | 9.9, 9.11 | — | dispersion of the excess criterion |
| `exp10_12/Cod/postE4` | 9.6 | — | the bounds of Observation O2 |
| `exp10_12/Cod/postE6` | 9.8 | — | orderings induced by the criterion |
| `exp10_12/Cod/figE1` | 9.3 | Figure 1 | redraws the figure from `E1.mat` |
| `exp10_12/Cod/figE3` | 9.5 | Figure 2 | redraws the figure from `E3.mat` |

The post-processing scripts read the `.mat` file of the experiment they
summarize and locate it themselves; see `bbFindMat.m`. Two dependencies are
worth stating:

* `baselinesE8` needs `exp08/Result/E8.mat`; `E9`, `E10`, `E10b` and `E11` need
  `exp06/Result/E6.mat`; `postE12_dispersao` needs
  `exp10_12/Result/E12.mat`. All four are included, so these scripts run
  without repeating the long simulations.
* `figE1` also needs `bbParams` from `exp01/Cod`, and `E1.mat`, which is **not**
  included (see below). Run `runE1` first; `run_all.m` already handles the path.

## Runtimes

Approximate, on a desktop machine, for the full runs. The quick flag reduces
each of E1 to E8 to about a minute.

| Script | Time |
|---|---|
| `runE1` … `runE5`, `runE7`, `runE8` | minutes each |
| `runE6` | exhaustive enumeration of 2^14 candidates; tens of minutes |
| `E9_pilotSelection` | tens of minutes |
| `E10_fullyEmpirical`, `E11_withinDesign` | tens of minutes each |
| `E10b_greedyOnly`, `postE4`, `postE6`, `postE12_dispersao` | seconds |
| `E12_nonGaussian` | about half an hour |
| `baselinesE8` | seconds; it reuses what `runE8` stored |

`runE8` itself is the long one: the full design study took about sixty hours,
and `E8.mat` is included so that it need not be repeated.

## Results included, and what is not

Each `Result` folder holds the text tables the scripts printed, so the numbers
in the paper can be checked against the output that produced them without
running anything, and the `.mat` files that other scripts consume as input.

Figures are not included, since they appear in the paper, and `E1.mat` is not
included because it is 44 MB; `runE1` regenerates it. The drivers do write both,
so `.gitignore` keeps them out of the repository rather than relying on anyone
remembering to.

## Seeds

Every experiment fixes its seeds explicitly, either through `rng` or through a
named `RandStream`. Two runs of the same script on the same MATLAB version
produce identical output. Where a replication index enters the seed, the
formula is in the script, so a single replication can be reproduced in
isolation.

## License

MIT; see `LICENSE`.
