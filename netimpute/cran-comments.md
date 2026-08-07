# cran-comments

## Submission

This is a new submission (netimpute 1.0.1).

## Test environments

* local Windows 11, R 4.6.1
* win-builder, R Under development (unstable) (2026-08-04 r90350 ucrt)
* (add win-builder release and macOS builder results before submitting)

## R CMD check results

Local Windows 11, R 4.6.1: 0 errors | 0 warnings | 0 notes.

win-builder (R-devel): 0 errors | 0 warnings | 1 note. The note is the
standard new-submission note:

```
* checking CRAN incoming feasibility ... NOTE
Maintainer: 'Robert W. Krause <robert.w.krause@mailbox.org>'

New submission
```

## Notes for the reviewer

* The examples of `netmice()`, `complete_netmice()`, `plot.netmids()` and
  `print.netmids()` are wrapped in `\donttest{}`. They are correct and
  runnable, but the first of them loads the imputation engine (`mice` and
  its dependencies), which alone takes several seconds in a fresh session;
  the imputation itself takes about one second at the sizes used. The
  example datasets were kept small for this reason.

* Where any example is reported above five seconds, the time is dominated
  by loading a dependency on first use (`mice`, or `sna` for the
  sna-provided centrality measures in `net_measures_full()`), not by the
  computation itself.

* The package writes no files, starts no processes, and changes no global
  options or the working directory. `ncores > 1` starts `future`
  multisession workers, and the corresponding test is skipped on CRAN.
