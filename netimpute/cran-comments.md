# cran-comments

## Submission

This is a new submission (netimpute 0.1.0).

## Test environments

* local Windows 11, R 4.6.1
* (add win-builder devel/release and R-hub results before submitting)

## R CMD check results

0 errors | 0 warnings | 0 notes

## Notes for the reviewer

* The examples of `netmice()`, `complete_netmice()`, `plot.netmids()` and
  `print.netmids()` are wrapped in `\donttest{}`. They are correct and
  runnable, but each one loads the imputation engine (`mice` and its
  dependencies) on first use, which alone takes more than five seconds in
  a fresh session; the imputation itself takes about one second at the
  sizes used.

* The package writes no files, starts no processes, and changes no global
  options or the working directory. `ncores > 1` starts `future`
  multisession workers, and the corresponding test is skipped on CRAN.
