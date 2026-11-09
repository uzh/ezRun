# Coding Style

## General Rules

Basic paradigms:

* follow the tidyverse-style: https://style.tidyverse.org/
* use the RStudio Addin to style code: https://github.com/r-lib/styler
* only use ASCII characters in all files; never use Umlaut and other characters
* in any file: All names and comments must be based on English
* useful suggestions from Sarah Johnson: https://gist.github.com/sj-io/3828d64d0969f2a0f05297e59e6c15ad

Formatting:
* The tab spacing must be 2

## Naming

* use camelCase for variable names
* use underscore_separated for function names
* recommended variable naming scheme:
** counting variables: `nSamples`, `nReads`
** use plural in vectors and matrices: `reads`

## Conventions

* always use curly braces also for very simple `if` statements
* use `library` instead of `require` since `library` gives an error while `require` gives only a warning if the package is missing
* the use of row names as IDs for matrices and data.frames is encouraged
* never access columns/rows by hard-coded indices, be aware that column or row ordering may change, use column or row names if possible
* never user `for (i in 1:length(x))`, if the length of `x` is zero, this fails, always use `for (i in seq_along(x))`
* For logging purposes, use `ezLog` whenever possible. Other functions which can be used to write to stdout and stderr, like `cat`, `print`, `ezWrite`, or `message`, should only be used when it is necessary to parse the logs after the fact.

## Functions

* Named functions should never rely on a variable of the parent environment. Only exception should be anonymous functions.

## Branches, Pull Requests and Versions

* prefix branch names with the type of change: `feat/` for new functionality, `fix/` for bug fixes, e.g. `feat/nfcore-generic-app`, `fix/ena-download-retry`; other prefixes such as `docs/` or `chore/` are fine
* start PR titles the same way: `feat: ...`, `fix: ...`, `feat(ScSeurat): ...`
* do not edit `Version` or `Date` in `DESCRIPTION`; a GitHub Action bumps them after every push to master: a merge from a `feat/` branch (or a `feat:` commit) bumps the minor version (3.23.4 -> 3.24.0), anything else bumps the patch version (3.23.4 -> 3.23.5)
* for a major release, set `Version` in the PR yourself; the Action then leaves it alone
