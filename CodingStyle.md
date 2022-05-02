# Coding Style

## General Rules
Basic paradigms:

* follow the tidyverse-style: https://style.tidyverse.org/
* use the RStudio Addin to style code: https://github.com/r-lib/styler
* only use ascii characters in all files; neve use Umlaut and other characters
* in any file: All names and comments must be based on English


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


## Functions

Named functions should never rely on a variable of the parent environment. Only exception should be anonymous functions.
