# DNAr
[![Build Status](https://travis-ci.org/DanielKneipp/DNAr.svg?branch=master)](https://travis-ci.org/DanielKneipp/DNAr)
[![codecov](https://codecov.io/gh/DanielKneipp/DNAr/branch/master/graph/badge.svg)](https://codecov.io/gh/DanielKneipp/DNAr)

This is a R package developed to simulate formal Chemical Reaction
Networks (CRN) and DNA circuits designed to have the same behavior
of CRN.

Currently there is only one approach to build arbitrary CRNs with DNA
implemented in this package. This approach is described by
Soloveichik et al. **[1]**.

## How to Install It

First, using the R console, install the `devtools` package:

```R
install.packages("devtools")
```

After that, just install the `DNAr` package:

```R
library(devtools)
devtools::install_github('DanielKneipp/DNAr')
```

## How to Use

The two main functions are `react()` and `react_4domain()` to
simulate formal CRNs and the ones based on DNA, respectively. To
more info about them, `help(react)` and `help(react_4domain)`
describes the functions and shows some usage examples

## References

- **[1]** Soloveichik, David, Georg Seelig, and Erik Winfree. "DNA
as a universal substrate for chemical kinetics." Proceedings of the
National Academy of Sciences 107.12 (2010): 5393-5398.
