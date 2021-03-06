% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/dsd.R
\name{get_dsd_0eApB_str}
\alias{get_dsd_0eApB_str}
\title{Instantiate a '0 -> A + B' module in the DSD script}
\usage{
get_dsd_0eApB_str(qi, qmax, CiA, CiB, Cmax, A_domains, B_domains)
}
\arguments{
\item{qi}{String representing a variable name for the \code{qi}
parameter of the module signature.}

\item{qmax}{String representing a variable name for the \code{qmax}
parameter of the module signature.}

\item{CiA}{String representing a variable name for the \code{CiA}
parameter of the module signature.}

\item{CiB}{String representing a variable name for the \code{CiB}
parameter of the module signature.}

\item{Cmax}{String representing a variable name for the \code{Cmax}
parameter of the module signature.}

\item{A_domains}{Vector of strings representing the species domains of A.}

\item{B_domains}{Vector of strings representing the species domains of B.}
}
\value{
A string representing the instantiation of a \code{r0_e_ApB()} module.
}
\description{
This function returns a string representing an addition of a
'0 -> A + B' reaction in the DSD script. It creates a \code{r0_e_ApB()} module
in the script, replacing all the parameter strings by the ones
specified in this function.
}
