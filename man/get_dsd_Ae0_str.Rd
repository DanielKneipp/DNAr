% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/dsd.R
\name{get_dsd_Ae0_str}
\alias{get_dsd_Ae0_str}
\title{Instantiate a 'A -> 0' module in the DSD script}
\usage{
get_dsd_Ae0_str(qi, CiA, Cmax, A_domains)
}
\arguments{
\item{qi}{String representing a variable name for the \code{qi}
parameter of the module signature.}

\item{CiA}{String representing a variable name for the \code{CiA}
parameter of the module signature.}

\item{Cmax}{String representing a variable name for the \code{Cmax}
parameter of the module signature.}

\item{A_domains}{Vector of strings representing the species domains of A.}
}
\value{
A string representing the instantiation of a \code{A_e_0()} module.
}
\description{
This function returns a string representing an addition of a
'A -> 0' reaction in the DSD script. It creates a \code{A_e_0()} module
in the script, replacing all the parameter strings by the ones
specified in this function.
}
