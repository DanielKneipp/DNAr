% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/io.R
\name{save_reactions_txt}
\alias{save_reactions_txt}
\title{Save the reactions in a formatted text file}
\usage{
save_reactions_txt(species, cis, reactions, kis, filename)
}
\arguments{
\item{species}{The species of all CRN.}

\item{cis}{The initial concentration of the species.}

\item{reactions}{The reactions that will be saved.}

\item{kis}{The rate constant of those reactions.}

\item{filename}{The name of the file (without extension) that will be
saved. The extension \code{.txt} will be added automatically.}
}
\description{
This function can be used for saving the reactions
with their rate constants in a formatted text file.
}
