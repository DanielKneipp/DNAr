% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/dsd.R
\name{dsd_4d_make_module_str}
\alias{dsd_4d_make_module_str}
\title{Get a DSD 4domain module string}
\usage{
dsd_4d_make_module_str(template_str, ...)
}
\arguments{
\item{template_str}{The template string of the module;}

\item{...}{The parameters to be setted on the module string.}
}
\value{
A string representing an instantiation of a module with all
parameters set.
}
\description{
Based on a template string, this function prepare it with the
the parameters passed and returns a string representing the
module with all the parameters set.
}
