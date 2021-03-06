% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/neuron_hjelmfelt.R
\name{update_neuron_input_hje}
\alias{update_neuron_input_hje}
\title{Update the neuron input according to the bindings}
\usage{
update_neuron_input_hje(neuron, bindings)
}
\arguments{
\item{neuron}{The neuron which its input will be replaced}

\item{bindings}{List of bindings}
}
\value{
The neuron received as input but with its input changed.
}
\description{
Use this function to update the neuron (returned by
\code{\link{get_neuron_hje}()}) input according to the bindings
that are connect to it. Once this function was used, it can't be used
on the same neuron again since the reaction that will replaced will not
exists.
}
\details{
This function will replace the reaction \code{C -> X1 + C} by multiple
reactions (one for each binding), replacing the \code{C} by the output species
of the binding.
}
