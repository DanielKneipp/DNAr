# DNAr is a program used to simulate formal Chemical Reaction Networks
# and the ones based on DNA.
# Copyright (C) 2017  Daniel Kneipp <danielv[at]dcc[dot]ufmg[dot]com[dot]br>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU Affero General Public License as published
# by the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU Affero General Public License for more details.
#
# You should have received a copy of the GNU Affero General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.


#' Get a neuron based on Hjelmfelt A. et al. `[1]` approach
#'
#' This function returns a neuron based on the model described by
#' Hjelmfelt A. et al. `[1]`.
#'
#' @param name      The neuron name. This will be used for the species names.
#' @param input_ci  A numeric value specifying the initial concentration of
#'                  the input
#'
#' @return  A CRN representing the neuron.
#'
#' @export
#'
#' @references
#'   - `[1]` \insertRef{hjelmfelt1991chemical}{DNAr}
get_neuron_hje <- function(name, input_ci) {
    # Define the species names based on the neuron name
    sn <- list(
        C  = jn('C', name),
        X1 = jn('X1', name),
        X2 = jn('X2', name),
        X3 = jn('X3', name),
        X4 = jn('X4', name),
        A  = jn('A', name),
        B  = jn('B', name)
    )

    # Define the CRN
    neuron <- list(
        species   = c(sn$C,  sn$X1,  sn$X3,  sn$A, sn$B),
        ci        = c(input_ci, 0,      0,      1,    0),
        reactions = c(jn(sn$C, ' -> ', sn$X1, ' + ', sn$C),
                      jn(sn$X1, ' + ', sn$C, ' -> ', sn$C),

                      jn(sn$X1, ' + ', sn$B, ' -> ', sn$A),
                      jn(sn$A, ' -> ', sn$X1, ' + ', sn$B),

                      jn(sn$X3, ' + ', sn$A, ' -> ', sn$B),
                      jn(sn$B, ' -> ', sn$X3, ' + ', sn$A),

                      jn(sn$X3, ' -> 0'),
                      jn('0 -> ', sn$X3)),
        ki        = c(1e2, 1, 5e4, 1, 5e4, 1, 1, 1e2),
        name      = name
    )

    return(neuron)
}

#' Get a binding of two neurons
#'
#' This function generates a CRN which represents a binding between two neurons
#' returned by `\link{get_neuron_hje}()`. The output of the neuron 1 will be
#' binded to the input of neuron 2.
#'
#' The binding is made by creating the reaction `E_{ij} + A_{j} -> C_{ij}`.
#'
#' @param neuron1          The CRN which defines the neuron 1
#' @param neuron2          The CRN which defines the neuron 2
#' @param enzyme_config    A list with the following parameters: \itemize{
#'                           \item `ci`: the initial concentration of the
#'                             enzyme;
#'                           \item `k`: a vector with the rate constants of
#'                             forward and backward reactions, respectively.
#'                         }
#' @param ci               The initial concentration of the C species (the
#' the binding output which will be used as input of neuron 2).
#' @param bind_inhibitory  If `TRUE`, the inhibitory signal `B` will be binded
#' instead of the activatory one `A`.
#'
#' @return  The binding CRN.
#'
#' @export
get_binding_hje <- function(
    neuron1,
    neuron2,
    enzyme_config,
    ci,
    bind_inhibitory = FALSE
) {
    # Get the neuron names
    n1n <- neuron1$name
    n2n <- neuron2$name

    # Define the species
    sn <- list(
        input = if(bind_inhibitory) jn('B', n1n) else jn('A', n1n),
        E     = jn('E', n2n, n1n),
        Cn1n2 = jn('C', n2n, n1n)
    )

    # Define the CRN
    binding <- list(
        species   = c(sn$E, sn$Cn1n2),
        ci        = c(enzyme_config$ci, ci),
        reactions = c(jn(sn$E, ' + ', sn$input, ' -> ', sn$Cn1n2),
                      jn(sn$Cn1n2, ' -> ', sn$E, ' + ', sn$input)),
        ki        = c(enzyme_config$k[[1]], enzyme_config$k[[2]])
    )

    return(binding)
}

#' Update the neuron input according to the bindings
#'
#' Use this function to update the neuron (returned by
#' `\link{get_neuron_hje}()`) input according to the bindings
#' that are connect to it. Once this function was used, it can't be used
#' on the same neuron again since the reaction that will replaced will not
#' exists.
#'
#' This function will replace the reaction `C -> X1 + C` by multiple
#' reactions (one for each binding), replacing the `C` by the output species
#' of the binding.
#'
#' @param neuron  The neuron which its input will be replaced
#' @param ...     The bindings
#'
#' @return  The neuron received as input but whith its input changed.
#'
#' @export
update_neuron_input_hje <- function(neuron, ...) {
    bindings <- list(...)

    # Updating species (remove the input species, since it will come
    # from the binding)
    old_input <- neuron$species[[1]]
    neuron$species <- neuron$species[2:length(neuron$species)]

    # Updating Cis (remove the Ci of the input species)
    neuron$ci <- neuron$ci[2:length(neuron$ci)]

    # Updating reactions
    new_reactions <- sapply(bindings, function(bind) {
        c(
            stringr::str_replace_all(
                neuron$reactions[[1]],
                old_input,
                bind$species[[2]]
            ),
            stringr::str_replace_all(
                neuron$reactions[[2]],
                old_input,
                bind$species[[2]]
            )
        )
    })
    neuron$reactions <- c(
        new_reactions,
        neuron$reactions[3:length(neuron$reactions)]
    )

    # Updating kis
    new_kis <- rep(neuron$ki[1:2], length(bindings))
    neuron$ki <- c(new_kis, neuron$ki[3:length(neuron$ki)])

    return(neuron)
}
