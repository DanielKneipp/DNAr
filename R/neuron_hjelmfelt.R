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
    # TODO: Check input

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
get_neuron_binding_hje <- function(
    neuron1,
    neuron2,
    enzyme_config,
    ci,
    bind_inhibitory = FALSE
) {
    # TODO: check the input variables

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
#' @param neuron    The neuron which its input will be replaced
#' @param bindings  List of bindings
#'
#' @return  The neuron received as input but whith its input changed.
#'
#' @export
update_neuron_input_hje <- function(neuron, bindings) {
    # TODO: Check input

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

get_neuron_generic_gate_hje <- function(
    input_neuron_names,
    output_neuron_name,
    input_neuron_cis,
    binding_enzyme_configs,
    binding_cis,
    binding_inhibit = NULL
) {
    #
    # Input checking
    #

    # Check if output_neuron_name is only a string
    assertthat::assert_that(
        is.character(output_neuron_name),
        msg = 'output_neuron_name must be a string'
    )

    # Check if the length of input_neuron_cis and input_neuron_names
    # are equal
    assertthat::assert_that(
        length(input_neuron_names) == length(input_neuron_cis),
        msg = 'The number of input initial concentrations (input_neuron_cis)
               must be equal to the number of input neurons
               (input_neuron_names)'
    )

    # Check if the length of binding_enzyme_configs and input_neuron_names
    # are equal
    assertthat::assert_that(
        length(input_neuron_names) == length(binding_enzyme_configs),
        msg = 'The number of enzyme configurations
               (binding_enzyme_configs) must be equal to the number of
               input neurons (input_neuron_names)'
    )

    # Check if the length of binding_cis and input_neuron_names
    # are equal
    assertthat::assert_that(
        length(input_neuron_names) == length(binding_cis),
        msg = 'The number of initial concentration of the bindings
        (binding_cis) must be equal to the number of
        input neurons (input_neuron_names)'
    )

    # If a binding_inhibit wasn't passed,
    if(is.null(binding_inhibit)) {
        binding_inhibit <- rep(FALSE, length(input_neuron_names))
    } else {
        # Check if the length of binding_inhibit and input_neuron_names
        # are equal
        assertthat::assert_that(
            length(input_neuron_names) == length(binding_inhibit),
            msg = 'The number of binding inhibition specification
            (binding_inhibit) must be equal to the number of
            input neurons (input_neuron_names)'
        )
    }

    # Check if there is any neuron name duplicate
    assertthat::assert_that(
        !any(duplicated(input_neuron_names)),
        msg = 'Input neuron names must be all unique'
    )

    #
    # CRN construction
    #

    # Get the CRNs of the input neurons
    input_neuron_crns <- mapply(function(neuron_name, neuron_ci) {
        get_neuron_hje(neuron_name, neuron_ci)
    }, input_neuron_names, input_neuron_cis, SIMPLIFY = FALSE)

    # Get the output neuron CRN
    output_neuron_crn <- get_neuron_hje(output_neuron_name, 0)

    # Get the bindings
    binding_crns <- mapply(function(neuron_crn, e_conf, ci, inhibit) {
        get_neuron_binding_hje(
            neuron_crn,
            output_neuron_crn,
            e_conf,
            ci,
            inhibit
        )
    }, input_neuron_crns, binding_enzyme_configs,
       binding_cis, binding_inhibit, SIMPLIFY = FALSE)

    # Update input of the output neuron according to the input neurons
    output_neuron_crn <- update_neuron_input_hje(
        output_neuron_crn,
        binding_crns
    )

    # Get the gate specification
    gate <- list(
        output_neuron_crn = output_neuron_crn,
        input_neuron_crns = input_neuron_crns,
        binding_crns = binding_crns
    )

    return(gate)
}

#' @export
get_crn_from_neuron_gate_hje <- function(gate) {
    # Combine all the gate CRNs
    gate_crn <- combine_crns(c(
        list(gate$output_neuron_crn),
        gate$input_neuron_crns,
        gate$binding_crns
    ))

    return(gate_crn)
}

#' Check gates specification
#'
check_gate_hje <- function(gate, gate_id_str) {
    # Check if the the gate is a list and it has length > 0
    assertthat::assert_that(
        is.list(gate) && length(gate) > 0,
        msg = paste('The gate', gate_id_str, 'must be a list with at least
                    three named elements (output_neuron_crn,
                    input_neuron_crns, binding_crns) in it')
    )

    # Check if output_neuron_crn is set
    assertthat::assert_that(
        is.list(gate$output_neuron_crn) &&
            length(gate$output_neuron_crn) > 0,
        msg = paste('The gate', gate_id_str, 'list must have an elemnt
                    with the name output_neuron_crn which is a list
                    with length > 0')
    )

    # Check if output_neuron_crn is set
    assertthat::assert_that(
        is.list(gate$input_neuron_crns) &&
            length(gate$input_neuron_crns) > 0,
        msg = paste('The gate', gate_id_str, 'list must have an elemnt
                    with the name input_neuron_crns which is a list
                    with length > 0')
    )

    # Check if output_neuron_crn is set
    assertthat::assert_that(
        is.list(gate$binding_crns) &&
            length(gate$binding_crns) > 0,
        msg = paste('The gate', gate_id_str, 'list must have an elemnt
                    with the name binding_crns which is a list
                    with length > 0')
    )
}

#' @export
get_gate_binding_hje <- function(gate1, gate2, input_neuron_idx, ...) {
    # Check the gate specification
    check_gate_hje(gate1, 'gate1')
    check_gate_hje(gate2, 'gate2')

    # Check if input_neuron_idx is whitin length(gate2$input_neuron_crns)
    assertthat::assert_that(
        input_neuron_idx <= length(gate2$input_neuron_crns),
        msg = 'input_neuron_idx is out of limits of gate2$input_neuron_crns'
    )

    # Get the input (signal) and ouput neurons (which will receive the signal)
    gate1_output_neuron <- gate1$output_neuron_crn
    gate2_input_neuron <- gate2$input_neuron_crns[[input_neuron_idx]]

    # Get the binding between these two neurons, which will be the binding
    # between the gates
    binding <- get_neuron_binding_hje(
        gate1_output_neuron,
        gate2_input_neuron,
        ...
    )

    return(binding)
}

#' @export
get_neuron_AND_gate <- function(gate_name, input_cis) {
    # Define the input and output names
    input1_name <- jn(gate_name, 'i1')
    input2_name <- jn(gate_name, 'i2')
    output_name <- jn(gate_name, 'o')

    # Get the gate
    gate <- get_neuron_generic_gate_hje(
        input_neuron_names = list(input1_name, input2_name),
        output_neuron_name = output_name,
        input_neuron_cis = input_cis,
        binding_enzyme_configs = list(
            list(ci = 1/3, k = c(200, 100)),
            list(ci = 1/3, k = c(200, 100))
        ),
        binding_cis = c(2/3, 2/3)
    )

    return(gate)
}
