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
#' binded to the input of neuron 2. The default values of `enzyme_config` and
#' `ci` define a binding which  is capable of activate the next neuron without
#' dependence of any other binding.
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
#'                         the binding output which will be used as input of
#'                         neuron 2).
#' @param bind_inhibitory  If `TRUE`, the inhibitory signal `B` will be binded
#'                         instead of the activatory one `A`.
#'
#' @return  The binding CRN.
#'
#' @export
get_neuron_binding_hje <- function(
    neuron1,
    neuron2,
    enzyme_config = list(ci = 2/3, k = c(200, 100)),
    ci = 4/3,
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
#' @return  The neuron received as input but with its input changed.
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

#' Get a generic neuron gate with n inputs and 1 output
#'
#' This function returns a generic neuron gate. Generic because the
#' number of inputs of this gate will be specified as an argument
#' (length of `input_neuron_names`).
#'
#' A neuron gate is a gate made with neurons, more specifically, 1 + n neurons,
#' where n is the number of inputs. Each input is a binary neuron that passes
#' its output to the last neuron, the output neuron. This neuron has as input
#' the sum of the output of all input neurons.
#'
#' @param gate_name               String specifying the gate name. Used for
#'                                identification purpose.
#' @param input_neuron_name       A list of string. Each string specifies
#' @param output_neuron_name      A string being the output neuron name
#' @param input_neuron_cis        A vector of number representing the
#'                                initial concentration of each neuron
#' @param binding_enzyme_configs  A list where each element represents
#'                                a enzyme configuration. Each configuration
#'                                is a list with the following parameters:
#'                                \itemize{
#'                                  \item `ci`: Initial concentration of the
#'                                    enzyme
#'                                  \item `k`: Vector with the rate constants
#'                                    of forward and backward reactions,
#'                                    respectively.
#'                                }
#' @param binding_cis             Numeric vector with the initial concentration
#'                                of the bindings
#' @param binding_inhibit         Boolean vector where each elements specifies
#'                                whether the binding with the same index will
#'                                be an inhibitory binding or not. If no vector
#'                                is passed, no binding will be inhibited.
#'
#' @return  A list representing the gate. This list comes with three parameters:
#'          \itemize{
#'            \item `output_neuron_crn`: A neuron CRN returned by
#'              `\link{get_neuron_hje}()`;
#'            \item `input_neuron_crns`: A list of neuron CRNs also returned by
#'              `\link{get_neuron_hje}()`;
#'            \item `binding_crns`: A list of binding CRNs returned by
#'              `\link{get_neuron_binding_hje}()`, which are the binding
#'              of each input neuron with the output neuron.
#'          }
#'
#' @seealso `\link{react}()` for more about the neuron behavior.
get_neuron_generic_gate_hje <- function(
    gate_name,
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

    # Check if gate_name is only a string
    assertthat::assert_that(
        is.character(gate_name),
        msg = 'gate_name must be a string'
    )

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

    # TODO: Check the content of binding_enzyme_configs

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
        name = gate_name,
        output_neuron_crn = output_neuron_crn,
        input_neuron_crns = input_neuron_crns,
        binding_crns = binding_crns
    )

    return(gate)
}

#' Get a CRN given a gate
#'
#' This function returns a list specifying a CRN given a list
#' specifying a gate as a parameter, combining all the CRNs of the gate
#' into one single CRN.
#'
#' @param gate  List representing the gate (returned by
#'              `\link{get_neuron_generic_gate_hje}()`, for example).
#'
#' @return  A CRN representing the `gate`.
#'
#' @export
get_crn_from_neuron_gate_hje <- function(gate) {
    # Combine all the gate CRNs
    gate_crn <- combine_crns(c(
        gate$input_neuron_crns,
        list(gate$output_neuron_crn),
        gate$binding_crns
    ))

    return(gate_crn)
}

#' Checks the gate specification
#'
#' This function is used to verify the list that represents a gate is
#' correct defined.
#'
#' @param gate         The list to be verified
check_neuron_gate_hje <- function(gate) {
    # Check if the the gate is a list and it has length > 0
    assertthat::assert_that(
        is.list(gate) && length(gate) > 0,
        msg = paste('The gate', 'must be a list with at least
                    four named elements (name, output_neuron_crn,
                    input_neuron_crns, binding_crns) in it')
    )

    # Check if the the gate has a name field
    assertthat::assert_that(
        is.character(gate$name),
        msg = paste('The gate', gate$name, 'list must have an element
                    called name which is a string representing the gate name')
    )

    # Check if output_neuron_crn is set
    assertthat::assert_that(
        is.list(gate$output_neuron_crn) &&
            length(gate$output_neuron_crn) > 0,
        msg = paste('The gate', gate$name, 'list must have an element
                    with the name output_neuron_crn which is a list
                    with length > 0')
    )

    # Check if output_neuron_crn is set
    assertthat::assert_that(
        is.list(gate$input_neuron_crns) &&
            length(gate$input_neuron_crns) > 0,
        msg = paste('The gate', gate$name, 'list must have an element
                    with the name input_neuron_crns which is a list
                    with length > 0')
    )

    # Check if output_neuron_crn is set
    assertthat::assert_that(
        is.list(gate$binding_crns) &&
            length(gate$binding_crns) > 0,
        msg = paste('The gate', gate$name, 'list must have an element
                    with the name binding_crns which is a list
                    with length > 0')
    )
}

#' Get the binding between two neuron gates
#'
#' This function returns a modified gate of `gate2` input. this modified
#' version has one of its input neurons binded to the output neuron of `gate2`.
#'
#' The binding is made using the function
#' `\link{add_gate_binding_on_neuron_hje}()`. The binding is a list with
#' the following structure:
#'  - `gate1_name`, `gate2_name`: = gate 1 and 2 names, respectively;
#'  - `input_neuron_idx`: the parameter `input_neuron_idx`,
#'  - `neuron_binding_crn`: binding CRN returned by
#'      `\link{get_neuron_binding_hje}()`
#'
#' @param gate1,gate2       The gates which will be binded
#' @param input_neuron_idx  The index of the input neuron of `gate2` which
#'                          will be binded.
#' @param ...               The parameters `enzyme_config`, `ci` and optionally
#'                          `bind_inhibitory`, which will be passed to
#'                          `\link{get_neuron_binding_hje}()`
#'
#' @return  The `gate2` received as input but with the specified input neuron
#'          binded to the `gate1` output.
#'
#' @seealso  `\link{get_neuron_binding_hje}()` for details about the parameter
#'           that it s expecting.
#'
#' @export
define_neuron_gate_binding_hje <- function(
    gate1,
    gate2,
    input_neuron_idx,
    ...
) {
    # Check the gate specification
    check_neuron_gate_hje(gate1)
    check_neuron_gate_hje(gate2)

    # Check if input_neuron_idx is within length(gate2$input_neuron_crns)
    assertthat::assert_that(
        input_neuron_idx <= length(gate2$input_neuron_crns),
        msg = 'input_neuron_idx is out of limits of gate2$input_neuron_crns'
    )

    # Get the input (signal) and output neurons (which will receive the signal)
    gate1_output_neuron <- gate1$output_neuron_crn
    gate2_input_neuron <- gate2$input_neuron_crns[[input_neuron_idx]]

    # Get the binding between these two neurons, which will be the binding
    # between the gates
    neuron_binding_crn <- get_neuron_binding_hje(
        gate1_output_neuron,
        gate2_input_neuron,
        ...
    )

    # Construct the gate binding structure
    binding <- list(
        gate1_name = gate1$name,
        gate2_name = gate2$name,
        input_neuron_idx = input_neuron_idx,
        neuron_binding_crn = neuron_binding_crn
    )

    # Update the neuron structure to add the bindings.
    # This will be used for the CRN construction
    gate2_input_neuron <- add_gate_binding_on_neuron_hje(
        gate2_input_neuron,
        binding
    )

    # Update the neuron on the gate
    gate2$input_neuron_crns[[input_neuron_idx]] <- gate2_input_neuron

    return(gate2)
}

#' Add a gate binding into a neuron
#'
#' Use this function to add a binding into a list called `gate_bindings`
#' within the neuron. this list will be create if it doesn't exists
#'
#' @param neuron   The neuron which will have a gate binding added.
#' @param binding  A list representing the binding. The structure
#'                 can be anything.
#'
#' @return  The neuron passed as parameter but with the binding in it.
add_gate_binding_on_neuron_hje <- function(neuron, binding) {
    # Creates a empty list called gate_bindings in the neuron
    # if there is no such list
    if(is.null(neuron$bindings)) {
        neuron$gate_bindings <- list()
    }

    # Add the binding
    neuron$gate_bindings[[length(neuron$gate_bindings) + 1]] <- binding

    return(neuron)
}

#' Check the gate binding structure
#'
#' This function checks if a given object has the structure of a gate binding.
#'
#' It is checked if the fields `gate1_name`, `gate2_name`, `input_neuron_idx`
#' `neuron_binding_crn` are all set and with the correct types.
#'
#' @param gate_binding  A list representing a gate binding.
check_neuron_gate_binding_hje <- function(gate_binding) {
    # Check if the gate is a list
    assertthat::assert_that(
        is.list(gate_binding) && length(gate_binding) > 0,
        msg = 'The gate binding must be a list and has at least the following
               parameters: gate1_name, gate2_name, input_neuron_idx,
               neuron_binding_crn')

    # Check if the gate has the gate1_name field and it is a string
    assertthat::assert_that(
        is.character(gate_binding$gate1_name),
        msg = 'The gate binding must have a string parameter called
               gate1_name'
    )

    # Check if the gate has the gate2_name field and it is a string
    assertthat::assert_that(
        is.character(gate_binding$gate2_name),
        msg = 'The gate binding must have a string parameter called
               gate2_name'
    )

    # Check if the gate has the input_neuron_idx field and it is a string
    assertthat::assert_that(
        is.numeric(gate_binding$input_neuron_idx),
        msg = 'The gate binding must have a number parameter called
               input_neuron_idx'
    )

    # Check if the gate has the neuron_binding_crn field and it is a list
    assertthat::assert_that(
        is.list(gate_binding$neuron_binding_crn),
        msg = 'The gate binding must have a list parameter called
               neuron_binding_crn (defined by define_neuron_gate_binding_hje())'
    )
}

#' Get a circuit from a list of gates
#'
#' This function checks if the gates and gate bindings are correctly defined
#' (using `\link{check_neuron_gate_hje}()` and `\link{check_neuron_gate_hje}()`
#' functions). Then, it structures all gates into one single list called
#' circuit. The circuit represents all the gates and their interactions,
#' defining one single object containing all the gate network.
#'
#' @param gates  A list of gates.
#'
#' @return  A circuit, which is a list with the following fields in it:
#'  - `gates`: The list of gates.
#'
#' @export
get_circuit_from_neuron_gates_hje <- function(gates) {
    # Check the structure of all gates and gate bindings
    for(gate in gates) {
        check_neuron_gate_hje(gate)

        # Check the gate binding structures
        for(neuron in gate$input_neuron_crns) {
            if(!is.null(neuron$gate_bindings)) {
                for(gate_binding in neuron$gate_bindings) {
                    check_neuron_gate_binding_hje(gate_binding)
                }
            }
        }
    }

    # Check if all gates have different names
    gate_names <- lapply(gates, '[[', 'name')
    assertthat::assert_that(
        !any(duplicated(gate_names)),
        msg = 'All gates must have unique names'
    )

    # Combine the gates and bindings into a list with specific structure
    circuit <- list(
        gates = gates
    )

    return(circuit)
}

#' Checks the circuit structure
#'
#' This function checks if the circuit passed as a parameter has the correct
#' structure.
#'
#' Currently, it is checking if `circuit` is a list and if it has a list
#' called `gates` in it.
#'
#' @param circuit  Circuit object to be checked.
check_neuron_circuit_hje <- function(circuit) {
    # Check if the circuit is a list
    assertthat::assert_that(
        is.list(circuit),
        msg = 'The circuit must be a list'
    )

    # Check if the circuit has the field gates and
    # and if it is a list
    assertthat::assert_that(
        is.list(circuit$gates),
        msg = 'The circuit list must have a list called gates'
    )
}

#' Get a CRN from a circuit
#'
#' This function returns a CRN given a circuit as a parameter.
#'
#' In this function, the gate bindings are correctly set, updating
#' the neuron which will receive the bindings as input. The bindings
#' within the gates are assumed to be already setted.
#'
#' @param circuit  The circuit (returned by
#'                 `\link{get_circuit_from_neuron_gates_hje}()`)
#'
#' @return  A list specifying the CRN.
#'
#' @export
get_crn_from_neuron_circuit_hje <- function(circuit) {
    # Check the circuit structure
    check_neuron_circuit_hje(circuit)

    # Update the gate inputs according to the gate bindings
    circuit$gates <- lapply(circuit$gates, function(gate) {
        new_neurons <- lapply(gate$input_neuron_crns, function(neuron) {
            # If there is no gate binding, the neuron won't change
            if(is.null(neuron$gate_bindings)){
                return(neuron)
            } else {
                new_neuron <- update_neuron_input_hje(
                    neuron,
                    lapply(neuron$gate_bindings, '[[', 'neuron_binding_crn')
                )
                return(new_neuron)
            }
        })

        # Update the input neurons
        gate$input_neuron_crns <- new_neurons
        return(gate)
    })

    # Get all CRNs (gates and bindings)
    crns <- list()
    for(gate in circuit$gates) {
        # Get the gate CRNs (neurons and bindings between neurons within gates)
        crns[[length(crns) + 1]] <- get_crn_from_neuron_gate_hje(gate)

        # Get the binding CRNs
        for(input_neuron in gate$input_neuron_crns) {
            for(gate_binding in input_neuron$gate_bindings) {
                crns[[length(crns) + 1]] <- gate_binding$neuron_binding_crn
            }
        }
    }

    # Combine all CRNs into one single CRN
    g_crn <- combine_crns(crns)

    return(g_crn)
}

#' Get a AND neuron gate
#'
#' This function returns an AND gate using neurons based on McCullogh-Pitts
#' model (neurons returned by `\link{get_neuron_hje}()`). The two input neurons
#' have the name `i1` and `i2`, and the output neuron has the `o` name. This
#' (with the gate name) determines the species names created for the CRNs
#' that represent them.
#'
#' @param gate_name        The name of the gate
#' @param input_cis        A numeric vector representing the initial
#'                         concentration of the inputs.
#' @param binding_inhibit  Boolean vector where each elements specifies
#'                         whether the input with the same index will
#'                         be an inhibitory binding or not. If no vector
#'                         is passed, no binding will be inhibited. This
#'                         parameter is passed to
#'                         `\link{get_neuron_generic_gate_hje}()`.
#'
#' @return  A list representing the gate. To transform this representation into
#'          a CRN (to simulate with `\link{react}()`, for example) you have
#'          to use the function `\link{get_crn_from_neuron_gate_hje}()`.
#'
#' @export
#'
#' @examples
#' g         <- DNAr::get_neuron_AND_gate_hje('_AND', c(0.5, 1.5))
#' and_crn   <- DNAr::get_crn_from_neuron_gate_hje(g)
#' and_crn$t <- seq(0, 1, length.out = 100)
#' b         <- do.call(DNAr::react, and_crn)
#' DNAr::plot_behavior(b, c('A_ANDo'))
get_neuron_AND_gate_hje <- function(gate_name, input_cis, binding_inhibit = NULL) {
    # Define the input and output names
    input1_name <- jn(gate_name, 'i1')
    input2_name <- jn(gate_name, 'i2')
    output_name <- jn(gate_name, 'o')

    # Get the gate
    gate <- get_neuron_generic_gate_hje(
        gate_name = gate_name,
        input_neuron_names = list(input1_name, input2_name),
        output_neuron_name = output_name,
        input_neuron_cis = input_cis,
        binding_enzyme_configs = list(
            list(ci = 1/3, k = c(2e5, 1e5)),
            list(ci = 1/3, k = c(2e5, 1e5))
        ),
        binding_cis = c(2/3, 2/3),
        binding_inhibit = binding_inhibit
    )

    return(gate)
}

#' Get a OR neuron gate
#'
#' This function returns an OR gate using neurons based on McCullogh-Pitts
#' model (neurons returned by `\link{get_neuron_hje}()`). The two input neurons
#' have the name `i1` and `i2`, and the output neuron has the `o` name. This
#' (with the gate name) determines the species names created for the CRNs
#' that represent them.
#'
#' @param gate_name        The name of the gate
#' @param input_cis        A numeric vector representing the initial
#'                         concentration of the inputs.
#' @param binding_inhibit  Boolean vector where each elements specifies
#'                         whether the input with the same index will
#'                         be an inhibitory binding or not. If no vector
#'                         is passed, no binding will be inhibited. This
#'                         parameter is passed to
#'                         `\link{get_neuron_generic_gate_hje}()`.
#'
#' @return  A list representing the gate. To transform this representation into
#'          a CRN (to simulate with `\link{react}()`, for example) you have
#'          to use the function `\link{get_crn_from_neuron_gate_hje}()`.
#'
#' @export
get_neuron_OR_gate_hje <- function(gate_name, input_cis, binding_inhibit = NULL) {
    # Define the input and output names
    input1_name <- jn(gate_name, 'i1')
    input2_name <- jn(gate_name, 'i2')
    output_name <- jn(gate_name, 'o')

    # Get the gate
    gate <- get_neuron_generic_gate_hje(
        gate_name = gate_name,
        input_neuron_names = list(input1_name, input2_name),
        output_neuron_name = output_name,
        input_neuron_cis = input_cis,
        binding_enzyme_configs = list(
            list(ci = 2/3, k = c(2e5, 1e5)),
            list(ci = 2/3, k = c(2e5, 1e5))
        ),
        binding_cis = c(4/3, 4/3),
        binding_inhibit = binding_inhibit
    )

    return(gate)
}
