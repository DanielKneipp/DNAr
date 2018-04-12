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


#
# 4-domain DNA approach
#

#' Check if the reaction is compatible with the ones
#' supported by \code{\link{react_4domain}()}.
#'
#' The reactions supported by \code{\link{react_4domain}()} are
#' only the uni or bimolecular. This function checks if the
#' \code{reaction} parameter meets this requirement, returning
#' \code{TRUE} if the reaction is supported, or \code{FALSE}
#' otherwise.
#'
#' @examples
#' DNAr:::check_reaction_4domain('A + B -> C')   # Should return TRUE
#' DNAr:::check_reaction_4domain('2A -> B')      # Should return TRUE
#' DNAr:::check_reaction_4domain('2A + B -> C')  # Should return FALSE
check_reaction_4domain <- function(reaction) {
    left_part <- get_first_part(reaction)
    sto <- get_stoichiometry_part(left_part)
    return(sto < 3)
}

#' Get buffer modules
#'
#' This function is used to add buffer modules according
#' to the theory described by Soloveichik D et al. `[1]`.
#' The parameters of this function follows the same
#' semantics of \code{\link{react_4domain}()}.
#'
#' @return \code{NULL} if no buffer modules were added. Otherwise
#' it returns a list with:
#'   - `lambda_1`      = lambda^{-1} value;
#'   - `new_species`   = vector with the new species added;
#'   - `new_cis`       = vector with the initial concentrations;
#'   - `new_reactions` = vector with the new reactions;
#'   - `new_ks`        = constant rate of the new reactions.
#'
#' @references
#'   - `[1]` \insertRef{soloveichik2010dna}{DNAr}
get_buff_modules <- function(reactions, ki, qmax, cmax) {
    sigmas <- list()
    bff_aux_species <- c('LS', 'HS', 'WS')
    #             LS   HS  WS
    bff_cis <- c(cmax, 0, cmax)

    # Calculate the sigma for each species.
    # Formation reactions (0 -> A) are ignored
    uni_count <- 0
    for(i in 1:length(reactions)) {
        reactants <- get_reactants(reactions[i])
        first_reactant <- reactants[[1]]

        if(is_bimolecular(reactions[[i]])) {
            if(is.null(sigmas[[first_reactant]])) {
                sigmas[first_reactant] <- ki[[i]]
            } else {
                sigmas[first_reactant] <- sigmas[[first_reactant]] + ki[[i]]
            }
        } else if(!is_formation(reactions[[i]])) {
            uni_count <- uni_count + 1
            if(is.null(sigmas[[first_reactant]])) {
                sigmas[first_reactant] <- 0
            }
        }
    }

    # There is no sigma (there is only formation or unimolecular reactions)
    if(length(sigmas) == 0) {
        return(NULL)
    }

    reaction_sigma <- max(unlist(sigmas))
    lambda_1 <- qmax / (qmax - reaction_sigma)

    new_ks <- c()
    new_bff_reactions <- c()
    new_species <- c()
    new_cis <- c()
    for(i in 1:length(sigmas)) {
        if(!(sigmas[[i]] == reaction_sigma)) {
            qs <- lambda_1 * (reaction_sigma - sigmas[[i]])
            aux_specs <- paste(bff_aux_species, as.character(i), sep = '')
            input_spec <- names(sigmas)[i]

            forward_reaction <- paste(input_spec, '+', aux_specs[1], '-->',
                                      aux_specs[2], '+', aux_specs[3])
            backward_reaction <- paste(aux_specs[2], '+', aux_specs[3], '-->',
                                       input_spec, '+', aux_specs[1])

            new_bff_reactions <- c(new_bff_reactions,
                                   forward_reaction,
                                   backward_reaction)
            new_ks <- c(new_ks, qs, qmax)
            new_species <- c(new_species, aux_specs[1],
                             aux_specs[2],
                             aux_specs[3])
            new_cis <- c(new_cis, bff_cis)
        }
    }

    # If even with bimolecular reactions, there is no buffer module
    # to add.
    if(length(new_ks) == 0) {
        return(NULL)
    } else {
        ret <- list(
            lambda_1 = lambda_1,
            new_species = new_species,
            new_cis = new_cis,
            new_reactions = new_bff_reactions,
            new_ks = new_ks
        )
        return(ret)
    }
}

#' Translate a formal CRN into a set of DNA based reactions according to the
#' approach described by Soloveichik D et al. `[1]`
#'
#' This function is used to simulate a chemical circuit based on DNA made
#' to behavior as expected from a CRN specified by the parameters. In another
#' words, given the CRN^* passed through the parameters, another \eqn{CRN_2} is
#' created based on reactions between strands of DNA. CRN_2 is simulated using
#' \code{\link{react}()}. The matrix behavior of CRN_2 is returned but only
#' of the species specified in the \code{species} parameter, the behavior of
#' the  auxiliary ones are not returned. The parameters of these functions
#' follows the same pattern of \code{\link{react}()}, only with some additions
#' required by this approach `[1]` (here named as 4-domain).
#'
#' @section Known limitations:
#'   - It only support uni or bimolecular reactions;
#'   - Because of \code{\link{react}()} known limitation, this function also
#'   doesn't support bidirectional reactions;
#'   - The species names `L`, `H`, `W`, `O`, `T`, `G`, `LS`,
#'   `HS`, `WS` with or without numbers after it are not supported because
#'   these are the reserved for the auxiliary ones. Ex.: `L2` and `LS2`
#'   are not supported but `LT` and `LT2` are.
#'
#' @param species      A vector with the species of the reaction. The order of
#'                     this vector is important because it will define the
#'                     column order of the returned behavior. The species names
#'                     `L[0-9]*`, `H[0-9]*`, `W[0-9]*`, `O[0-9]*`,
#'                     `T[0-9]*`, `G[0-9]*`, `LS[0-9]*`, `HS[0-9]*`,
#'                     `WS[0-9]*` are not supported. For more information
#'                     about this, see the Section of **Known limitations**.
#' @param ci           A vector specifying the initial concentrations of the
#'                     \code{species} specified, in order.
#' @param reactions    A vector with the reactions of the CRN^*.
#' @param ki           A vector defining the constant rate of each reaction
#'                     in \code{reactions}, in order.
#' @param qmax         Maximum rate constant for the auxiliary reactions.
#' @param cmax         Maximum initial concentration for the auxiliary species.
#' @param alpha,beta   Rescaling parameters.
#' @param t            A vector specifying the time interval. Each value
#'                     would be a specific time point.
#' @param auto_buffer  With the default value of `TRUE`, this specifies if
#'                     buffer modules should be generated automatically.
#' @param verbose      Be verbose and print information about the integration
#'                     process with `deSolve::diagnostics.deSolve()`. Default
#'                     value is `FALSE`
#' @param ...          Parameters passed to `deSolve::ode()`.
#'
#' @return A list with the attributes `behavior`, `species`, `ci`, `reactions`
#' and `ki`. These attributes are:
#'   - `behavior`: A matrix with each line being a specific point in the time
#'                 and each column but the first being the concentration of a
#'                 species. The first column is the time interval. The
#'                 behavior of the auxiliary species are also returned.
#'  - `species`  : A vector with all the species used in the reactions.
#'  - `ci`       : The initial concentration of each species.
#'  - `reactions`: All the reactions computed, including the ones generated
#'                 according to the 4-domain approach.
#'  - `ki`       : The rate constants of the reactions.
#'
#' @export
#'
#' @references
#'   - `[1]` \insertRef{soloveichik2010dna}{DNAr}
#'
#' @example demo/main_4domain.R
react_4domain <- function(
    species,
    ci,
    reactions,
    ki,
    qmax,
    cmax,
    alpha,
    beta,
    t,
    auto_buffer = TRUE,
    verbose = FALSE,
    ...
) {
    reactions <- check_crn(species, ci, reactions, ki, t)

    for(r in reactions) {
        if(!check_reaction_4domain(r)) {
            stop(paste('Failed to process reaction', r))
        }
    }

    new_reactions <- c()
    new_species <- c(species)
    new_ks <- c()
    new_cis <- c(ci * beta)

    uni_aux_species <- c('G', 'O', 'T')
    bi_aux_species <- c('L', 'H', 'W', 'O', 'T')

    # Get buffer modules
    buffer_stuff <- NULL
    if(auto_buffer) {
        buffer_stuff <- get_buff_modules(reactions, ki, qmax, cmax)
    }

    # Change cis according to the lambda^{-1} factor
    if(!is.null(buffer_stuff)) {
        for(i in 1:length(new_cis)) {
            new_cis[i] <- new_cis[[i]] * buffer_stuff$lambda_1
        }
    }

    for(i in 1:length(reactions)) {
        new_reactions_for_i <- c()
        new_species_for_i <- c()
        new_ks_for_i <- c()
        new_cis_for_i <- c()

        if(is_bimolecular(reactions[i])) {
            aux <- paste(bi_aux_species, as.character(i), sep = '')

            left_part <- get_first_part(reactions[i])
            l_p_specs <- get_species(left_part)

            # For a species, get its count on left.
            # If the count is 2 (2A), transform it into 2 molecules (A + A).
            # If there is only one species and the reaction is bimolecular,
            # this unique molecule needs to be duplicated
            if(length(l_p_specs) == 1) {
                l_p_specs <- c(l_p_specs, l_p_specs)
            }

            # Get the products string.
            # This time we don't need to expand 2A to A + A
            right_part <- get_second_part(reactions[i])

            # Combine the molecules into reactions
            new_reactions_for_i <- c(paste(l_p_specs[1], '+', aux[1], '-->',
                                           aux[2], '+', aux[3]),
                                     paste(aux[2], '+', aux[3], '-->',
                                           l_p_specs[1], '+', aux[1]),
                                     paste(l_p_specs[2], '+', aux[2], '-->',
                                           aux[4]))

            # Set the new species, initial concentrations and ks
            new_species_for_i <- c(aux[1], aux[2], aux[3], aux[4])
            #                    L    H     B    O
            new_cis_for_i <- c(cmax, 0.0, cmax, 0.0)

            # Recalculate qis according to the buffer module theory
            qi_with_buff <- ki[i] / (alpha * beta)
            if(!is.null(buffer_stuff)) {
                qi_with_buff <- qi_with_buff * buffer_stuff$lambda_1
            }
            new_ks_for_i <- c(qi_with_buff, qmax, qmax)

            # If the right part is only a 0, do not add the last reaction,
            # otherwise:
            if(!isempty_part(right_part)) {
                new_reactions_for_i <- c(new_reactions_for_i,
                                         paste(aux[4], '+', aux[5], '-->',
                                               right_part))
                new_species_for_i <- c(new_species_for_i, aux[5])
                #                                   T
                new_cis_for_i <- c(new_cis_for_i, cmax)
                new_ks_for_i <- c(new_ks_for_i, qmax)
            }
        } else {
            aux <- paste(uni_aux_species, as.character(i), sep = '')

            left_part <- get_first_part(reactions[i])
            l_p_specs <- get_species(left_part)
            right_part <- get_second_part(reactions[i])

            # If the only 'species' is 0, doesn't add the 0 to the reaction
            if(is_formation(reactions[[i]])) {
                new_reactions_for_i <- c(paste(aux[1], '-->', aux[2]))
            } else {
                new_reactions_for_i <- c(paste(l_p_specs, '+', aux[1], '-->',
                                               aux[2]))
            }
            new_species_for_i <- c(aux[1], aux[2])
            #                    G    O
            new_cis_for_i <- c(cmax, 0.0)
            qi_with_buff <- ki[i] / alpha / cmax
            if(!is.null(buffer_stuff)) {
                qi_with_buff <- qi_with_buff * buffer_stuff$lambda_1
            }
            new_ks_for_i <- c(qi_with_buff)

            if(!isempty_part(right_part)) {
                new_reactions_for_i <- c(new_reactions_for_i,
                                         paste(aux[2], '+', aux[3], '-->',
                                               right_part))
                new_species_for_i <- c(new_species_for_i, aux[3])
                #                                   T
                new_cis_for_i <- c(new_cis_for_i, cmax)
                new_ks_for_i <- c(new_ks_for_i, qmax)
            }
        }

        new_species <- c(new_species, new_species_for_i)
        new_reactions <- c(new_reactions, new_reactions_for_i)
        new_ks <- c(new_ks, new_ks_for_i)
        new_cis <- c(new_cis, new_cis_for_i)
    }

    # Add buffer stuff
    if(!is.null(buffer_stuff)) {
        new_species <- c(new_species, buffer_stuff$new_species)
        new_reactions <- c(new_reactions, buffer_stuff$new_reactions)
        new_ks <- c(new_ks, buffer_stuff$new_ks)
        new_cis <- c(new_cis, buffer_stuff$new_cis)
    }

    # Run the reaction
    b <- react(
        species   = new_species,
        ci        = new_cis,
        reactions = new_reactions,
        ki        = new_ks,
        t         = t,
        verbose   = verbose,
        ...
    )

    # Arrange the data to be returned in a list
    result <- list(
        behavior  = b,
        species   = new_species,
        ci        = new_cis,
        reactions = new_reactions,
        ki        = new_ks
    )

    # Return the behavior of all species (including the auxiliary ones),
    # initial concentrations, reactions and rate constants.
    return(result)
}

#' Automatically modify a CRN to make it compatible with the DNA simulation
#'
#' This function automatically set the `qmax`, `cmax`, `alpha`, `beta`, and
#' optionally the parameters `method` and `t` (duration time).
#'
#' @param crn     The CRN to be modified.
#' @param t       Optional parameter used to redefine the duration time.
#' @param method  Optional parameter which sets the solver method.
#'
#' @return  The modified CRN with the nwe parameters, ready to be used with
#'          the function \code{\link{react_4domain}()}.
#'
#' @export
update_crn_4domain <- function(crn, t = NULL, method = NULL) {
    crn$qmax <- max(crn$ki) * 1e4
    crn$cmax <- max(crn$ci) * 1e4
    crn$alpha <- 1
    crn$beta <- 1
    if(!is.null(method)){
        crn$method <- method
    }
    if(!is.null(t)) {
        crn$t <- t
    }

    return(crn)
}
