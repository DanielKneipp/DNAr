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
        } else {
            uni_count <- uni_count + 1
            if(is.null(sigmas[[first_reactant]])) {
                sigmas[first_reactant] <- 0
            }
        }
    }

    # There is no bimolecular reactions
    if(uni_count == length(reactions)) {
        return(NULL)
    }

    reaction_sigma <- max(unlist(sigmas))
    lambda_1 <- qmax / (qmax - reaction_sigma)

    new_ks <- c()
    new_bff_reactions <- c()
    new_species <- c()
    new_cis <- c()
    for(i in 1:length(sigmas)) {
        if(!sigmas[[i]] == reaction_sigma) {
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

    ret <- list(
        lambda_1 = lambda_1,
        new_species = new_species,
        new_cis = new_cis,
        new_reactions = new_bff_reactions,
        new_ks = new_ks
    )
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
#' @param species     A vector with the species of the reaction. The order of
#'                    this vector is important because it will define the
#'                    column order of the returned behavior. The species names
#'                    `L[0-9]*`, `H[0-9]*`, `W[0-9]*`, `O[0-9]*`,
#'                    `T[0-9]*`, `G[0-9]*`, `LS[0-9]*`, `HS[0-9]*`,
#'                    `WS[0-9]*` are not supported. For more information
#'                    about this, see the Section of **Known limitations**.
#' @param ci          A vector specifying the initial concentrations of the
#'                    \code{species} specified, in order.
#' @param reactions   A vector with the reactions of the CRN^*.
#' @param ki          A vector defining the constant rate of each reaction
#'                    in \code{reactions}, in order.
#' @param qmax        Maximum rate constant for the auxiliary reactions.
#' @param cmax        Maximum initial concentration for the auxiliary species.
#' @param alpha,beta  Rescaling parameters.
#' @param t           A vector specifying the time interval. Each value
#'                    would be a specific time point.
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
    t
) {
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
    buffer_stuff <- get_buff_modules(reactions, ki, qmax, cmax)

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

            new_reactions_for_i <- c(paste(l_p_specs, '+', aux[1], '-->',
                                           aux[2]))
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
        t         = t
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

# TODO: Doc this function
get_dsd_header_str <- function(time, species_to_plot) {
    # Set template string
    s <- 'directive duration %t points 1000
directive simulation deterministicstiff
directive concentration M
directive compilation infinite'

    # Set time
    stringr::str_replace(s, '%t', as.character(max(time)))

    # Add plot def template
    paste(s, 'directive plot', sep = '\n')

    # Add the species names to plot
    for(i in 1:length(species_to_plot)) {
        spec <- paste(species_to_plot[[i]], '()', sep = '')
        paste(s, spec, sep = ' ')

        # Check if it is the last species to plot
        if(i != length(species_to_plot)) {
            # If it isn't, add a ';'
            paste(s, ';', sep = '')
        }
    }

    return(s)
}

# TODO: Doc this function
get_dsd_species_str <- function(species_name, domains) {
    # Set the string template
    s <- 'def %name() = (Signal(%domains))'

    # Set the species name
    stringr::str_replace(s, '%name', species_name)

    # Construct the domains string
    domains_str <- paste(domains, collapse = ', ')

    # Set the domains string
    stringr::str_replace(s, '%domains', domains_str)

    return(s)
}

# TODO: Doc this function
get_dsd_4domain_modules_str <- function() {
    s <- '(* Input and output *)
def Signal(unk, i1, i2, i3) = <unk i1^ i2 i3^>

(* Unimolecular to one product reaction *)
def G_1(ia, ib, ic, unko1, o1a) = {ia^*}[ib ic^]<unko1 o1a^>
def T_1(ic, unko, oa, ob, oc) = {ic^*}[unko oa^]<ob oc^>

(* Intermediate states *)
def Waste1_uni(unki, ia, ib, ic) = <unki>[ia^ ib ic^]
def Waste2_1(ib, ic, unko, oa) = <ib>[ic^ unko oa^]
def P_1(ib, ic, unko, oa) = <ib ic^ unko oa^>

def A_e_B(
    qi, qmax, CiA, CiB, Cmax,
    unki, ia, ib, ic,
    unko, oa, ob, oc
) = (
      CiA  * Signal(unki, ia, ib, ic)
    | CiB  * Signal(unko, oa, ob, oc)
    | Cmax * G_1(ia, ib, ic, unko, oa)
    | Cmax * T_1(ic, unko, oa, ob, oc)
    | rxn Signal(unki, ia, ib, ic) + G_1(ia, ib, ic, unko, oa) ->{qi} Waste1_uni(unki, ia, ib, ic) + P_1(ib, ic, unko, oa)
    | rxn P_1(ib, ic, unko, oa) + T_1(ic, unko, oa, ob, oc) ->{qmax} Waste2_1(ib, ic, unko, oa) + Signal(unko, oa, ob, oc)
)

(* Unimolecular to two products reaction *)
def G_2(ia, ib, ic, unko1, o1a, unko2, o2a) = {ia^*}[ib ic^]<unko1 o1a^ unko2 o2a^>
def T_2(ic, unko1, o1a, o1b, o1c, unko2, o2a, o2b, o2c) = {ic^*}[unko1 o1a^]<o1b o1c^>:[unko2 o2a^]<o2b o2c^>

(* Intermediate states *)
def Waste2_2(ib, ic, unko1, o1a, unko2, o2a) = <ib>[ic^ unko1 o1a^ unko2 o2a^]
def P_2(ib, ic, unko1, o1a, unko2, o2a) = <ib ic^ unko1 o1a^ unko2 o2a^>

def A_e_BpC(
    qi, qmax, CiA, CiB, CiC, Cmax,
    unki, ia, ib, ic,
    unko1, o1a, o1b, o1c,
    unko2, o2a, o2b, o2c
) = (
      CiA  * Signal(unki, ia, ib, ic)
    | CiB  * Signal(unko1, o1a, o1b, o1c)
    | CiC  * Signal(unko2, o2a, o2b, o2c)
    | Cmax * G_2(ia, ib, ic, unko1, o1a, unko2, o2a)
    | Cmax * T_2(ic, unko1, o1a, o1b, o1c, unko2, o2a, o2b, o2c)
    | rxn Signal(unki, ia, ib, ic) + G_2(ia, ib, ic, unko1, o1a, unko2, o2a) ->{qi} Waste1_uni(unki, ia, ib, ic) + P_2(ib, ic, unko1, o1a, unko2, o2a)
    | rxn P_2(ib, ic, unko1, o1a, unko2, o2a) + T_2(ic, unko1, o1a, o1b, o1c, unko2, o2a, o2b, o2c) ->{qmax} Waste2_2(ib, ic, unko1, o1a, unko2, o2a) + Signal(unko1, o1a, o1b, o1c) + Signal(unko2, o2a, o2b, o2c)
)

(* Bimolecular to one product reaction *)
def L_1(i1a, i1b, i1c, i2a, i2b, i2c, unko, oa) = {i1a^*}[i1b i1c^ i2a^]:[i2b i2c^]<unko oa^>
def B(i1b, i1c, i2a) = <i1b i1c^ i2a^>

(* Intermediate states *)
def H_1(unki1, i1a, i1b, i1c, i2a, i2b, i2c, unko, oa) = <unki1>[i1a^ i1b i1c^]:{i2a^*}[i2b i2c^]<unko oa^>
def Waste1_bi(unki1, unki2, i1a, i1b, i1c, i2a, i2b, i2c) = <unki1>[i1a^ i1b i1c^]:<unki2>[i2a^ i2b i2c^]

def ApB_e_C(
    qi, qmax, CiA, CiB, CiC, Cmax,
    unka, i1a, i1b, i1c,
    unkb, i2a, i2b, i2c,
    unkc, c1, c2, c3
) = (
      CiA  * Signal(unka, i1a, i1b, i1c)
    | CiB  * Signal(unkb, i2a, i2b, i2c)
    | CiC  * Signal(unkc, c1, c2, c3)
    | Cmax * L_1(i1a, i1b, i1c, i2a, i2b, i2c, unkc, c1)
    | Cmax * B(i1b, i1c, i2a)
    | Cmax * T_1(i2c, unkc, c1, c2 ,c3)
    | rxn Signal(unka, i1a, i1b, i1c) + L_1(i1a, i1b, i1c, i2a, i2b, i2c, unkc, c1) <->{qi}{qmax} H_1(unka, i1a, i1b, i1c, i2a, i2b, i2c, unkc, c1) + B(i1b, i1c, i2a)
    | rxn Signal(unkb, i2a, i2b, i2c) + H_1(unka, i1a, i1b, i1c, i2a, i2b, i2c, unkc, c1) ->{qmax} Waste1_bi(unka, unkb, i1a, i1b, i1c, i2a, i2b, i2c) + P_1(ib, ic, unkc, c1)
    | rxn P_1(ib, ic, unkc, c1) + T_1(i2c, unkc, c1, c2 ,c3) ->{qmax} Waste2_1(i2b, i2c, unkc, c1) + Signal(unkc, c1, c2, c3)
)

(* Bimolecular to two products reaction *)
def L_2(i1a, i1b, i1c, i2a, i2b, i2c, unko1, o1a, unko2, o2a) = {i1a^*}[i1b i1c^ i2a^]:[i2b i2c^]<unko1 o1a^ unko2 o2a^>

(* Intermediate states *)
def H_2(unki1, i1a, i1b, i1c, i2a, i2b, i2c, unko1, o1a, unko2, o2a) = <unki1>[i1a^ i1b i1c^]:{i2a^*}[i2b i2c^]<unko1 o1a^ unko2 o2a^>

def ApB_e_CpD(
    qi, qmax, CiA, CiB, CiC, CiD, Cmax,
    unka, i1a, i1b, i1c,
    unkb, i2a, i2b, i2c,
    unkc, o1a, o1b, o1c,
    unkd, o2a, o2b, o2c
) = (
      CiA  * Signal(unka, i1a, i1b, i1c)
    | CiB  * Signal(unkb, i2a, i2b, i2c)
    | CiC  * Signal(unkc, o1a, o1b, o1c)
    | CiD  * Signal(unkd, o2a, o2b, o2c)
    | Cmax * L_2(i1a, i1b, i1c, i2a, i2b, i2c, unkc, o1a, unkd, o2a)
    | Cmax * B(i1b, i1c, i2a)
    | Cmax * T_2(i2c, unkc, o1a, o1b, o1c, unkd, o2a, o2b, o2c)
    | rxn Signal(unka, i1a, i1b, i1c) + L_2(i1a, i1b, i1c, i2a, i2b, i2c, unkc, o1a, unkd, o2a) <->{qi}{qmax} H_2(unka, i1a, i1b, i1c, i2a, i2b, i2c, unkc, o1a, unkd, o2a) + B(i1b, i1c, i2a)
    | rxn Signal(unkb, i2a, i2b, i2c) + H_2(unka, i1a, i1b, i1c, i2a, i2b, i2c, unkc, o1a, unkd, o2a) ->{qmax} Waste1_bi(unka, unkb, i1a, i1b, i1c, i2a, i2b, i2c) + P_2(i2b, i2c, unkc, o1a, unkd, o2a)
    | rxn P_2(i2b, i2c, unkc, o1a, unkd, o2a) + T_2(i2c, unkc, o1a, o1b, o1c, unkd, o2a, o2b, o2c) ->{qmax} Waste2_2(i2b, i2c, unkc, o1a, unkd, o2a) + Signal(unkc, o1a, o1b, o1c) + Signal(unkd, o2a, o2b, o2c)
)

(* Buffer module *)
def LS(ia, ib, ic, d) = {ia^*}[ib ic^ d^]
def BS(ib, ic, d) = <ib ic^ d^>

(* Intermediate states *)
def HS(unki, ia, ib, ic, d) = <unki>[ia^ ib ic^]{d^*}

def Buff(
    qs, qmax, Cmax, Cii, d,
    unki, ia, ib, ic
) = (
      Cii * Signal(unki, ia, ib, ic)
    | Cmax * LS(ia, ib, ic, d)
    | Cmax * BS(ib, ic, d)
    | rxn Signal(unki, ia, ib, ic) + LS(ia, ib, ic, d) <->{qs}{qmax} HS(unki, ia, ib, ic, d) + BS(ib, ic, d)
)

(* Degradation reactions *)
def A_e_0(qi, CiA, Cmax, unki, ia, ib, ic) = new unko new oa (
      CiA  * Signal(unki, ia, ib, ic)
    | Cmax * G_1(ia, ib, ic, unko, oa)
    | rxn Signal(unki, ia, ib, ic) + G_1(ia, ib, ic, unko, oa) ->{qi} Waste1_uni(unki, ia, ib, ic) + P_1(ib, ic, unko, oa)
)

def ApB_e_0(
    qi, qmax, CiA, CiB, Cmax,
    unka, i1a, i1b, i1c,
    unkb, i2a, i2b, i2c
) = new unko new oa (
      CiA  * Signal(unka, i1a, i1b, i1c)
    | CiB  * Signal(unkb, i2a, i2b, i2c)
    | Cmax * L_1(i1a, i1b, i1c, i2a, i2b, i2c, unko, oa)
    | Cmax * B(i1b, i1c, i2a)
    | rxn Signal(unka, i1a, i1b, i1c) + L_1(i1a, i1b, i1c, i2a, i2b, i2c, unko, oa) <->{qi}{qmax} H_1(unka, i1a, i1b, i1c, i2a, i2b, i2c, unko, oa) + B(i1b, i1c, i2a)
    | rxn Signal(unkb, i2a, i2b, i2c) + H_1(unka, i1a, i1b, i1c, i2a, i2b, i2c, unko, oa) ->{qmax} Waste1_bi(unka, unkb, i1a, i1b, i1c, i2a, i2b, i2c) + P_1(ib, ic, unko, oa)
)'

    return(s)
}

# TODO: Doc this function
get_dsd_AeB_str <- function(
    qi, qmax, CiA, CiB, Cmax,
    unki, ia, ib, ic,
    unko, oa, ob, oc
) {
    # TODO: Implement this function
}

# TODO: Doc this function
get_dsd_AeBpC_str <- function(
    qi, qmax, CiA, CiB, CiC, Cmax,
    unki, ia, ib, ic,
    unko1, o1a, o1b, o1c,
    unko2, o2a, o2b, o2c
) {
    # TODO: Implement this function
}

# TODO: Doc this function
get_dsd_ApBeC_str <- function(
    qi, qmax, CiA, CiB, CiC, Cmax,
    unka, i1a, i1b, i1c,
    unkb, i2a, i2b, i2c,
    unkc, c1, c2, c3
) {
    # TODO: Implement this function
}

# TODO: Doc this function
get_dsd_ApBeCpD_str <- function(
    qi, qmax, CiA, CiB, CiC, CiD, Cmax,
    unka, i1a, i1b, i1c,
    unkb, i2a, i2b, i2c,
    unkc, o1a, o1b, o1c,
    unkd, o2a, o2b, o2c
) {
    # TODO: Implement this function
}

# TODO: Doc this function
get_dsd_buff_str <- function(
    qs, qmax, Cmax, Cii, d,
    unki, ia, ib, ic
) {

}

# TODO: Doc this function
get_dsd_Ae0_str <- function(
    qi, CiA, Cmax, unki, ia, ib, ic
) {
    # TODO: Implement this function
}

# TODO: Doc this function
get_dsd_ApBe0_str <- function(

) {
    # TODO: Implement this function
}

# TODO: Doc this function
get_dsd_def_str <- function(key, val) {
    # String template
    s <- 'def %k = %v'

    # Set the key and value
    s <- stringr::str_replace(s, '%k', key)
    s <- stringr::str_replace(s, '%v', as.character(val))

    return(s)
}

# TODO: Doc this function
save_dsd_script <- function(
    species,
    ci,
    reactions,
    ki,
    qmax,
    cmax,
    alpha,
    beta,
    t,
    filename
) {
    # Check if all reactions attend 4-domain requirements
    for(r in reactions) {
        if(!check_reaction_4domain(r)) {
            stop(paste('Failed to process reaction', r))
        }
    }

    # Initialize new reaction configs
    new_cis <- c(ci * beta)
    new_kis <- c()

    # Get lambda value and buffer k
    buffer_stuff <- get_buff_modules(reactions, ki, qmax, cmax)

    # Change cis according to the lambda^{-1} factor
    if(!is.null(buffer_stuff)) {
        for(i in 1:length(new_cis)) {
            new_cis[i] <- new_cis[[i]] * buffer_stuff$lambda_1
        }
    }

    # Set the script header
    script_str <- paste(get_dsd_header_str(t, species), '\n\n', sep = '')

    # Add 4-domain modules
    script_str <- paste(
        script_str, get_dsd_4domain_modules_str(), '\n\n', sep= ''
    )

    # Set Cmax and qmax
    script_str <- paste(script_str, get_dsd_def_str('Cmax', cmax), sep = '\n')
    script_str <- paste(script_str, get_dsd_def_str('qmax', qmax), sep = '\n')

    # Set the initial concentration of the species
    for(i in 1:length(species)) {
        key <- paste('Ci', species[[i]], sep = '')
        script_str <- paste(
            script_str, get_dsd_def_str(key, new_cis[[i]]), sep = '\n'
        )
    }

    # Set the rate constants
    script_str <- paste(script_str, '\n', sep = '')
    for(i in 1:length(reactions)) {
        if(is_bimolecular(reactions[i])) {
            # Recalculate ki according to the buffer module theory
            k <- ki[i] / (alpha * beta)
            if(!is.null(buffer_stuff)) {
                k <- k * buffer_stuff$lambda_1
            }
            new_kis <- c(new_kis, k)

            # Add k to the script
            script_str <- paste(
                script_str, get_dsd_def_str(
                    paste('k', as.character(i), sep = ''), k
                ), sep = '\n'
            )
        } else {
            # Recalculate ki according to the buffer module theory
            k <- ki[i] / alpha / cmax
            if(!is.null(buffer_stuff)) {
                k <- k * buffer_stuff$lambda_1
            }
            new_kis <- c(new_kis, k)

            # Add k to the script
            script_str <- paste(
                script_str, get_dsd_def_str(
                    paste('k', as.character(i), sep = ''), k
                ), sep = '\n'
            )
        }
    }

    # Set buffer rate constants
    bff_kis <- c()
    # If there is buffer modules to add
    if(!is.null(buffer_stuff)) {
        # Set a counter for naming purposes
        count <- 1

        # Iterate over the rate constant values that aren't cmax
        # (all the even indexes are cmax rate constants)
        for(i in seq(1, length(buffer_stuff$new_ks), by = 2)) {
            # Get the k for the buffer module
            bff_kis <- c(bff_kis, buffer_stuff$new_ks[[i]])

            # Add the k to the script (with the name kb)
            script_str <- paste(
                script_str, get_dsd_def_str(
                    paste('kb', as.character(count), sep = ''),
                    buffer_stuff$new_ks[[i]]
                ), sep = '\n'
            )

            count <- count + 1
        }
    }

    # Set species with their domains
    script_str <- paste(script_str, '\n', sep = '')
    domain_counter <- 1
    species_domains <- list()
    for(spec in species) {
        # Each species has 4 domains
        species_domains[spec] <- c(
            paste('d', as.character(domain_counter)),
            paste('d', as.character(domain_counter + 1)),
            paste('d', as.character(domain_counter + 2)),
            paste('d', as.character(domain_counter + 3))
        )

        # Set the species definition to the script
        script_str <- paste(script_str, get_dsd_species_str(
                spec, species_domains[spec]
            ), sep = '\n'
        )

        # Add the offset to the counter
        domain_counter <- domain_counter + 4
    }

    # TODO: Add reactions to the script
    script_str <- paste(script_str, '\n', sep = '')
}
