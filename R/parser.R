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


#' Get the first part of a reaction.
#'
#' Given a reaction like 'A + B -> C', this function
#' returns 'A + B '.
get_first_part <- function(react_str) {
    return(sub('->.*', '', react_str))
}

#' Get the second part of a reaction.
#'
#' Given a reaction like 'A + B -> C', this function
#' returns ' C'.
get_second_part <- function(react_str) {
    return(sub('.*->', '', react_str))
}

#' Check if a reactions bimolecular.
#'
#' Give a reaction, this functions checks if it is
#' a bimolecular reaction.
#'
#' @examples
#' DNAr:::is_bimolecular('2A -> B')     # Should return TRUE
#' DNAr:::is_bimolecular('A + B -> C')  # Should return TRUE
#' DNAr:::is_bimolecular('A -> B')      # Should return FALSE
is_bimolecular <- function(react_str) {
    first_part <- get_first_part(react_str)
    return(get_stoichiometry_part(first_part) == 2)
}

#' Check if part of a reaction is equal to 0.
#'
#' This function is useful when you have a reaction and
#' you want to check if it is something like 'A -> 0' (something to nothing).
#' To do this you get the second part of the reaction with
#' \code{\link{get_second_part}()} and then use this function with the second
#' part as the parameter.
#'
#' @examples
#' sp <- DNAr:::get_second_part('A -> 0')
#' DNAr:::isempty_part(sp)  # Should return TRUE
isempty_part <- function(react_part) {
    return(!is.na(suppressWarnings(as.numeric(react_part))) &&
               as.numeric(react_part) == 0)
}

#' Get the count of occurrences of a given species in a reaction part.
#'
#' With the reaction part 'A + B ', for instance, and \code{one_species}
#' begin equal to 'A', this function would return 1. On the case of '2A ' it
#' would return 2.
get_onespecies_count <- function(one_species, reaction_part) {
    m <- gregexpr(paste('\\b[1-9]*', one_species, '\\b', sep = ''), reaction_part)
    matches <- regmatches(reaction_part, m)
    nums <- array(0, length(matches[[1]]))
    if(length(nums) > 0) {
        for(i in 1:length(nums)){
            nums[i] <- as.numeric(sub(one_species, '', matches[[1]][i]))
            if(is.na(nums[i])) {
                nums[i] <- 1
            }
        }
    }
    return(nums)
}

#' Get the stoichiometry of a specific species in a reaction.
#'
#' This function uses \code{\link{get_onespecies_count}()} in the left
#' and right part of a reaction to get the stoichiometry of a species
#' in a reaction.
#'
#' @return A list with \code{left_sto} being the stoichiometry of a
#' species in the left side of a reaction, and \code{right_sto} being
#' the same but for the right side of the reaction.
#'
#' @examples
#' # It should return list(left_sto = 1, right_sto = 2)
#' DNAr:::get_stoichiometry_onespecies('A', 'A + B -> 2A')
#' # It should return list(left_sto = 0, right_sto = 1)
#' DNAr:::get_stoichiometry_onespecies('A', 'B -> A + B')
get_stoichiometry_onespecies <- function(one_species, reaction) {
    f_p <- get_first_part(reaction)
    s_p <- get_second_part(reaction)

    f_p_n <- get_onespecies_count(one_species, f_p)
    s_p_n <- get_onespecies_count(one_species, s_p)

    r <- list(left_sto = sum(f_p_n), right_sto = sum(s_p_n))
    return(r)
}

#' It returns the stoichiometry of part of a reaction.
#'
#' This function is used when you want to know the count of molecules
#' in part of a reaction.
#'
#' @examples
#' DNAr:::get_stoichiometry_part('A + B ')  # It should return 2
#' DNAr:::get_stoichiometry_part('2A ')     # It should return 2
#' DNAr:::get_stoichiometry_part(' C')      # It should return 1
get_stoichiometry_part <- function(reaction_part) {
    matches <- stringr::str_match_all(reaction_part,
                                      '([1-9])*([a-zA-Z0-9_]+)')[[1]]
    count <- 0
    for(i in 1:dim(matches)[1]) {
        n <- as.numeric(matches[i,2])
        if(is.na(n)) {
            n <- 1
        }
        count <- count + n
    }
    return(count)
}

#' Get the stoichiometry of a reaction
#'
#' Use this function to get the stoichiometry of the left and
#' right part of a reaction.
#'
#' @return a list with left_sto and right_sto being the stoichiometry
#' of the left and right part respectively.
#'
#' @examples
#' # Returns list(left_sto = 2, right_sto = 1)
#' DNAr:::get_stoichiometry_all('A + B -> C')
#' # Returns list(left_sto = 1, right_sto = 2)
#' DNAr:::get_stoichiometry_all('A -> 2B')
get_stoichiometry_all <- function(reaction) {
    f_p <- get_first_part(reaction)
    s_p <- get_second_part(reaction)
    r <- list(left_sto = get_stoichiometry_part(f_p),
              right_sto = get_stoichiometry_part(s_p))
    return(r)
}

#' Remove the stoichiometry of a species string.
#'
#' Use this function to get rid of the stoichiometry of a
#' species string.
#'
#' @return The species string without the stoichiometry (number)
#' specifying the number of molecules
#'
#' @examples
#' DNAr:::remove_stoichiometry('2A')   # Returns 'A'
#' DNAr:::remove_stoichiometry('2A2')  # Returns 'A2', The other numbers
#'                                     # are considered part of the species name
remove_stoichiometry <- function(species) {
    no_sto_spec <- c()
    for(i in 1:length(species)) {
        s <- sub('^[1-9]*', '', species[i])
        no_sto_spec <- c(no_sto_spec, s)
    }
    return(no_sto_spec)
}

#' Get the species of a reaction part
#'
#' Given part of a reaction, this function returns
#' the species of it without the stoichiometry.
#'
#' @examples
#' DNAr:::get_species('A + 2B')  # Should return c('A', 'B')
get_species <- function(reaction_part) {
    specs <- strsplit(reaction_part, '[^a-zA-Z0-9_]')[[1]]
    specs <- specs[specs != '']
    specs <- remove_stoichiometry(specs)
    return(specs)
}

#' Get the reactants of a given reaction
#'
#' This function returns the reactants of a reactions,
#' removing their stoichiometry.
#'
#' @examples
#' DNAr:::get_reactants('A + B -> C')  # Returns c('A', 'B')
#' DNAr:::get_reactants('2A -> B')     # Returns c('A')
get_reactants <- function(reaction) {
    f_p <- get_first_part(reaction)
    reactants <- get_species(f_p)
    return(reactants)
}

#' Get the products of a given reaction
#'
#' This function returns the products of a reactions,
#' removing their stoichiometry. If the reaction is of the type
#' 'A -> 0', an empty string '' is returned.
#'
#' @examples
#' DNAr:::get_products('A + B -> C')   # Returns c('C')
#' DNAr:::get_products('2A -> B + C')  # Returns c('B', 'C')
#' DNAr:::get_products('A -> 0')       # Returns c('')
get_products <- function(reaction) {
    s_p <- get_second_part(reaction)
    if(isempty_part(s_p)) {
        return(c(''))
    }
    products <- get_species(s_p)
    return(products)
}

#' Check which species are reactants in a given reaction
#'
#' It is used to check which of the given species are
#' reactants in a reaction, returning a vector with the
#' species indexes in \code{species} that are reactants
#' in \code{reaction}.
#'
#' @return A vector filled with indexes specifying the
#' species that are in a reaction as a reactant.
#'
#' @examples
#' # Should return c(1, 2)
#' DNAr:::reactants_in_reaction(c('A', 'B', 'C'), 'A + B -> C')
#' # Should return c(1)
#' DNAr:::reactants_in_reaction(c('A', 'B', 'C'), '2A -> B + C')
#' # Should return c(1, 3)
#' DNAr:::reactants_in_reaction(c('A', 'B', 'C'), 'A + C -> B')
reactants_in_reaction <- function(species, reaction) {
    words <- get_reactants(reaction)
    r <- c()
    for(i in 1:length(species)) {
        if(any(words == species[i])) {
            r <- c(r, i)
        }
    }
    return(r)
}
