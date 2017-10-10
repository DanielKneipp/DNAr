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


#' Calculate the root-mean-square error of two data sets
#'
#' This function is used to calculate the root-mean-square error (RMSE)
#' of two data sets.
#'
#' @param data1  A numerical vector representing one of the data sets.
#' @param data2  A numerical vector representing the other data set.
#'
#' @return  The RMSE value.
#'
#' @export
rmse <- function(data1, data2) {
    return(sqrt(mean((data2 - data1)^2)))
}

#' Calculate the normalized root-mean-square error of two data sets
#'
#' Use this function to calculate the normalized root-mean-square error
#' (NRMSE) of two data sets. In this measures, a distinction between the
#' data sets is needed since the normalization is made using only one of
#' the data sets. For convention, one of the data sets is called simulated
#' data set, and the other one is called observed data set.
#'
#' @param sim_data  The simulated data set
#' @param obs_data  The observed data set (used in the normalization)
#'
#' @return  The NRMSE measure.
#'
#' @export
nrmse <- function(sim_data, obs_data) {
    rmse_num <- rmse(sim_data, obs_data)
    return(rmse_num / (1 + max(obs_data) - min(obs_data)))
}

#' Compare the behavior of two reactions
#'
#' Use this function to compare the behavior of two reactions to see
#' the similarity between them for each species. The normalized
#' root-mean-square error (NRMSE) measure (\code{\link{nrmse}()}) is used
#' to make this comparison.
#'
#' For convention, one of the behaviors is called
#' simulated behavior, and the other one is called observed behavior. This
#' differentiation is important because the order that you pass the
#' behaviors impacts in the result.
#'
#' The normalization made by the NRMSE
#' uses the observed values (`bhv_obs` parameter) only, consequently,
#' `compare_behaviors_nrmse(data1, data2)` results in different measures
#' than `compare_behaviors_nrmse(data2, data1)`.
#'
#' @param bhv_sim  The simulated behavior.
#' @param bhv_obs  The observed behavior (used in the normalization).
#'
#' @return  A data frame with the same columns of the behaviors and one row.
#'          Each value is the NRMSE of that species.
#'
#' @export
compare_behaviors_nrmse <- function(bhv_sim, bhv_obs) {
    # Create an empty data frame with the same columns of the behaviors,
    # but the time column
    result <- data.frame(matrix(nrow = 1, ncol = dim(bhv_sim)[2] - 1))
    column_names <- names(bhv_sim)[2:length(names(bhv_sim))]
    colnames(result) <- column_names

    # Calculate the NRMSE for each column
    for(i in column_names) {
        result[i] <- nrmse(bhv_sim[i], bhv_obs[i])
    }

    # Return the result data frame
    return(result)
}

#' This function returns the concentration derivative of each species
#'
#' This function can be used for study what is impacting each species and
#' how much. this is useful to analyse medium size (dozens of reactions) CRNs.
#' all parameters follows the parameters of \code{\link{react}()}, except
#' the optional `time_point` and `behavior`. If a `time_point` is passed,
#' a `behavior` must be passed as well. If both parameters are set,
#' this functions returns the concentration of each species
#' at a specific point in time within the derivative.
#'
#' @param behavior    The data returned by \code{\link{react}()}.
#' @param time_point  An index (representing a point in time) used to access
#'                    a specific line of `behavior`.
#'
#' @return A data frame with the derivatives. To access the derivative
#'         of a species `'A'`, you just have to access `df['A']`.
#'
#' @export
analyze_behavior <- function(
    species,
    ci,
    reactions,
    ki,
    time_point = NULL,
    behavior = NULL
) {
    # Helper function to concat strings
    jn <- function(...) { paste(..., sep = '') }

    # Check if behavior exists (in case of a time_point has been passed)
    if(!is.null(time_point)) {
        assertthat::assert_that(!is.null(behavior))
    }

    # Get stoichiometry information
    sto_info <- get_M(reactions, species)
    sto_prod <- t(sto_info$prod)
    sto_react <- t(sto_info$react)
    # Get the transpose of the M matrix
    Mt <- t(sto_info$M)

    # Set the output data frame
    df <- data.frame(matrix(nrow = 1, ncol = length(species)))
    names(df) <- species

    for(i in 1:length(species)) {
        # Set the left part of the derivative equation
        s <- jn('d[', species[i], ']/dt = ')
        for(j in 1:length(reactions)) {
            # Set the k with stoichiometry
            k <- ki[j] * Mt[i,j]

            # Go to the next reactions if this one has k = 0
            # (this reaction doesn't impact this species)
            if(k == 0) {
                next
            }

            if(j != 1 && substr(s, nchar(s) - 1, nchar(s)) != '= ') {
                s <- jn(s, ' + ')
            }
            s <- jn(s, '(', k)

            # Get the reactant names or values
            reactants <- get_reactants(reactions[j])
            for(reactant in reactants) {
                reactant_idx <- match(reactant, species)

                # If the reactant is not in the species list
                if(is.na(reactant_idx)) {
                    break
                }

                # Set exponent
                react_exp <- sto_react[reactant_idx, j]

                # Set concentration wit exponent
                if(is.null(time_point)) {
                    s <- jn(s, ' * [', reactant, ']')
                } else {
                    s <- jn(s, ' * ', behavior[time_point, reactant]^react_exp,
                            '[', reactant, ']')
                }
                if(react_exp > 1) {
                    s <- jn(s, '^', react_exp)
                }
            }

            s <- jn(s, ')')
        }

        if(substr(s, nchar(s) - 1, nchar(s)) == '= ') {
            s <- jn(s, '0')
        }

        df[i] <- s
    }

    return(df)
}

