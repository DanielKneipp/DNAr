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
    return(rmse_num / (max(obs_data) - min(obs_data)))
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
#' @param bhv_sim             The simulated behavior.
#' @param bhv_obs             The observed behavior (used in the normalization).
#' @param ignore_time_column  Ignore the time column, the first column of
#'                            a behavior
#'
#' @return  A data frame with the same columns of the behaviors and one row.
#'          Each value is the NRMSE of that species.
#'
#' @export
compare_behaviors_nrmse <- function(bhv_sim, bhv_obs, ignore_time_column = T) {
    column_names <- names(bhv_sim)

    if(ignore_time_column) {
        column_names <- column_names[2:length(column_names)]
    }

    # Create an empty data frame with the same columns of the behaviors,
    # but the time column
    result <- data.frame(matrix(nrow = 1, ncol = length(column_names)))

    colnames(result) <- column_names

    # Calculate the NRMSE for each column
    for(i in column_names) {
        result[i] <- nrmse(bhv_sim[[i]], bhv_obs[[i]])
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
#' @param behavior     The data returned by \code{\link{react}()}.
#' @param time_points  A vector of indexes (representing multiple points in
#'                     time) used for access lines of `behavior`.
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
    time_points = NULL,
    behavior = NULL
) {
    # Helper function to concat strings
    jn <- function(...) { paste(..., sep = '') }

    # Check if behavior exists (in case of a time_point has been passed)
    if(!is.null(time_points)) {
        assertthat::assert_that(!is.null(behavior))

        # Check if all time_points are within the behavior data
        assertthat::assert_that(
            max(time_points) <= dim(behavior)[[1]],
            msg = 'All the time points must be within the behavior data.'
        )
    }

    # Get stoichiometry information
    sto_info <- get_M(reactions, species)
    sto_prod <- t(sto_info$prod)
    sto_react <- t(sto_info$react)
    # Get the transpose of the M matrix
    Mt <- t(sto_info$M)

    # Get a list of data frames (one for each time point)
    # If no time point was passed, only one data frame should
    # be instantiated
    n_df <- 1
    if(!is.null(time_points)) {
        n_df <- length(time_points)
    }
    df_list <- lapply(rep(1, n_df), function(nothing) {
        df <- data.frame(matrix(nrow = 1, ncol = length(species)))
        names(df) <- species
        return(df)
    })

    for(t in (if(is.null(time_points)) 1:1 else 1:length(time_points))) {
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
                    if(is.null(time_points)) {
                        s <- jn(s, ' * [', reactant, ']')
                    } else {
                        s <- jn(s, ' * ',
                                behavior[time_points[[t]], reactant]^react_exp,
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

            df_list[[t]][i] <- s
        }
    }

    # Return the data frame if there is only oen time point in the list
    if(length(df_list) == 1) {
        return(df_list[[1]])
    } else {
        return(df_list)
    }
}

#' Evaluate a derivative returned by
#' \code{\link{analyze_behavior}()}
#'
#' If \code{\link{analyze_behavior}()} was used with a `behavior` and
#' `time_point`, you can use this function to calculate the result of
#' the derivative.
#'
#' @param derivative  Derivative with concentration values returned by
#'                    \code{\link{analyze_behavior}()}.
#'
#' @return A numeric value representing the result of the derivative.
#'
#' @export
eval_derivative <- function(derivative) {
    # Get the part after the '='
    right_part <- stringr::str_split(derivative, '=')[[1]][2]

    # Remove the concentration names (with exponent, if they have)
    express <- stringr::str_replace_all(
        right_part,
        '\\[[^\\[\\]]*\\](\\^[0-9]+)?',
        ''
    )

    # Evaluate expression
    return(eval(parse(text = express)))
}

#' Evaluate subexpressions of a derivative returned by
#' \code{\link{analyze_behavior}()}
#'
#' This functions works like `\link{eval_derivative}()`, but instead of
#' evaluating the entire derivative, it will evaluate parts of it
#' (delimited by `()`).
#'
#' @param derivative  Derivative with concentration values returned by
#'                    \code{\link{analyze_behavior}()}.
#'
#' @return A numeric named vector value representing the result of each part
#'         of the the derivative. The names of the results are the evaluated
#'         subexpressions.
#'
#' @export
eval_derivative_part <- function(derivative) {
    # Get the part after the '='
    right_part <- stringr::str_split(derivative, '=')[[1]][2]

    # Calculating the result of each subexpression (within `()``)
    exps <- stringr::str_match_all(right_part, '\\(.+?\\)')[[1]]
    exp_results <- sapply(exps, function(exp) {
        # Remove the concentration names (with exponent, if they have)
        exp <- stringr::str_replace_all(
            exp,
            '\\[[^\\[\\]]*\\](\\^[0-9]+)?',
            ''
        )

        # Evaluate the subexpression
        eval(parse(text = exp))
    })

    # Return a vector with the results
    return(exp_results)
}
