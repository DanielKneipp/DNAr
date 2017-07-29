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
    m <- gregexpr(paste('[1-9]*', one_species, sep = ''), reaction_part)
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

#' Simulate a CRN
#'
#' This is the function used to actually simulate the chemical reaction network.
#' Given the CRN specifications, it returns the behavior of the reaction.
#'
#' @section Known limitations:
#'   \itemize{
#'     \item It only supports uni or bimolecular reactions;
#'     \item Bidirectional reactions (e.g.: 'A <--> B') are not supported
#'   yet (you have to describe them with two separated reactions).
#'   }
#'
#' @param species    A vector with the species of the reaction. The order of
#'                   this vector is important because it will define the
#'                   column order of the returned behavior.
#' @param ci         A vector specifying the initial concentrations of the
#'                   \code{species} specified, in order.
#' @param reactions  A vector with the reactions of the CRN.
#' @param ki         A vector defining the constant rate of each reaction
#'                   in \code{reactions}, in order.
#' @param t          A vector specifying the time interval. Each value
#'                   would be a specific time point.
#'
#' @return A data frame with each line being a specific point in the time
#'         and each column but the first being the concentration of a
#'         species. The first column is the time interval. The column names
#'         are filled with the species's names.
#'
#' @export
#'
#' @example demo/main_crn.R
react <- function(species, ci, reactions, ki, t) {
    products <- matrix(data = 0,
                       nrow = length(reactions),
                       ncol = length(species))
    reactants <- matrix(data = 0,
                        nrow = length(reactions),
                        ncol = length(species))

    for(i in 1:length(reactions)) {
        for(j in 1:length(species)) {
            stoc <- get_stoichiometry_onespecies(species[j], reactions[i])
            products[i,j] <- stoc$right_sto
            reactants[i,j] <- stoc$left_sto
        }
    }

    M <- products - reactants
    Mt <- t(M)

    fx <- function(t, y, parms) {
        dy <- numeric(length(y))

        v <- matrix(data = 0, nrow = length(reactions), ncol = 1)
        for(i in 1:length(reactions)) {
            s_is_reac <- reactants_in_reaction(species, reactions[i])
            v[i,1] <- ki[i] * prod(y[s_is_reac])
        }

        dy <- Mt %*% v

        return(list(dy))
    }

    result <- deSolve::ode(times = t, y = ci, func = fx, parms = NULL)

    # Convert double matrix to dataframe
    df_result <- data.frame(result)
    names(df_result) <- c('time', species)

    return(df_result)
}

#' `ggplot2` theme developed by Max Woolf.
#'
#' This is a simple theme for `ggplot2` package designed by Max Woolf.
#'
#' @references
#'   - [Website with the theme.](http://minimaxir.com/2015/02/ggplot-tutorial/)
fte_theme <- function() {
    # Generate the colors for the chart procedurally with RColorBrewer
    palette <- RColorBrewer::brewer.pal("Greys", n=9)
    color.background = palette[2]
    color.grid.major = palette[3]
    color.axis.text = palette[6]
    color.axis.title = palette[7]
    color.title = palette[9]

    # Begin construction of chart
    ggplot2::theme_bw(base_size=9) +

        # Set the entire chart region to a light gray color
        ggplot2::theme(panel.background = ggplot2::element_rect(
            fill = color.background, color = color.background
        )) +
        ggplot2::theme(plot.background = ggplot2::element_rect(
            fill = color.background, color = color.background
        )) +
        ggplot2::theme(panel.border = ggplot2::element_rect(
            color = color.background
        )) +

        # Format the grid
        ggplot2::theme(panel.grid.major = ggplot2::element_line(
            color = color.grid.major, size = .25
        )) +
        ggplot2::theme(panel.grid.minor = ggplot2::element_blank()) +
        ggplot2::theme(axis.ticks = ggplot2::element_blank()) +

        # Format the legend, but hide by default
        ggplot2::theme(legend.position="none") +
        ggplot2::theme(legend.background = ggplot2::element_rect(
            fill = color.background
        )) +
        ggplot2::theme(legend.text = ggplot2::element_text(
            size = 7,color=color.axis.title
        )) +

        # Set title and axis labels, and format these and tick marks
        ggplot2::theme(plot.title = ggplot2::element_text(
            color = color.title, size = 10, vjust = 1.25
        )) +
        ggplot2::theme(axis.text.x = ggplot2::element_text(
            size = 7, color = color.axis.text
        )) +
        ggplot2::theme(axis.text.y = ggplot2::element_text(
            size = 7, color = color.axis.text
        )) +
        ggplot2::theme(axis.title.x = ggplot2::element_text(
            size = 8, color = color.axis.title, vjust=0
        )) +
        ggplot2::theme(axis.title.y = ggplot2::element_text(
            size = 8, color = color.axis.title, vjust=1.25
        )) +

        # Plot margins
        ggplot2::theme(plot.margin = ggplot2::unit(
            c(0.35, 0.2, 0.3, 0.35), "cm")
        )
}

#' Plot the behavior of a CRN.
#'
#' This function plots the behavior returned by \code{\link{react}()}.
#'
#' @param behavior          Behavior returned by \code{\link{react}()}
#' @param species           The vector with the species that should be plotted.
#'                          If no one is specified, than all of them are
#'                          plotted.
#' @param x_label,y_label   Label name of the x and y axis, respectively.
#'                          The default are 'Time' and 'Concentration'.
#' @param legend_name       Name of the legend. The default is 'Species'.
#' @param save_file_name    Name file that the plot should be saved in.
#'                          Currently files with .pdf and .png were tested but
#'                          it should support any extension supported by
#'                          \code{\link[ggplot2]{ggsave}()}. If no file name
#'                          is specified, the won't be saved.
#'
#' @return  The object returned by \code{\link[ggplot2]{ggplot}()}, so you
#'          can modify the plot or save it in a different way.
#'
#' @export
plot_behavior <- function(
    behavior,
    species = NULL,
    x_label = 'Time',
    y_label = 'Concentration',
    legend_name = 'Species',
    save_file_name = NULL
) {
    # If no species was specified, pick all of them.
    if(is.null(species)) {
        species <- names(behavior)[2:dim(behavior)[2]]
    }

    # Convert the data frame to the proper format
    df <- behavior[,c('time', species)]
    dfm <- reshape2::melt(df, id.vars = 'time')

    # Show the plot
    print(
        g <- ggplot2::ggplot(dfm, ggplot2::aes(
            time, value, color = variable
        )) +
        ggplot2::geom_line(size = 1.3) +
        ggplot2::theme_minimal(base_size = 18) +
        ggplot2::labs(x = x_label, y = y_label, color = legend_name) +
        ggplot2::scale_color_brewer(palette="Dark2")
    )

    # if a name file was specified, save the plot there.
    if(!is.null(save_file_name)) {
        ggplot2::ggsave(save_file_name, dpi = 300, width = 6, height = 4.5)
    }

    # Return the plot object
    return(g)
}

#' Save the reactions in a formatted text file
#'
#' This function can be used for saving the reactions
#' with their rate constants in a formatted text file.
#'
#' @param reactions  The reactions that will be saved.
#' @param kis        The rate constant of those reactions.
#' @param filename   The name of the file (without extension) that will be
#'                   saved. The extension `.txt` will be added automatically.
#'
#' @export
save_reactions_txt <- function(reactions, kis, filename) {
    # Create a data frame with the reactions and rate constants
    data <- data.frame(reactions, kis)
    # Rename the data frame columns
    names(data) <- c('Reaction', 'Rate Constant')
    # Get the file name
    filename <- paste(filename, '.txt', sep = '')
    # Store the data frame
    capture.output(
        print(data, row.names = FALSE, print.gap = 4),
        file = filename
    )
}
