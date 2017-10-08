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

    # Create the plot
    g <- ggplot2::ggplot(dfm, ggplot2::aes(
        time, value, color = variable
    )) +
        ggplot2::geom_line(size = 1.3) +
        ggplot2::theme_minimal(base_size = 18) +
        ggplot2::labs(x = x_label, y = y_label, color = legend_name) +
        ggplot2::scale_color_brewer(palette="Dark2")

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
#' @param species    The species of all CRN.
#' @param cis        The initial concentration of the species.
#' @param reactions  The reactions that will be saved.
#' @param kis        The rate constant of those reactions.
#' @param filename   The name of the file (without extension) that will be
#'                   saved. The extension `.txt` will be added automatically.
#'
#' @export
#'
#' @importFrom utils capture.output
save_reactions_txt <- function(species, cis, reactions, kis, filename) {
    # Create two data frames, one with the reactions and rate constants and
    # the other with the species and initial concentrations.
    spec_data <- data.frame(species, cis)
    react_data <- data.frame(reactions, kis)

    # Rename the data frame columns
    names(react_data) <- c('Reaction', 'Rate Constant')
    names(spec_data) <- c('Species', 'Initial Concentration')

    # Get the file name
    filename <- paste(filename, '.txt', sep = '')

    # Store the data frame
    capture.output(
        print(spec_data, row.names = FALSE, print.gap = 4),
        file = filename
    )
    cat('\n', file = filename, append = TRUE)
    capture.output(
        print(react_data, row.names = FALSE, print.gap = 4),
        file = filename,
        append = TRUE
    )
}

#' Save a reaction behavior
#'
#' This function saves a reaction behavior returned by \code{\link{react}()}
#' in a csv file.
#'
#' @param behavior  The reaction behavior returned by \code{\link{react}()}.
#' @param filename  The name of the file that will be saved. You don't need
#'                  to specify the extension of the file, the `.csv` extension
#'                  will be added automatically at the end of the name
#'                  specified.
#'
#' @export
#'
#' @importFrom utils write.csv
save_behavior_csv <- function(behavior, filename) {
    filename <- paste(filename, '.csv', sep = '')
    write.csv(behavior, file = filename, row.names = FALSE)
}

#' Load a reaction behavior
#'
#' This function loads a reaction behavior saved by
#' \code{\link{save_behavior_csv}()}.
#'
#' @param filename  The name of the file that will be loaded. You don't need
#'                  to specify the extension of the file, the `.csv` extension
#'                  will be added automatically at the end of the name
#'                  specified.
#'
#' @return  A data frame with the reaction behavior.
#'
#' @export
#'
#' @importFrom utils read.csv
load_behavior_csv <- function(filename) {
    filename <- paste(filename, '.csv', sep = '')
    data <- read.csv(filename, header = TRUE)
    return(data)
}
