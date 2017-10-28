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


#' Get the DSD script header
#'
#' Use this function to get the header of the DSD script.
#' This header comes with simulation time, ODE solver,
#' the concentration unit measure, type of compilation and
#' the species to be plotted.
#'
#' @param time             The time sequence that simulation will occur.
#'                         this parameter can be exactly the same used in
#'                         the parameter `t` of the \code{\link{react}()}
#'                         or \code{\link{react_4domain}()} functions.
#' @param species_to_plot  A vector with the species names to be plotted.
#'
#' @return  A string with the DSD script header.
get_dsd_header_str <- function(time, species_to_plot) {
    # Set template string
    s <- 'directive duration %t points 1000
directive simulation deterministicstiff
directive concentration M
directive compilation infinite'

    # Set time
    s <- stringr::str_replace(s, '%t', as.character(max(time)))

    # Add plot def template
    s <- paste(s, 'directive plot', sep = '\n')

    # Add the species names to plot
    for(i in 1:length(species_to_plot)) {
        spec <- paste(species_to_plot[[i]], '()', sep = '')
        s <- paste(s, spec, sep = ' ')

        # Check if it is the last species to plot
        if(i != length(species_to_plot)) {
            # If it isn't, add a ';'
            s <- paste(s, ';', sep = '')
        }
    }

    return(s)
}

#' Get the DSD definition of a species
#'
#' Use this function to get a string representing the definition
#' of a species.
#'
#' @param species_name  The species name.
#' @param domains       A vector of four string representing the
#'                      domain names of this species.
#'
#' @return  A string with the species definition.
get_dsd_species_str <- function(species_name, domains) {
    # Set the string template
    s <- 'def %name() = (Signal(%domains))'

    # Set the species name
    s <- stringr::str_replace(s, '%name', species_name)

    # Construct the domains string
    domains_str <- paste(domains, collapse = ', ')

    # Set the domains string
    s <- stringr::str_replace(s, '%domains', domains_str)

    return(s)
}

#' Get the 4-domain modules for the DSD script.
#'
#' Use this function to get a string with all 4-domain
#' modules for the DSD script.
#'
#' @return String with the 4-domain modules.
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
def W(i1b, i1c, i2a) = <i1b i1c^ i2a^>

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
    | Cmax * W(i1b, i1c, i2a)
    | Cmax * T_1(i2c, unkc, c1, c2 ,c3)
    | rxn Signal(unka, i1a, i1b, i1c) + L_1(i1a, i1b, i1c, i2a, i2b, i2c, unkc, c1) <->{qi}{qmax} H_1(unka, i1a, i1b, i1c, i2a, i2b, i2c, unkc, c1) + W(i1b, i1c, i2a)
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
    | Cmax * W(i1b, i1c, i2a)
    | Cmax * T_2(i2c, unkc, o1a, o1b, o1c, unkd, o2a, o2b, o2c)
    | rxn Signal(unka, i1a, i1b, i1c) + L_2(i1a, i1b, i1c, i2a, i2b, i2c, unkc, o1a, unkd, o2a) <->{qi}{qmax} H_2(unka, i1a, i1b, i1c, i2a, i2b, i2c, unkc, o1a, unkd, o2a) + W(i1b, i1c, i2a)
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
    | Cmax * W(i1b, i1c, i2a)
    | rxn Signal(unka, i1a, i1b, i1c) + L_1(i1a, i1b, i1c, i2a, i2b, i2c, unko, oa) <->{qi}{qmax} H_1(unka, i1a, i1b, i1c, i2a, i2b, i2c, unko, oa) + W(i1b, i1c, i2a)
    | rxn Signal(unkb, i2a, i2b, i2c) + H_1(unka, i1a, i1b, i1c, i2a, i2b, i2c, unko, oa) ->{qmax} Waste1_bi(unka, unkb, i1a, i1b, i1c, i2a, i2b, i2c) + P_1(ib, ic, unko, oa)
)

(* Formation reactions *)
(* These reaction get unstable when G_1 or G_2 run out *)
def r0_e_A(
    qi, qmax, CiA, Cmax,
    unko, oa, ob, oc
) = new unki new ia new ib new ic (
      CiA  * Signal(unko, oa, ob, oc)
    | Cmax * G_1(ia, ib, ic, unko, oa)
    | Cmax * T_1(ic, unko, oa, ob, oc)
    | rxn G_1(ia, ib, ic, unko, oa) ->{qi} Waste1_uni(unki, ia, ib, ic) + P_1(ib, ic, unko, oa)
    | rxn P_1(ib, ic, unko, oa) + T_1(ic, unko, oa, ob, oc) ->{qmax} Waste2_1(ib, ic, unko, oa) + Signal(unko, oa, ob, oc)
)

def r0_e_ApB(
    qi, qmax, CiA, CiB, Cmax,
    unko1, o1a, o1b, o1c,
    unko2, o2a, o2b, o2c
) = new unki new ia new ib new ic (
      CiA  * Signal(unko1, o1a, o1b, o1c)
    | CiB  * Signal(unko2, o2a, o2b, o2c)
    | Cmax * G_2(ia, ib, ic, unko1, o1a, unko2, o2a)
    | Cmax * T_2(ic, unko1, o1a, o1b, o1c, unko2, o2a, o2b, o2c)
    | rxn G_2(ia, ib, ic, unko1, o1a, unko2, o2a) ->{qi} Waste1_uni(unki, ia, ib, ic) + P_2(ib, ic, unko1, o1a, unko2, o2a)
    | rxn P_2(ib, ic, unko1, o1a, unko2, o2a) + T_2(ic, unko1, o1a, o1b, o1c, unko2, o2a, o2b, o2c) ->{qmax} Waste2_2(ib, ic, unko1, o1a, unko2, o2a) + Signal(unko1, o1a, o1b, o1c) + Signal(unko2, o2a, o2b, o2c)
)'

    return(s)
}

#' Replaces markers on a template string by other strings
#'
#' This function can be used to replace markers in a string to other
#' strings according to a pattern.
#'
#' @param template_str   The template string (string with the markers).
#' @param marker_pattern  The pattern of the markers. It is a regular
#'                        expression rule.
#' @param objs            A vector with the strings that will replace the
#'                        markers. The replacement is made in order, that is,
#'                        the first markers will be replaced by the first
#'                        strings specified in this vector. This vector
#'                        must have a number of strings less or equal
#'                        the number of markers in the template string.
#'
#' @return  The template string with the markers replaced.
replace_module_markers <- function(template_str, marker_pattern, objs) {
    # Iterate over the objects (variable names), replacing the markers
    # by them.
    for(obj_str in objs) {
        template_str <- stringr::str_replace(
            template_str, marker_pattern, obj_str
        )
    }

    return(template_str)
}

#' Get a DSD 4domain module string
#'
#' Based on a template string, this function prepare it with the
#' the parameters passed and returns a string representing the
#' module with all the parameters set.
#'
#' @param template_str  The template string of the module;
#' @param ...           The parameters to be setted on the module string.
#'
#' @return  A string representing an instantiation of a module with all
#'          parameters set.
dsd_4d_make_module_str <- function(template_str, ...) {
    # Put all variables in one vector
    data <- c(...)

    # Replace the marker by the actual variable names
    s <- replace_module_markers(template_str, '%[a-zA-Z0-9]+', data)

    return(s)
}

#' Get the arguments of the current function
#'
#' Use this function to get the arguments of the current function.
#'
#' @return  A list with the name of the arguments being the name of the
#'          element, and the argument's content being the content of
#'          the element.
arg_list <- function() {
    # Get the parameters of this function
    parms <- as.list(match.call(
        definition = sys.function(sys.parent()),
        call = sys.call(sys.parent()),
        envir = parent.frame(2L)
    ))

    # Remove the first element which is the name of the function it self
    parms <- parms[2:length(parms)]

    # Evaluate languages
    pf <- parent.frame(2L)
    parms <- lapply(parms, function(p) {
        if(is.language(p)) {
            eval(p, envir = pf)
        } else {
            p
        }
    })

    return(parms)
}

#' Instantiate a 'A -> B' module in the DSD script
#'
#' This function returns a string representing an addition of a
#' 'A -> B' reaction in the DSD script. It creates a `A_e_B()` module
#' in the script, replacing all the parameter strings by the ones
#' specified in this function.
#'
#' @param qi         String representing a variable name for the `qi`
#'                   parameter of the module signature.
#' @param qmax       String representing a variable name for the `qmax`
#'                   parameter of the module signature.
#' @param CiA        String representing a variable name for the `CiA`
#'                   parameter of the module signature.
#' @param CiB        String representing a variable name for the `CiB`
#'                   parameter of the module signature.
#' @param Cmax       String representing a variable name for the `Cmax`
#'                   parameter of the module signature.
#' @param A_domains  Vector of strings representing the species domains of A.
#' @param B_domains  Vector of strings representing the species domains of B.
#'
#' @return  A string representing the instantiation of a `A_e_B()` module.
get_dsd_AeB_str <- function(
    qi, qmax, CiA, CiB, Cmax,
    A_domains,
    B_domains
) {
    # String template
    s <- 'A_e_B(
    %qi, %qmax, %CiA, %CiB, %Cmax,
    %unki, %ia, %ib, %ic,
    %unko, %oa, %ob, %oc
)'
    # Get dsd module string
    do.call(dsd_4d_make_module_str, c(s, arg_list()))
}

#' Instantiate a 'A -> B + C' module in the DSD script
#'
#' This function returns a string representing an addition of a
#' 'A -> B + C' reaction in the DSD script. It creates a `A_e_BpC()` module
#' in the script, replacing all the parameter strings by the ones
#' specified in this function.
#'
#' @param qi         String representing a variable name for the `qi`
#'                   parameter of the module signature.
#' @param qmax       String representing a variable name for the `qmax`
#'                   parameter of the module signature.
#' @param CiA        String representing a variable name for the `CiA`
#'                   parameter of the module signature.
#' @param CiB        String representing a variable name for the `CiB`
#'                   parameter of the module signature.
#' @param CiC        String representing a variable name for the `CiC`
#'                   parameter of the module signature.
#' @param Cmax       String representing a variable name for the `Cmax`
#'                   parameter of the module signature.
#' @param A_domains  Vector of strings representing the species domains of A.
#' @param B_domains  Vector of strings representing the species domains of B.
#' @param C_domains  Vector of strings representing the species domains of C.
#'
#' @return  A string representing the instantiation of a `A_e_BpC()` module.
get_dsd_AeBpC_str <- function(
    qi, qmax, CiA, CiB, CiC, Cmax,
    A_domains,
    B_domains,
    C_domains
) {
    # String template
    s <- 'A_e_BpC(
    %qi, %qmax, %CiA, %CiB, %CiC, %Cmax,
    %unki, %ia, %ib, %ic,
    %unko1, %o1a, %o1b, %o1c,
    %unko2, %o2a, %o2b, %o2c
)'

    # Get dsd module string
    do.call(dsd_4d_make_module_str, c(s, arg_list()))
}

#' Instantiate a 'A + B -> C' module in the DSD script
#'
#' This function returns a string representing an addition of a
#' 'A + B -> C' reaction in the DSD script. It creates a `ApB_e_C()` module
#' in the script, replacing all the parameter strings by the ones
#' specified in this function.
#'
#' @param qi         String representing a variable name for the `qi`
#'                   parameter of the module signature.
#' @param qmax       String representing a variable name for the `qmax`
#'                   parameter of the module signature.
#' @param CiA        String representing a variable name for the `CiA`
#'                   parameter of the module signature.
#' @param CiB        String representing a variable name for the `CiB`
#'                   parameter of the module signature.
#' @param CiC        String representing a variable name for the `CiC`
#'                   parameter of the module signature.
#' @param Cmax       String representing a variable name for the `Cmax`
#'                   parameter of the module signature.
#' @param A_domains  Vector of strings representing the species domains of A.
#' @param B_domains  Vector of strings representing the species domains of B.
#' @param C_domains  Vector of strings representing the species domains of C.
#'
#' @return  A string representing the instantiation of a `ApB_e_C()` module.
get_dsd_ApBeC_str <- function(
    qi, qmax, CiA, CiB, CiC, Cmax,
    A_domains,
    B_domains,
    C_domains
) {
    # String template
    s <- 'ApB_e_C(
    %qi, %qmax, %CiA, %CiB, %CiC, %Cmax,
    %unka, %i1a, %i1b, %i1c,
    %unkb, %i2a, %i2b, %i2c,
    %unkc, %c1, %c2, %c3
)'

    # Get dsd module string
    do.call(dsd_4d_make_module_str, c(s, arg_list()))
}

#' Instantiate a 'A + B -> C + D' module in the DSD script
#'
#' This function returns a string representing an addition of a
#' 'A + B -> C + D' reaction in the DSD script. It creates a `ApB_e_CpD()`
#' module in the script, replacing all the parameter strings by the ones
#' specified in this function.
#'
#' @param qi         String representing a variable name for the `qi`
#'                   parameter of the module signature.
#' @param qmax       String representing a variable name for the `qmax`
#'                   parameter of the module signature.
#' @param CiA        String representing a variable name for the `CiA`
#'                   parameter of the module signature.
#' @param CiB        String representing a variable name for the `CiB`
#'                   parameter of the module signature.
#' @param CiC        String representing a variable name for the `CiC`
#'                   parameter of the module signature.
#' @param CiD        String representing a variable name for the `CiD`
#'                   parameter of the module signature.
#' @param Cmax       String representing a variable name for the `Cmax`
#'                   parameter of the module signature.
#' @param A_domains  Vector of strings representing the species domains of A.
#' @param B_domains  Vector of strings representing the species domains of B.
#' @param C_domains  Vector of strings representing the species domains of C.
#' @param D_domains  Vector of strings representing the species domains of D.
#'
#' @return  A string representing the instantiation of a `ApB_e_CpD()` module.
get_dsd_ApBeCpD_str <- function(
    qi, qmax, CiA, CiB, CiC, CiD, Cmax,
    A_domains,
    B_domains,
    C_domains,
    D_domains
) {
    # String template
    s <- 'ApB_e_CpD(
    %qi, %qmax, %CiA, %CiB, %CiC, %CiD, %Cmax,
    %unka, %i1a, %i1b, %i1c,
    %unkb, %i2a, %i2b, %i2c,
    %unkc, %o1a, %o1b, %o1c,
    %unkd, %o2a, %o2b, %o2c
)'

    # Get dsd module string
    do.call(dsd_4d_make_module_str, c(s, arg_list()))
}

#' Instantiate a buffer module in the DSD script
#'
#' This function returns a string representing an addition of
#' a buffer reaction module in the DSD script. It creates a `Buff()`
#' module in the script, replacing all the parameter strings by the ones
#' specified in this function.
#'
#' @param qs         String representing a variable name for the `qs`
#'                   parameter of the module signature.
#' @param qmax       String representing a variable name for the `qmax`
#'                   parameter of the module signature.
#' @param Cii        String representing a variable name for the `Cii`
#'                   parameter of the module signature.
#' @param d          String representing a variable name for the `d`
#'                   parameter of the module signature.
#' @param Cmax       String representing a variable name for the `Cmax`
#'                   parameter of the module signature.
#' @param domains    Vector of strings representing the domains of
#'                   the input species.
#'
#' @return  A string representing the instantiation of a `Buff()` module.
get_dsd_buff_str <- function(
    qs, qmax, Cmax, Cii, d,
    domains
) {
    # Template string
    s <- 'Buff(
    %qs, %qmax, %Cmax, %Cii, %d,
    %unki, %ia, %ib, %ic
)'
    # Get dsd module string
    do.call(dsd_4d_make_module_str, c(s, arg_list()))
}

#' Instantiate a 'A -> 0' module in the DSD script
#'
#' This function returns a string representing an addition of a
#' 'A -> 0' reaction in the DSD script. It creates a `A_e_0()` module
#' in the script, replacing all the parameter strings by the ones
#' specified in this function.
#'
#' @param qi         String representing a variable name for the `qi`
#'                   parameter of the module signature.
#' @param CiA        String representing a variable name for the `CiA`
#'                   parameter of the module signature.
#' @param Cmax       String representing a variable name for the `Cmax`
#'                   parameter of the module signature.
#' @param A_domains  Vector of strings representing the species domains of A.
#'
#' @return  A string representing the instantiation of a `A_e_0()` module.
get_dsd_Ae0_str <- function(
    qi, CiA, Cmax,
    A_domains
) {
    # Template string
    s <- 'A_e_0(
    %qi, %CiA, %Cmax, %unki, %ia, %ib, %ic
)'

    # Get dsd module string
    do.call(dsd_4d_make_module_str, c(s, arg_list()))
}

#' Instantiate a 'A + B -> 0' module in the DSD script
#'
#' This function returns a string representing an addition of a
#' 'A + B -> 0' reaction in the DSD script. It creates a `ApB_e_0()` module
#' in the script, replacing all the parameter strings by the ones
#' specified in this function.
#'
#' @param qi         String representing a variable name for the `qi`
#'                   parameter of the module signature.
#' @param qmax       String representing a variable name for the `qmax`
#'                   parameter of the module signature.
#' @param CiA        String representing a variable name for the `CiA`
#'                   parameter of the module signature.
#' @param CiB        String representing a variable name for the `CiB`
#'                   parameter of the module signature.
#' @param Cmax       String representing a variable name for the `Cmax`
#'                   parameter of the module signature.
#' @param A_domains  Vector of strings representing the species domains of A.
#' @param B_domains  Vector of strings representing the species domains of B.
#'
#' @return  A string representing the instantiation of a `ApB_e_0()` module.
get_dsd_ApBe0_str <- function(
    qi, qmax, CiA, CiB, Cmax,
    A_domains,
    B_domains
) {
    # Template string
    s <- 'ApB_e_0(
    %qi, %qmax, %CiA, %CiB, %Cmax,
    %unka, %i1a, %i1b, %i1c,
    %unkb, %i2a, %i2b, %i2c
)'

    # Get dsd module string
    do.call(dsd_4d_make_module_str, c(s, arg_list()))
}

#' Instantiate a '0 -> A' module in the DSD script
#'
#' This function returns a string representing an addition of a
#' '0 -> A' reaction in the DSD script. It creates a `r0_e_A()` module
#' in the script, replacing all the parameter strings by the ones
#' specified in this function.
#'
#' @param qi         String representing a variable name for the `qi`
#'                   parameter of the module signature.
#' @param qmax       String representing a variable name for the `qmax`
#'                   parameter of the module signature.
#' @param CiA        String representing a variable name for the `CiA`
#'                   parameter of the module signature.
#' @param Cmax       String representing a variable name for the `Cmax`
#'                   parameter of the module signature.
#' @param A_domains  Vector of strings representing the species domains of A.
#'
#' @return  A string representing the instantiation of a `r0_e_A()` module.
get_dsd_0eA_str <- function(
    qi, qmax, CiA, Cmax,
    A_domains
) {
    # Template string
    s <- 'r0_e_A(
    %qi, %qmax, %CiA, %Cmax, %unko, %oa, %ob, %oc
)'

    # Get dsd module string
    do.call(dsd_4d_make_module_str, c(s, arg_list()))
}

#' Instantiate a '0 -> A + B' module in the DSD script
#'
#' This function returns a string representing an addition of a
#' '0 -> A + B' reaction in the DSD script. It creates a `r0_e_ApB()` module
#' in the script, replacing all the parameter strings by the ones
#' specified in this function.
#'
#' @param qi         String representing a variable name for the `qi`
#'                   parameter of the module signature.
#' @param qmax       String representing a variable name for the `qmax`
#'                   parameter of the module signature.
#' @param CiA        String representing a variable name for the `CiA`
#'                   parameter of the module signature.
#' @param CiB        String representing a variable name for the `CiB`
#'                   parameter of the module signature.
#' @param Cmax       String representing a variable name for the `Cmax`
#'                   parameter of the module signature.
#' @param A_domains  Vector of strings representing the species domains of A.
#' @param B_domains  Vector of strings representing the species domains of B.
#'
#' @return  A string representing the instantiation of a `r0_e_ApB()` module.
get_dsd_0eApB_str <- function(
    qi, qmax, CiA, CiB, Cmax,
    A_domains,
    B_domains
) {
    # Template string
    s <- 'r0_e_ApB(
    %qi, %qmax, %CiA, %CiB, %Cmax,
    %unko1, %o1a, %o1b, %o1c,
    %unko2, %o2a, %o2b, %o2c
)'

    # Get dsd module string
    do.call(dsd_4d_make_module_str, c(s, arg_list()))
}

#' Get a definition string for the DSD script
#'
#' This function returns a string representing a definition in the
#' DSD script.
#'
#' A definition is made by the following construction:
#' `def *key* = *value*`. passing a key and value to this function,
#' it returns a definition construction replacing the *key* and *value*
#' by the parameters `key` and `value` specified.
#'
#' @param key  String representing the key of the definition.
#' @param val  String representing the value of the definition.
#'
#' @return The definition string.
get_dsd_def_str <- function(key, val) {
    # String template
    s <- 'def %k = %v'

    # Set the key and value
    s <- stringr::str_replace(s, '%k', key)
    s <- stringr::str_replace(s, '%v', as.character(val))

    return(s)
}

#' Get the available 4-domain modules for DSD script
#'
#' This function returns a matrix with all available 4-domain modules
#' for the DSD script. The matrix is organized based on the reaction
#' stoichiometry.
#'
#' If a reaction has 2 reactants and 1 product,
#' will be in the position 3,2. The formation reactions are in the
#' line 1, and the degradation reactions are in the column 1 ([[,1]]).
#'
#' The position [[1,1]] is filled with `NA` because it represents an
#' invalid reaction (`0 -> 0`).
#'
#' @return  The matrix of functions which can be accesses with [[i,j]]
dsd_4d_modules <- function() {
    modules <- matrix(list(), nrow = 3, ncol = 3)

    # 0 -> 0 is not a valid reaction
    modules[[1,1]] <- NA
    # 0 -> A
    modules[[1,2]] <- get_dsd_0eA_str
    # 0 -> A + B
    modules[[1,3]] <- get_dsd_0eApB_str
    # A -> 0
    modules[[2,1]] <- get_dsd_Ae0_str
    # A -> B
    modules[[2,2]] <- get_dsd_AeB_str
    # A -> B + C
    modules[[2,3]] <- get_dsd_AeBpC_str
    # A + B -> 0
    modules[[3,1]] <- get_dsd_ApBe0_str
    # A + B -> C
    modules[[3,2]] <- get_dsd_ApBeC_str
    # A + B -> C + D
    modules[[3,3]] <- get_dsd_ApBeCpD_str

    return(modules)
}

#' Selects a module based on the reaction stoichiometry
#'
#' This function returns a 4-domain DSD module based on the number
#' of products and reactants (stoichiometry) of the reaction.
#'
#' @param n_reac,n_prod  Number of reactants and products.
#'
#' @return  A list with two parameters:
#'           - `mod_func`: The module function obtained from
#'             `\link{dsd_4d_modules}()`;
#'           - `idx`: A vector with two numbers specifying the position
#'             used to get the module function.
select_dsd_4d_module <- function(n_reac, n_prod) {
    # Return the module function and a reference index
    list(
        mod_func = dsd_4d_modules()[[n_reac + 1, n_prod + 1]],
        idx = c(n_reac + 1, n_prod + 1)
    )
}

#' Binds the module function with its parameters.
#'
#' Use this function to bind a module function to its parameters.
#'
#' @param module           The module function;
#' @param k_idx            The k index to name the k for the module;
#' @param species_domains  The list with the species domains previously
#'                         generated
#' @param reactants        Reactants of the reaction. If a reaction is
#'                         something like `2A -> B`, a list with
#'                         two separated A's `('A', 'A')` must be passed.
#'                         If a reaction is of formation, `NULL` must be passed
#' @param products         The reaction products. It has the same requirements
#'                         of `reactants`.
#'
#' @return  A function with no arguments. When this function is called, a
#'          string representing the modules with its arguments defined will
#'          be returned.
prepare_dsd_4d_module <- function(
    module,
    k_idx,
    species_domains,
    reactants,
    products
) {
    # Get the arguments needed by the module
    args <- names(formals(module))

    # Function to append at the end of a list
    append <- function(l, var, idx = NULL) {
        # If idx is null, set it to the last position of the list
        if(is.null(idx)) {
            idx <- length(l) + 1
        }

        # Add the element
        l[[idx]] <- var

        return(l)
    }

    # Helper to add a parameter to the parms list if the parameter
    # is needed by the module
    add_parm <- function(parms, param_val, param_name = NULL) {
        # Check if the parameter is needed
        if(param_name %in% args) {
            # Add the parameter to the list
            parms <- append(parms, param_val, param_name)
        }

        return(parms)
    }

    # Define other parameters needed by some modules
    other_parms <- list(
        list(name = 'qi', val = paste('k', as.character(k_idx), sep = '')),
        list(name = 'qmax', val = 'qmax'),
        list(name = 'Cmax', val = 'Cmax')
    )

    # Extract the Ci parameters
    ci_parms <- lapply(c(reactants, products), function(species) {
        paste('Ci', species, sep = '')
    })

    # Extract the domain parameters
    domain_parms <- lapply(c(reactants, products), function(species) {
        species_domains[[species]]
    })

    # Define the list of parameters
    parms <- list()
    for (p in other_parms) {
        parms <- add_parm(parms, p$val, p$name)
    }
    parms <-c(parms, ci_parms, domain_parms)

    # Make the module a function without arguments
    function() {
        do.call(module, parms)
    }
}

#' Returns a string with the 4-domain dsd module
#'
#' Use this function to get a string representing a 4-domain
#' dsd module for your reaction.
#'
#' @param reaction         The reaction for which a dsd module will be
#'                         instantiated
#' @param species_domains  A list of lists in which there is a list of domains
#'                         for each species.
#' @param k_idx            The index of k, used to name the k variable.
#'
#' @return  A string representing the dsd module instantiated.
dsd_4d_module_str <- function(reaction, species_domains, k_idx) {
    # Get The number of reactants and products
    stoichiometry <- get_stoichiometry_all(reaction)

    # get_stoichiometry_all consider 0 a species (0 -> A has a left
    # stoichiometry of 1). To select the correct dsd module, we have
    # to change the stoichiometry to 0 in this cases
    if(is_formation(reaction)) {
        stoichiometry$left_sto <- 0
    }

    # Do the same for degradation reactions
    if(is_degradation(reaction)) {
        stoichiometry$right_sto <- 0
    }

    # Get the species names of reactants and products
    # The reaction is of formation or degradation, the reactants
    # or products will be NULL, respectively
    reactants <- NULL
    if(!is_formation(reaction)) {
        reactants <- get_reactants(reaction)
        reactants <- expand_species(
            reactants,
            function(species) {
                get_stoichiometry_onespecies(species, reaction)$left_sto
            }
        )
    }
    products <- NULL
    if(!is_degradation(reaction)) {
        products <- get_products(reaction)
        products <- expand_species(
            products,
            function(species) {
                get_stoichiometry_onespecies(species, reaction)$right_sto
            }
        )
    }

    # Get the correct module
    module <- select_dsd_4d_module(
        stoichiometry$left_sto,
        stoichiometry$right_sto
    )

    # Get the module with all the parameter already set
    mod <- prepare_dsd_4d_module(
        module$mod_func,
        k_idx,
        species_domains,
        reactants,
        products
    )

    # Run the module function which will return the string module
    mod()
}

#' Export a DSD script for a CRN
#'
#' Use this function to export a CRN to a DSD script using the 4-domain
#' approach of implementing CRNs using DNA. The parameters are the same
#' of \code{\link{react_4domain}()}, except for the `filename`, which is
#' the name of the file that will be exported and should be specified with
#' an extension, since there is no standard extension for DSD scripts.
#'
#' @export
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
    # Check the CRN
    reactions <- check_crn(species, ci, reactions, ki, t)

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
            script_str, get_dsd_def_str(
                key, format(new_cis[[i]], nsmall = 1)
            ), sep = '\n'
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

    # Set the variable names of the buffer rate constants
    bff_kis_str <- c()
    # If there is buffer modules to add
    if(!is.null(buffer_stuff)) {
        # Set a counter for naming purposes
        count <- 1

        # Iterate over the rate constant values that aren't cmax
        # (all the even indexes are cmax rate constants)
        for(i in seq(1, length(buffer_stuff$new_ks), by = 2)) {
            # Set buffer k variable name
            bff_k_str <- paste('kb', as.character(count), sep = '')

            # Get the k for the buffer module
            bff_kis_str <- c(bff_kis_str, bff_k_str)

            # Add the k to the script (with the name kb)
            script_str <- paste(
                script_str, get_dsd_def_str(
                    bff_k_str,
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
        species_domains[[spec]] <- c(
            paste('d', as.character(domain_counter), sep = ''),
            paste('d', as.character(domain_counter + 1), sep = ''),
            paste('d', as.character(domain_counter + 2), sep = ''),
            paste('d', as.character(domain_counter + 3), sep = '')
        )

        # Set the species definition to the script
        script_str <- paste(script_str, get_dsd_species_str(
            spec, species_domains[[spec]]
        ), sep = '\n')

        # Add the offset to the counter
        domain_counter <- domain_counter + 4
    }

    # Add reactions into another variable to be added
    # to the script afterwards
    modules_str <- '\n\n(\n\n'

    # Set counters for each input species to count how many times
    # they are used on the modules. This counters will be used
    # to change the initial concentration ('Ci') of each species.
    species_mod_counter <- list()
    for(spec in species) {
        species_mod_counter[[spec]] <- 0
    }

    for(i in 1:length(reactions)) {
        # Get the species names of reactants and products
        reactants <- unique(get_reactants(reactions[[i]]))
        products <- unique(get_products(reactions[[i]]))

        # Update the counters according to the species stoichiometry.
        unique_specs_in_react <- unique(c(reactants, products))
        for(spec in unique_specs_in_react) {
            num <- get_stoichiometry_onespecies(spec, reactions[[i]])
            num <- Reduce('+', num)
            species_mod_counter[[spec]] <- species_mod_counter[[spec]] + num
        }

        # Set module str
        mod_str <- dsd_4d_module_str(reactions[[i]], species_domains, i)

        # Put the current module string together with the other ones
        modules_str <- paste(modules_str, mod_str, sep = '')

        # Prepare the script for new modules, if necessary
        if(i != length(reactions)) {
            modules_str <- paste(modules_str, ' |\n', sep = '')
        }
    }

    # Add buffer reactions (if necessary) to the script
    if(!is.null(buffer_stuff)) {
        # Prepare the script to receive more modules
        modules_str <- paste(modules_str, ' |\n', sep = '')

        # Set independent index for the rate constants
        ii <- 1

        # Iterate over the first reaction of each buffer module
        for(i in seq(1, length(buffer_stuff$new_reactions), by = 2)) {
            # Get the reactants of the first reaction of the buffer module
            reactants <- get_reactants(buffer_stuff$new_reactions[[i]])

            # Get the first reactant, must be one of the input species
            first_reactant <- reactants[[1]]

            # Add one to the counter of this species
            species_mod_counter[[first_reactant]] <-
                species_mod_counter[[first_reactant]] + 1

            # Get the domains of this species
            spec_domains <- species_domains[[first_reactant]]

            # Get buff str module
            bff_mod_str <- get_dsd_buff_str(
                bff_kis_str[[ii]], 'qmax', 'Cmax',
                paste('Ci', first_reactant, sep = ''),
                paste('d', as.character(domain_counter), sep = ''),
                spec_domains
            )

            # Add the buffer module
            modules_str <- paste(modules_str, bff_mod_str, sep = '')

            # Add 1 to the independent counter and the domain counter
            ii <- ii + 1
            domain_counter <- domain_counter + 1

            # If there is more modules prepare the script to it
            if(i != (length(buffer_stuff$new_reactions) - 1)) {
                modules_str <- paste(modules_str, ' |\n', sep = '')
            }
        }
    }

    modules_str <- paste(modules_str, '\n\n)\n', sep = '')

    # Change the species Cis according ot their counters
    script_str <- paste(script_str, '\n\n', sep = '')
    for(spec in species) {
        # Get the variable name of the initial concentration
        var_name <- paste('Ci', spec, sep = '')

        # Redefine the variable in the script, dividing its
        # value by the number specified in the counter
        script_str <- paste(script_str, get_dsd_def_str(
            key = var_name,
            val = paste(
                var_name, '/',
                format(species_mod_counter[[spec]], nsmall = 1),
                sep = ''
            )
        ), sep = '\n')
    }

    # Add the modules into the script
    script_str <- paste(script_str, modules_str, sep = '')

    # Save the script into the file
    cat(script_str, file = filename, sep = '')
}
