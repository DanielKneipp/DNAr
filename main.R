library(deSolve)

get_first_part <- function(react_str) {
    return(sub('->.*', '', react_str))
}

get_second_part <- function(react_str) {
    return(sub('.*->', '', react_str))
}

is_bimolecular <- function(react_str) {
    first_part <- get_first_part(react_str)
    return(grepl('\\+', first_part))
}

get_stoichiometry <- function(species, reaction_part) {
    m <- gregexpr(paste('[1-9]*', species, sep = ''), reaction_part)
    matches <- regmatches(reaction_part, m)
    nums <- array(0, length(matches[[1]]))
    if(length(nums) > 0) {
        for(i in 1:length(nums)){
            nums[i] <- as.numeric(sub(species, '', matches[[1]][i]))
            if(is.na(nums[i])) {
                nums[i] <- 1
            }
        }
    }
    return(nums)
}

get_species_flow<- function(species, reactions) {
    flow <- c()
    for(i in 1:length(reactions)) {
        f_p <- get_first_part(reactions[i])
        s_p <- get_second_part(reactions[i])
        
        f_p_n <- get_stoichiometry(species, f_p)
        s_p_n <- get_stoichiometry(species, s_p)
        
        r <- list(react = i, left_sto = f_p_n, right_sto = s_p_n)
        flow <- c(flow, list(r))
    }
    return(flow)
}

# Reaction must be unimolecular or bimolecular
#  E.g.: A + B -> C
#        2A -> B
react <- function(species, ci, reactions, qi, qmax, cmax) {
    n_reactions <- length(reactions)
    
    # Get reaction type, number of auxuliary reactions and Js 
    # (for the diferential equations computation)
    n_aux_spec <- 0
    n_js <- 0
    reaction_type <- numeric(n_reactions)
    for(i in 1:n_reactions) {
        if(is_bimolecular(reactions[i])) {
            reaction_type[i] = 2
            n_aux_spec <- n_aux_spec + 5
            n_js <- n_js + 3
        } else {
            reaction_type[i] = 1
            n_aux_spec <- n_aux_spec + 3
            n_js <- n_js + 2
        }
    }
    
    y <- numeric(length(species) + n_aux_spec)
    Js <- numeric(n_js)
    
    y_js_map <- array(list(), length(y))
    
    for(i in 1:length(species)) {
        flow <- get_species_flow(species[i], reactions)
        y_js_map[i] <- flow
    }
    
    for(i in length(species):length(y)) {
    }
}

behaviour <- react(
    species   = c('A','B','C'),
    ci        = c(1e3, 1e3, 0), 
    reactions = c('A + B -> C'),
    qi        = c(1e-7),
    qmax      = c(1e-3),
    cmax      = c(1e5)
)
